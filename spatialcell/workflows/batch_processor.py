#!/usr/bin/env python3
"""
Batch Processing Module for Spatialcell Pipeline

This module enables efficient batch processing of multiple samples through the complete
Spatialcell workflow, with parallel execution, progress tracking, and comprehensive
error handling and reporting.

Author: Xinyan
License: MIT
"""

import os
import sys
import argparse
import logging
import multiprocessing as mp
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import json
import yaml
from datetime import datetime
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

# Add module paths
sys.path.append(str(Path(__file__).parent.parent))

from .pipeline_config import PipelineConfig, validate_pipeline_inputs
from .complete_pipeline import run_complete_pipeline


def setup_batch_logging(output_dir: str) -> None:
    """
    Setup logging for batch processing.
    
    Args:
        output_dir (str): Base output directory for logs
    """
    log_dir = Path(output_dir) / "batch_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    
    log_file = log_dir / f"batch_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ],
        force=True
    )
    
    logging.info(f"Batch processing logging initialized")
    logging.info(f"Log file: {log_file}")


class BatchConfig:
    """
    Configuration manager for batch processing operations.
    """
    
    def __init__(self,
                 base_output_dir: str,
                 max_parallel: int = 4,
                 continue_on_error: bool = True,
                 generate_summary: bool = True):
        """
        Initialize batch configuration.
        
        Args:
            base_output_dir (str): Base directory for all batch outputs
            max_parallel (int): Maximum number of parallel processes
            continue_on_error (bool): Continue processing other samples if one fails
            generate_summary (bool): Generate comprehensive batch summary
        """
        self.base_output_dir = Path(base_output_dir)
        self.max_parallel = max_parallel
        self.continue_on_error = continue_on_error
        self.generate_summary = generate_summary
        
        # Create base output directory
        self.base_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup batch-specific directories
        self.batch_logs_dir = self.base_output_dir / "batch_logs"
        self.batch_reports_dir = self.base_output_dir / "batch_reports"
        self.failed_configs_dir = self.base_output_dir / "failed_configs"
        
        for dir_path in [self.batch_logs_dir, self.batch_reports_dir, self.failed_configs_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
    
    def get_sample_output_dir(self, sample_name: str) -> Path:
        """Get output directory for specific sample."""
        return self.base_output_dir / "samples" / sample_name
    
    def save_config(self, config_path: str) -> None:
        """Save batch configuration to file."""
        config_data = {
            'base_output_dir': str(self.base_output_dir),
            'max_parallel': self.max_parallel,
            'continue_on_error': self.continue_on_error,
            'generate_summary': self.generate_summary,
            'created_at': datetime.now().isoformat()
        }
        
        with open(config_path, 'w') as f:
            yaml.dump(config_data, f, default_flow_style=False)


def load_sample_configs(config_dir: str) -> List[Tuple[str, PipelineConfig]]:
    """
    Load all sample configuration files from directory.
    
    Args:
        config_dir (str): Directory containing sample configuration files
        
    Returns:
        list: List of (sample_name, config) tuples
        
    Raises:
        ValueError: If no valid configuration files found
    """
    config_dir = Path(config_dir)
    
    if not config_dir.exists():
        raise ValueError(f"Configuration directory not found: {config_dir}")
    
    configs = []
    
    # Look for YAML and JSON config files
    for config_file in config_dir.glob("*.yaml"):
        try:
            config = PipelineConfig.from_file(config_file)
            sample_name = config.sample_name or config_file.stem
            configs.append((sample_name, config))
            logging.info(f"Loaded config for sample: {sample_name}")
        except Exception as e:
            logging.warning(f"Failed to load config {config_file}: {e}")
    
    for config_file in config_dir.glob("*.yml"):
        try:
            config = PipelineConfig.from_file(config_file)
            sample_name = config.sample_name or config_file.stem
            configs.append((sample_name, config))
            logging.info(f"Loaded config for sample: {sample_name}")
        except Exception as e:
            logging.warning(f"Failed to load config {config_file}: {e}")
    
    for config_file in config_dir.glob("*.json"):
        try:
            config = PipelineConfig.from_file(config_file)
            sample_name = config.sample_name or config_file.stem
            configs.append((sample_name, config))
            logging.info(f"Loaded config for sample: {sample_name}")
        except Exception as e:
            logging.warning(f"Failed to load config {config_file}: {e}")
    
    if not configs:
        raise ValueError(f"No valid configuration files found in {config_dir}")
    
    logging.info(f"Loaded {len(configs)} sample configurations")
    return configs


def create_batch_configs_from_template(template_config: str,
                                     sample_info: List[Dict[str, Any]],
                                     output_dir: str) -> List[str]:
    """
    Create batch configuration files from template and sample information.
    
    Args:
        template_config (str): Path to template configuration file
        sample_info (list): List of sample-specific parameter dictionaries
        output_dir (str): Output directory for generated configs
        
    Returns:
        list: List of paths to generated configuration files
    """
    template = PipelineConfig.from_file(template_config)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    generated_configs = []
    
    for sample_params in sample_info:
        # Create copy of template
        sample_config = PipelineConfig.from_dict(template.to_dict())
        
        # Update with sample-specific parameters
        for key, value in sample_params.items():
            if hasattr(sample_config, key):
                setattr(sample_config, key, value)
        
        # Save sample configuration
        config_path = output_dir / f"{sample_config.sample_name}_config.yaml"
        sample_config.save(config_path, format='yaml')
        generated_configs.append(str(config_path))
        
        logging.info(f"Generated config for sample: {sample_config.sample_name}")
    
    return generated_configs


def process_single_sample(args: Tuple[str, PipelineConfig, BatchConfig]) -> Dict[str, Any]:
    """
    Process a single sample through the complete pipeline.
    
    Args:
        args (tuple): (sample_name, config, batch_config) tuple
        
    Returns:
        dict: Processing result summary
    """
    sample_name, config, batch_config = args
    
    start_time = time.time()
    result = {
        'sample_name': sample_name,
        'start_time': datetime.now().isoformat(),
        'success': False,
        'error_message': None,
        'output_dir': None,
        'processing_time': 0,
        'pipeline_steps': []
    }
    
    try:
        # Update config with batch-specific output directory
        sample_output_dir = batch_config.get_sample_output_dir(sample_name)
        config.output_dir = str(sample_output_dir)
        
        # Save sample-specific config
        config_file = sample_output_dir / f"{sample_name}_config.yaml"
        config_file.parent.mkdir(parents=True, exist_ok=True)
        config.save(config_file)
        
        # Run complete pipeline
        logging.info(f"Starting pipeline for sample: {sample_name}")
        run_complete_pipeline(str(config_file))
        
        result['success'] = True
        result['output_dir'] = str(sample_output_dir)
        result['pipeline_steps'] = ['ROI Extraction', 'SVG Conversion', 'Spatial Segmentation', 
                                  'Cell Annotation', 'Visualization']
        
        logging.info(f"Sample {sample_name} completed successfully")
        
    except Exception as e:
        error_msg = f"Pipeline failed for sample {sample_name}: {str(e)}"
        logging.error(error_msg)
        logging.error(traceback.format_exc())
        
        result['error_message'] = error_msg
        
        # Save failed configuration for debugging
        if hasattr(batch_config, 'failed_configs_dir'):
            failed_config_path = batch_config.failed_configs_dir / f"{sample_name}_failed_config.yaml"
            config.save(failed_config_path)
    
    finally:
        result['processing_time'] = time.time() - start_time
        result['end_time'] = datetime.now().isoformat()
    
    return result


def run_batch_processing(config_dir: str,
                        batch_config: BatchConfig,
                        validate_only: bool = False) -> Dict[str, Any]:
    """
    Run batch processing for multiple samples.
    
    Args:
        config_dir (str): Directory containing sample configuration files
        batch_config (BatchConfig): Batch processing configuration
        validate_only (bool): Only validate configurations without processing
        
    Returns:
        dict: Batch processing summary
    """
    setup_batch_logging(str(batch_config.base_output_dir))
    
    logging.info("=" * 80)
    logging.info("SPATIALCELL BATCH PROCESSING")
    logging.info("=" * 80)
    
    # Load sample configurations
    sample_configs = load_sample_configs(config_dir)
    
    logging.info(f"Loaded {len(sample_configs)} sample configurations")
    logging.info(f"Maximum parallel processes: {batch_config.max_parallel}")
    
    # Validation phase
    logging.info("\nValidating sample configurations...")
    valid_configs = []
    validation_errors = []
    
    for sample_name, config in sample_configs:
        try:
            validate_pipeline_inputs(config)
            valid_configs.append((sample_name, config))
            logging.info(f"✓ {sample_name}: Configuration valid")
        except Exception as e:
            error_msg = f"✗ {sample_name}: {str(e)}"
            validation_errors.append(error_msg)
            logging.error(error_msg)
    
    if validation_errors:
        logging.warning(f"Found {len(validation_errors)} configuration errors")
        if not batch_config.continue_on_error:
            raise ValueError("Validation failed. Set continue_on_error=True to proceed with valid configs.")
    
    if not valid_configs:
        raise ValueError("No valid configurations found")
    
    logging.info(f"Proceeding with {len(valid_configs)} valid configurations")
    
    if validate_only:
        return {
            'total_samples': len(sample_configs),
            'valid_samples': len(valid_configs),
            'validation_errors': validation_errors,
            'validated_only': True
        }
    
    # Processing phase
    batch_start_time = time.time()
    
    # Prepare arguments for parallel processing
    process_args = [(sample_name, config, batch_config) for sample_name, config in valid_configs]
    
    results = []
    
    if batch_config.max_parallel == 1:
        # Sequential processing
        logging.info("Running sequential processing...")
        for args in process_args:
            result = process_single_sample(args)
            results.append(result)
    else:
        # Parallel processing
        logging.info(f"Running parallel processing with {batch_config.max_parallel} processes...")
        
        with ProcessPoolExecutor(max_workers=batch_config.max_parallel) as executor:
            # Submit all jobs
            future_to_sample = {
                executor.submit(process_single_sample, args): args[0] 
                for args in process_args
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_sample):
                sample_name = future_to_sample[future]
                try:
                    result = future.result()
                    results.append(result)
                    
                    if result['success']:
                        logging.info(f"✓ Completed: {sample_name} ({result['processing_time']:.1f}s)")
                    else:
                        logging.error(f"✗ Failed: {sample_name}")
                        
                except Exception as e:
                    logging.error(f"✗ Exception in {sample_name}: {e}")
                    results.append({
                        'sample_name': sample_name,
                        'success': False,
                        'error_message': str(e),
                        'processing_time': 0
                    })
    
    # Generate batch summary
    total_time = time.time() - batch_start_time
    successful_samples = [r for r in results if r['success']]
    failed_samples = [r for r in results if not r['success']]
    
    batch_summary = {
        'batch_start_time': datetime.now().isoformat(),
        'total_processing_time': total_time,
        'total_samples': len(sample_configs),
        'valid_samples': len(valid_configs),
        'successful_samples': len(successful_samples),
        'failed_samples': len(failed_samples),
        'success_rate': len(successful_samples) / len(valid_configs) if valid_configs else 0,
        'validation_errors': validation_errors,
        'processing_results': results,
        'batch_config': {
            'max_parallel': batch_config.max_parallel,
            'continue_on_error': batch_config.continue_on_error,
            'base_output_dir': str(batch_config.base_output_dir)
        }
    }
    
    # Save batch summary
    if batch_config.generate_summary:
        summary_file = batch_config.batch_reports_dir / f"batch_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(summary_file, 'w') as f:
            json.dump(batch_summary, f, indent=2)
        
        logging.info(f"Batch summary saved to: {summary_file}")
    
    # Log final results
    logging.info("\n" + "=" * 80)
    logging.info("BATCH PROCESSING COMPLETED")
    logging.info("=" * 80)
    logging.info(f"Total samples: {len(sample_configs)}")
    logging.info(f"Valid configurations: {len(valid_configs)}")
    logging.info(f"Successful: {len(successful_samples)}")
    logging.info(f"Failed: {len(failed_samples)}")
    logging.info(f"Success rate: {batch_summary['success_rate']:.1%}")
    logging.info(f"Total time: {total_time:.1f} seconds")
    
    if failed_samples:
        logging.info("\nFailed samples:")
        for result in failed_samples:
            logging.info(f"  - {result['sample_name']}: {result.get('error_message', 'Unknown error')}")
    
    return batch_summary


def generate_batch_report(batch_summary: Dict[str, Any], output_path: str) -> None:
    """
    Generate detailed batch processing report.
    
    Args:
        batch_summary (dict): Batch processing summary
        output_path (str): Output path for report
    """
    with open(output_path, 'w') as f:
        f.write("# Spatialcell Batch Processing Report\n\n")
        
        f.write(f"**Processing Date**: {batch_summary['batch_start_time']}\n")
        f.write(f"**Total Processing Time**: {batch_summary['total_processing_time']:.1f} seconds\n")
        f.write(f"**Success Rate**: {batch_summary['success_rate']:.1%}\n\n")
        
        f.write("## Summary\n\n")
        f.write(f"- Total samples: {batch_summary['total_samples']}\n")
        f.write(f"- Valid configurations: {batch_summary['valid_samples']}\n")
        f.write(f"- Successful: {batch_summary['successful_samples']}\n")
        f.write(f"- Failed: {batch_summary['failed_samples']}\n\n")
        
        if batch_summary['validation_errors']:
            f.write("## Validation Errors\n\n")
            for error in batch_summary['validation_errors']:
                f.write(f"- {error}\n")
            f.write("\n")
        
        f.write("## Processing Results\n\n")
        f.write("| Sample | Status | Time (s) | Error |\n")
        f.write("|--------|--------|----------|-------|\n")
        
        for result in batch_summary['processing_results']:
            status = "✓ Success" if result['success'] else "✗ Failed"
            time_str = f"{result['processing_time']:.1f}"
            error = result.get('error_message', '')[:50] + '...' if result.get('error_message') else ''
            f.write(f"| {result['sample_name']} | {status} | {time_str} | {error} |\n")
        
        f.write("\n## Configuration\n\n")
        f.write("```yaml\n")
        f.write(yaml.dump(batch_summary['batch_config'], default_flow_style=False))
        f.write("```\n")


def main():
    """Main entry point for batch processing."""
    parser = argparse.ArgumentParser(
        description="Batch process multiple samples through Spatialcell pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all configs in directory
  python batch_processor.py --config_dir ./sample_configs/ --output_dir ./batch_output/

  # Parallel processing with 8 processes
  python batch_processor.py --config_dir ./configs/ --output_dir ./output/ --max_parallel 8

  # Validate configurations only
  python batch_processor.py --config_dir ./configs/ --validate_only

  # Continue processing even if some samples fail
  python batch_processor.py --config_dir ./configs/ --output_dir ./output/ --continue_on_error
        """
    )
    
    parser.add_argument("--config_dir", required=True,
                        help="Directory containing sample configuration files")
    parser.add_argument("--output_dir", required=True,
                        help="Base output directory for batch processing")
    parser.add_argument("--max_parallel", type=int, default=4,
                        help="Maximum number of parallel processes (default: 4)")
    parser.add_argument("--continue_on_error", action='store_true',
                        help="Continue processing other samples if one fails")
    parser.add_argument("--validate_only", action='store_true',
                        help="Only validate configurations without processing")
    parser.add_argument("--no_summary", action='store_true',
                        help="Skip generating batch summary report")
    
    args = parser.parse_args()
    
    try:
        # Create batch configuration
        batch_config = BatchConfig(
            base_output_dir=args.output_dir,
            max_parallel=args.max_parallel,
            continue_on_error=args.continue_on_error,
            generate_summary=not args.no_summary
        )
        
        # Save batch config
        batch_config.save_config(batch_config.batch_reports_dir / "batch_config.yaml")
        
        # Run batch processing
        summary = run_batch_processing(
            args.config_dir,
            batch_config,
            args.validate_only
        )
        
        # Generate detailed report
        if not args.no_summary and not args.validate_only:
            report_path = batch_config.batch_reports_dir / "detailed_report.md"
            generate_batch_report(summary, report_path)
            print(f"Detailed report saved to: {report_path}")
        
        # Print final summary
        if args.validate_only:
            print(f"\nValidation completed:")
            print(f"  Valid configurations: {summary['valid_samples']}/{summary['total_samples']}")
        else:
            print(f"\nBatch processing completed:")
            print(f"  Success rate: {summary['success_rate']:.1%}")
            print(f"  Results in: {args.output_dir}")
        
        return 0 if summary.get('success_rate', 0) > 0.5 else 1
        
    except Exception as e:
        print(f"Batch processing failed: {e}")
        return 1


if __name__ == "__main__":
    exit(main())# Version: 1.0.7
