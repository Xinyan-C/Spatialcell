#!/usr/bin/env python3
"""
Complete End-to-End Pipeline for Spatialcell Analysis

This script orchestrates the entire Spatialcell workflow from ROI extraction
through spatial segmentation to cell type annotation and visualization.

Built upon:
- TopAct framework: https://gitlab.com/kfbenjamin/topact.git
- bin2cell package: https://github.com/Teichlab/bin2cell.git

Workflow:
1. ROI coordinate extraction (optional)
2. SVG to NPZ conversion  
3. Spatial segmentation with bin2cell
4. Cell type classification with TopAct
5. Comprehensive visualization

Author: Xinyan
License: MIT
"""

import os
import sys
import argparse
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
import yaml
from datetime import datetime

# Add module paths
sys.path.append(str(Path(__file__).parent.parent))

from ..utils.roi_extractor import extract_roi_coordinates
from ..preprocessing.svg_to_npz import convert_svg_to_npz
from ..spatial_segmentation.spatial_processor import process_spatial_data
from ..cell_annotation.annotation_processor import process_sample_annotation
from .pipeline_config import PipelineConfig, validate_pipeline_inputs


def setup_pipeline_logging(output_dir: str, sample_name: str) -> None:
    """
    Setup comprehensive logging for the complete pipeline.
    
    Args:
        output_dir (str): Output directory for logs
        sample_name (str): Sample identifier for log naming
    """
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    
    log_file = log_dir / f"pipeline_{sample_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    # Configure root logger
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ],
        force=True
    )
    
    logging.info(f"Pipeline logging initialized for sample: {sample_name}")
    logging.info(f"Log file: {log_file}")


def run_roi_extraction(config: PipelineConfig) -> str:
    """
    Extract ROI coordinates from Loupe Browser exports.
    
    Args:
        config (PipelineConfig): Pipeline configuration
        
    Returns:
        str: Path to generated ROI coordinate file
        
    Raises:
        RuntimeError: If ROI extraction fails
    """
    if not config.enable_roi_extraction:
        logging.info("ROI extraction disabled, using existing coordinate file")
        return config.roi_coordinate_file
    
    logging.info("Starting ROI coordinate extraction...")
    
    try:
        roi_ranges = extract_roi_coordinates(
            sample_name=config.sample_name,
            sample_dir=config.loupe_export_dir,
            output_path=config.roi_coordinate_file
        )
        
        logging.info(f"ROI extraction completed: {len(roi_ranges)} regions extracted")
        return config.roi_coordinate_file
        
    except Exception as e:
        raise RuntimeError(f"ROI extraction failed: {e}")


def run_svg_conversion(config: PipelineConfig) -> str:
    """
    Convert QuPath SVG exports to NPZ format.
    
    Args:
        config (PipelineConfig): Pipeline configuration
        
    Returns:
        str: Path to generated NPZ file
        
    Raises:
        RuntimeError: If SVG conversion fails
    """
    if not config.enable_svg_conversion:
        logging.info("SVG conversion disabled, using existing NPZ file")
        return config.npz_labels_file
    
    logging.info("Starting SVG to NPZ conversion...")
    
    try:
        num_objects = convert_svg_to_npz(
            svg_path=config.svg_file,
            height=config.image_height,
            width=config.image_width,
            output_path=config.npz_labels_file,
            verbose=config.verbose
        )
        
        logging.info(f"SVG conversion completed: {num_objects} objects converted")
        return config.npz_labels_file
        
    except Exception as e:
        raise RuntimeError(f"SVG conversion failed: {e}")


def run_spatial_segmentation(config: PipelineConfig, roi_file: str, npz_file: str) -> Tuple[str, str]:
    """
    Run spatial segmentation using bin2cell integration.
    
    Args:
        config (PipelineConfig): Pipeline configuration
        roi_file (str): Path to ROI coordinate file
        npz_file (str): Path to NPZ labels file
        
    Returns:
        tuple: (adata_file_path, cdata_file_path)
        
    Raises:
        RuntimeError: If spatial segmentation fails
    """
    logging.info("Starting spatial segmentation with bin2cell...")
    
    # Create segmentation output directory
    seg_output_dir = Path(config.output_dir) / "spatial_segmentation"
    seg_output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Create arguments object for spatial segmentation
        class SegmentationArgs:
            def __init__(self):
                self.path = config.visium_data_path
                self.source_image_path = config.source_image_path
                self.region_file = roi_file
                self.npz_path = npz_file
                self.output_dir = str(seg_output_dir)
                self.sample = config.sample_name
                self.prob_thresh = config.stardist_prob_thresh
                self.nms_thresh = config.stardist_nms_thresh
                self.labels_key = config.labels_key
                self.algorithm = config.expansion_algorithm
                self.max_bin_distance = config.max_bin_distance
                self.volume_ratio = config.volume_ratio
                self.k = config.nearest_neighbors
                self.subset_pca = config.subset_pca
        
        seg_args = SegmentationArgs()
        
        # Run spatial segmentation
        adata, cdata = process_spatial_data(seg_args)
        
        # Define output file paths
        data_dir = seg_output_dir / "Data"
        adata_file = data_dir / f"{config.sample_name}_2um.h5ad"
        cdata_file = data_dir / f"{config.sample_name}_b2c.h5ad"
        
        logging.info(f"Spatial segmentation completed")
        logging.info(f"Spot-level data: {adata_file}")
        logging.info(f"Cell-level data: {cdata_file}")
        
        return str(adata_file), str(cdata_file)
        
    except Exception as e:
        raise RuntimeError(f"Spatial segmentation failed: {e}")


def run_cell_annotation(config: PipelineConfig, roi_file: str, cdata_file: str) -> str:
    """
    Run cell type annotation using TopAct framework.
    
    Args:
        config (PipelineConfig): Pipeline configuration
        roi_file (str): Path to ROI coordinate file
        cdata_file (str): Path to cell-level data from segmentation
        
    Returns:
        str: Path to annotation results directory
        
    Raises:
        RuntimeError: If cell annotation fails
    """
    if not config.enable_cell_annotation:
        logging.info("Cell annotation disabled, skipping...")
        return ""
    
    logging.info("Starting cell type annotation with TopAct...")
    
    # Create annotation output directory
    annotation_output_dir = Path(config.output_dir) / "cell_annotation"
    annotation_output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Create arguments object for annotation
        class AnnotationArgs:
            def __init__(self):
                self.sample = config.sample_name
                self.out_dir = str(annotation_output_dir)
                self.expr_path = config.hd_expression_path
                self.pos_path = config.hd_positions_path
                self.source_image_path = config.source_image_path
                self.roi_file = roi_file
                self.clf_path = config.classifier_path
                self.bin2cell_dir = cdata_file
                self.min_scale = config.min_scale
                self.max_scale = config.max_scale
                self.num_proc = config.num_processes
                self.threads = str(config.num_threads)
                self.mem_gb = config.memory_limit_gb
                self.labels = config.cell_labels_column
        
        annotation_args = AnnotationArgs()
        
        # Run cell annotation
        process_sample_annotation(annotation_args)
        
        logging.info(f"Cell annotation completed")
        logging.info(f"Results saved to: {annotation_output_dir}")
        
        return str(annotation_output_dir)
        
    except Exception as e:
        raise RuntimeError(f"Cell annotation failed: {e}")


def run_visualization(config: PipelineConfig, annotation_dir: str, roi_file: str) -> str:
    """
    Generate comprehensive visualizations of results.
    
    Args:
        config (PipelineConfig): Pipeline configuration
        annotation_dir (str): Path to annotation results
        roi_file (str): Path to ROI coordinate file
        
    Returns:
        str: Path to visualization output directory
        
    Raises:
        RuntimeError: If visualization fails
    """
    if not config.enable_visualization:
        logging.info("Visualization disabled, skipping...")
        return ""
    
    logging.info("Starting comprehensive visualization...")
    
    # Create visualization output directory
    viz_output_dir = Path(config.output_dir) / "visualizations"
    viz_output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Run visualization script
        cmd = [
            sys.executable, 
            str(Path(__file__).parent.parent / "cell_annotation" / "annotation_visualizer.py"),
            "--sample", config.sample_name,
            "--sd_dir", annotation_dir,
            "--outfile_dir", annotation_dir,
            "--clf_dir", str(Path(config.classifier_path).parent),
            "--roi_file", roi_file,
            "--output_dir", str(viz_output_dir),
            "--color_scheme", config.color_scheme,
            "--point_size", str(config.point_size),
            "--rename_cell_types", str(config.rename_cell_types)
        ]
        
        # Add background image if provided
        if config.background_image_path:
            cmd.extend(["--background_image", config.background_image_path])
        
        # Add verbose flag if enabled
        if config.verbose:
            cmd.append("--verbose")
        
        # Run visualization
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        logging.info("Visualization completed successfully")
        logging.info(f"Results saved to: {viz_output_dir}")
        
        return str(viz_output_dir)
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Visualization command failed: {e.stderr}")
        raise RuntimeError(f"Visualization failed: {e}")
    except Exception as e:
        raise RuntimeError(f"Visualization failed: {e}")


def generate_pipeline_report(config: PipelineConfig, results: Dict[str, str]) -> str:
    """
    Generate comprehensive pipeline execution report.
    
    Args:
        config (PipelineConfig): Pipeline configuration
        results (dict): Dictionary of step results and output paths
        
    Returns:
        str: Path to generated report
    """
    report_dir = Path(config.output_dir) / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    
    report_file = report_dir / f"pipeline_report_{config.sample_name}.md"
    
    with open(report_file, 'w') as f:
        f.write(f"# Spatialcell Pipeline Report\n\n")
        f.write(f"**Sample**: {config.sample_name}\n")
        f.write(f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"**Output Directory**: {config.output_dir}\n\n")
        
        f.write("## Configuration\n\n")
        f.write("```yaml\n")
        f.write(config.to_yaml())
        f.write("```\n\n")
        
        f.write("## Execution Results\n\n")
        for step, result in results.items():
            f.write(f"- **{step}**: {result}\n")
        
        f.write("\n## Output Structure\n\n")
        f.write("```\n")
        f.write(f"{config.output_dir}/\n")
        f.write("├── spatial_segmentation/\n")
        f.write("├── cell_annotation/\n")
        f.write("├── visualizations/\n")
        f.write("├── reports/\n")
        f.write("└── logs/\n")
        f.write("```\n")
        
        f.write("\n## Next Steps\n\n")
        f.write("1. Review visualization outputs in `visualizations/`\n")
        f.write("2. Examine detailed logs in `logs/`\n")
        f.write("3. Analyze cell-level data in `spatial_segmentation/Data/`\n")
        f.write("4. Explore classification results in `cell_annotation/`\n")
    
    logging.info(f"Pipeline report generated: {report_file}")
    return str(report_file)


def run_complete_pipeline(config_file: str) -> None:
    """
    Execute the complete Spatialcell pipeline.
    
    Args:
        config_file (str): Path to pipeline configuration file
        
    Raises:
        RuntimeError: If any pipeline step fails
    """
    # Load configuration
    config = PipelineConfig.from_file(config_file)
    
    # Setup logging
    setup_pipeline_logging(config.output_dir, config.sample_name)
    
    # Validate inputs
    validate_pipeline_inputs(config)
    
    logging.info("=" * 80)
    logging.info("SPATIALCELL COMPLETE PIPELINE EXECUTION")
    logging.info("=" * 80)
    logging.info(f"Sample: {config.sample_name}")
    logging.info(f"Output: {config.output_dir}")
    
    results = {}
    
    try:
        # Step 1: ROI Extraction (if enabled)
        logging.info("\n" + "=" * 40)
        logging.info("STEP 1: ROI COORDINATE EXTRACTION")
        logging.info("=" * 40)
        roi_file = run_roi_extraction(config)
        results["ROI Extraction"] = roi_file
        
        # Step 2: SVG Conversion (if enabled)
        logging.info("\n" + "=" * 40)
        logging.info("STEP 2: SVG TO NPZ CONVERSION")
        logging.info("=" * 40)
        npz_file = run_svg_conversion(config)
        results["SVG Conversion"] = npz_file
        
        # Step 3: Spatial Segmentation
        logging.info("\n" + "=" * 40)
        logging.info("STEP 3: SPATIAL SEGMENTATION")
        logging.info("=" * 40)
        adata_file, cdata_file = run_spatial_segmentation(config, roi_file, npz_file)
        results["Spatial Segmentation"] = f"adata: {adata_file}, cdata: {cdata_file}"
        
        # Step 4: Cell Annotation (if enabled)
        logging.info("\n" + "=" * 40)
        logging.info("STEP 4: CELL TYPE ANNOTATION")
        logging.info("=" * 40)
        annotation_dir = run_cell_annotation(config, roi_file, cdata_file)
        if annotation_dir:
            results["Cell Annotation"] = annotation_dir
        
        # Step 5: Visualization (if enabled)
        logging.info("\n" + "=" * 40)
        logging.info("STEP 5: VISUALIZATION")
        logging.info("=" * 40)
        viz_dir = run_visualization(config, annotation_dir, roi_file)
        if viz_dir:
            results["Visualization"] = viz_dir
        
        # Generate final report
        logging.info("\n" + "=" * 40)
        logging.info("GENERATING PIPELINE REPORT")
        logging.info("=" * 40)
        report_file = generate_pipeline_report(config, results)
        results["Report"] = report_file
        
        logging.info("\n" + "=" * 80)
        logging.info("PIPELINE EXECUTION COMPLETED SUCCESSFULLY")
        logging.info("=" * 80)
        
        for step, result in results.items():
            logging.info(f"{step}: {result}")
            
    except Exception as e:
        logging.error(f"\nPIPELINE EXECUTION FAILED: {e}")
        logging.error("Check logs for detailed error information")
        raise


def main():
    """Main entry point for complete pipeline execution."""
    parser = argparse.ArgumentParser(
        description="Execute complete Spatialcell pipeline from configuration file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline
  python complete_pipeline.py --config sample_config.yaml

  # Run with custom output directory
  python complete_pipeline.py --config sample_config.yaml --output_dir /custom/output/

Configuration file should be in YAML format with all required parameters.
See examples/ directory for sample configuration files.
        """
    )
    
    parser.add_argument("--config", required=True,
                        help="Path to pipeline configuration file (YAML or JSON)")
    parser.add_argument("--output_dir", default=None,
                        help="Override output directory from config")
    parser.add_argument("--validate_only", action='store_true',
                        help="Only validate configuration without running pipeline")
    
    args = parser.parse_args()
    
    try:
        # Load and validate configuration
        config = PipelineConfig.from_file(args.config)
        
        # Override output directory if provided
        if args.output_dir:
            config.output_dir = args.output_dir
        
        if args.validate_only:
            print("Configuration validation...")
            validate_pipeline_inputs(config)
            print("Configuration is valid!")
            return
        
        # Run complete pipeline
        run_complete_pipeline(args.config)
        
        print("\nPipeline completed successfully!")
        print(f"Results available in: {config.output_dir}")
        
    except Exception as e:
        print(f"Pipeline failed: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())# Version: 1.0.7
