#!/usr/bin/env python3
"""
RDS Data Processor for TopAct Classifier Training

This module provides Python interface for extracting training data from Seurat RDS files
to prepare data for TopAct classifier training.

Author: Xinyan
License: Apache 2.0
"""

import os
import subprocess
import argparse
import logging
from pathlib import Path
from typing import List, Optional, Dict


def setup_logging(verbose: bool = False) -> None:
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(levelname)s - %(message)s'
    )


def check_r_dependencies() -> bool:
    """Check if required R packages are installed."""
    try:
        # Check R installation
        result = subprocess.run(['R', '--version'], 
                              capture_output=True, text=True, check=True)
        logging.info("R installation found")
        
        # Check required packages
        r_check_script = """
        packages <- c('Seurat', 'Matrix', 'argparse')
        missing <- packages[!packages %in% installed.packages()[,'Package']]
        if(length(missing) > 0) {
            cat('Missing packages:', paste(missing, collapse=', '), '\\n')
            quit(status=1)
        } else {
            cat('All required packages available\\n')
        }
        """
        
        result = subprocess.run(['R', '--slave', '-e', r_check_script],
                              capture_output=True, text=True)
        
        if result.returncode == 0:
            logging.info("All R dependencies available")
            return True
        else:
            logging.error(f"Missing R packages: {result.stdout}")
            return False
            
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logging.error(f"R dependency check failed: {e}")
        return False


def extract_training_data_from_rds(rds_file: str,
                                  output_dir: str,
                                  samples: Optional[List[str]] = None,
                                  celltype_col: str = "celltype_CS",
                                  verbose: bool = False) -> Dict[str, str]:
    """
    Extract training data from Seurat RDS file using R script.
    
    Args:
        rds_file (str): Path to Seurat RDS file
        output_dir (str): Output directory for training data
        samples (List[str], optional): Sample names to extract
        celltype_col (str): Cell type annotation column name
        verbose (bool): Enable verbose logging
        
    Returns:
        Dict[str, str]: Paths to generated files
        
    Raises:
        FileNotFoundError: If RDS file doesn't exist
        RuntimeError: If R script execution fails
    """
    setup_logging(verbose)
    
    # Validate inputs
    rds_path = Path(rds_file)
    if not rds_path.exists():
        raise FileNotFoundError(f"RDS file not found: {rds_file}")
    
    if not check_r_dependencies():
        raise RuntimeError("Required R dependencies not available. Please install R and required packages.")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Get R script path
    script_dir = Path(__file__).parent
    r_script_path = script_dir / "rds_extractor.R"
    
    if not r_script_path.exists():
        raise FileNotFoundError(f"R script not found: {r_script_path}")
    
    # Prepare R command
    cmd = [
        'Rscript', str(r_script_path),
        '--rds_file', str(rds_path.absolute()),
        '--output_dir', str(output_path.absolute()),
        '--celltype_col', celltype_col
    ]
    
    if samples:
        cmd.extend(['--sample'] + samples)
    
    logging.info(f"Executing R script: {' '.join(cmd)}")
    
    # Execute R script
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=str(script_dir))
        logging.info("R script execution completed successfully")
        if verbose:
            logging.debug(f"R script output: {result.stdout}")
        
    except subprocess.CalledProcessError as e:
        logging.error(f"R script execution failed: {e.stderr}")
        raise RuntimeError(f"RDS extraction failed: {e.stderr}")
    
    # Check generated files
    generated_files = {}
    
    if samples:
        for sample in samples:
            prefix = f"{sample}_all"
            files_dict = {
                'counts': str(output_path / f"{prefix}_counts.mtx"),
                'genes': str(output_path / f"{prefix}_genes.csv"),
                'barcodes': str(output_path / f"{prefix}_barcodes.csv"),
                'metadata': str(output_path / f"{prefix}_meta_data.csv"),
                'celltypes': str(output_path / f"{prefix}_celltypes.csv")
            }
            # Verify files exist
            existing_files = {k: v for k, v in files_dict.items() if Path(v).exists()}
            generated_files[sample] = existing_files
    else:
        files_dict = {
            'counts': str(output_path / "all_samples_counts.mtx"),
            'genes': str(output_path / "all_samples_genes.csv"),
            'barcodes': str(output_path / "all_samples_barcodes.csv"),
            'metadata': str(output_path / "all_samples_meta_data.csv"),
            'celltypes': str(output_path / "all_samples_celltypes.csv")
        }
        # Verify files exist
        existing_files = {k: v for k, v in files_dict.items() if Path(v).exists()}
        generated_files['all_samples'] = existing_files
    
    logging.info(f"Training data extracted to: {output_dir}")
    return generated_files


def main():
    """Main entry point for RDS data extraction."""
    parser = argparse.ArgumentParser(
        description="Extract training data from Seurat RDS files for TopAct classifier training",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract all samples
  python rds_processor.py --rds_file data.rds --output_dir ./training_data/

  # Extract specific time points
  python rds_processor.py --rds_file data.rds --output_dir ./training_data/ \\
    --samples E14.5 E18.5 P3

  # Custom cell type column
  python rds_processor.py --rds_file data.rds --output_dir ./training_data/ \\
    --celltype_col custom_celltype --verbose
        """
    )
    
    parser.add_argument('--rds_file', required=True,
                       help='Path to Seurat RDS file')
    parser.add_argument('--output_dir', required=True,
                       help='Output directory for training data')
    parser.add_argument('--samples', nargs='+', default=None,
                       help='Sample names to extract (e.g., E14.5 E18.5 P3)')
    parser.add_argument('--celltype_col', default='celltype_CS',
                       help='Cell type annotation column name (default: celltype_CS)')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    try:
        generated_files = extract_training_data_from_rds(
            args.rds_file,
            args.output_dir,
            args.samples,
            args.celltype_col,
            args.verbose
        )
        
        print("✓ Training data extraction completed successfully!")
        print(f"✓ Output directory: {args.output_dir}")
        
        if args.samples:
            for sample in args.samples:
                if sample in generated_files:
                    print(f"✓ Sample {sample}: {len(generated_files[sample])} files generated")
        else:
            if 'all_samples' in generated_files:
                print(f"✓ All samples: {len(generated_files['all_samples'])} files generated")
            
    except Exception as e:
        logging.error(f"Extraction failed: {e}")
        print(f"✗ Training data extraction failed: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())