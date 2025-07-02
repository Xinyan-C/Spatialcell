#!/usr/bin/env python3
"""
ROI Coordinate Extraction Tool for Spatialcell Pipeline

This tool extracts Region of Interest (ROI) coordinates from Loupe Browser 
exported CSV files and generates coordinate range files for downstream analysis.

Author: Xinyan
License: MIT
"""

import pandas as pd
import argparse
import os
import logging
from pathlib import Path
from typing import Tuple, Optional, Dict, List


def setup_logging(verbose=False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def calculate_roi_range(roi_csv_path: str, total_coords: pd.DataFrame, 
                       sample_name: str, roi_label: str) -> Tuple[Optional[float], ...]:
    """
    Calculate coordinate range for a specified ROI rectangle.
    
    Args:
        roi_csv_path (str): Path to ROI CSV file from Loupe Browser
        total_coords (pd.DataFrame): DataFrame containing all coordinates
        sample_name (str): Sample name
        roi_label (str): ROI label (e.g., CS1, CS2, WT1)
        
    Returns:
        tuple: (x_min, x_max, y_min, y_max) coordinate ranges
    """
    try:
        if not Path(roi_csv_path).exists():
            logging.warning(f"ROI file not found: {roi_csv_path}")
            return None, None, None, None
            
        # Load ROI CSV file
        roi_data = pd.read_csv(roi_csv_path)
        
        # Extract barcodes with values in the second column (index 1), indicating selection
        selected_barcodes = roi_data[roi_data.iloc[:, 1].notna()]['Barcode'].tolist()
        
        if not selected_barcodes:
            logging.warning(f"No selected barcodes found in {roi_label}")
            return None, None, None, None
        
        # Filter total_coords for selected barcodes
        selected_coords = total_coords[total_coords['Barcode'].isin(selected_barcodes)]
        
        if selected_coords.empty:
            logging.warning(f"No matching barcodes found for {roi_label}")
            return None, None, None, None
        
        # Calculate min and max coordinates
        x_min = selected_coords['X Coordinate'].min()
        x_max = selected_coords['X Coordinate'].max()
        y_min = selected_coords['Y Coordinate'].min()
        y_max = selected_coords['Y Coordinate'].max()
        
        # Log results
        logging.info(f"{sample_name} - {roi_label} Rectangle coordinates:")
        logging.info(f"  X: {x_min} - {x_max}")
        logging.info(f"  Y: {y_min} - {y_max}")
        
        return x_min, x_max, y_min, y_max
        
    except Exception as e:
        logging.error(f"Error processing {sample_name} - {roi_label}: {e}")
        return None, None, None, None


def detect_roi_files(sample_dir: str, sample_name: str) -> List[Tuple[str, str]]:
    """
    Automatically detect ROI CSV files in the sample directory.
    
    Args:
        sample_dir (str): Sample directory path
        sample_name (str): Sample name
        
    Returns:
        list: List of (roi_file_path, roi_label) tuples
    """
    roi_files = []
    sample_dir_path = Path(sample_dir)
    
    # Common ROI patterns
    roi_patterns = [
        ('CS1', f"{sample_name}-CS1.csv"),
        ('CS2', f"{sample_name}-CS2.csv"), 
        ('CS3', f"{sample_name}-CS3.csv"),
        ('WT1', f"{sample_name}-WT1.csv"),
        ('WT2', f"{sample_name}-WT2.csv"),
        ('WT3', f"{sample_name}-WT3.csv"),
    ]
    
    for roi_label, filename in roi_patterns:
        file_path = sample_dir_path / filename
        if file_path.exists():
            roi_files.append((str(file_path), roi_label))
            logging.info(f"Found ROI file: {filename}")
    
    # Also look for other CSV files that might be ROIs
    for csv_file in sample_dir_path.glob(f"{sample_name}-*.csv"):
        if csv_file.name not in [pattern[1] for pattern in roi_patterns]:
            # Extract ROI label from filename
            roi_label = csv_file.stem.replace(f"{sample_name}-", "")
            roi_files.append((str(csv_file), roi_label))
            logging.info(f"Found additional ROI file: {csv_file.name}")
    
    return roi_files


def extract_roi_coordinates(sample_name: str, sample_dir: str, 
                          output_path: Optional[str] = None,
                          roi_files: Optional[List[Tuple[str, str]]] = None) -> Dict:
    """
    Extract ROI coordinates from Loupe Browser exported CSV files.
    
    Args:
        sample_name (str): Sample name
        sample_dir (str): Directory containing sample files
        output_path (str, optional): Output file path for coordinate ranges
        roi_files (list, optional): List of (roi_file_path, roi_label) tuples
        
    Returns:
        dict: Dictionary containing ROI coordinate ranges
    """
    sample_dir_path = Path(sample_dir)
    
    # Check if sample directory exists
    if not sample_dir_path.exists():
        raise FileNotFoundError(f"Sample directory not found: {sample_dir}")
    
    # Load total coordinates file
    total_coords_path = sample_dir_path / f"{sample_name}-all.csv"
    if not total_coords_path.exists():
        raise FileNotFoundError(f"Total coordinates file not found: {total_coords_path}")
    
    try:
        total_coords = pd.read_csv(total_coords_path, 
                                 usecols=['Barcode', 'X Coordinate', 'Y Coordinate'])
        logging.info(f"Loaded {len(total_coords)} total coordinates")
    except Exception as e:
        raise ValueError(f"Failed to load {total_coords_path}: {e}")
    
    # Auto-detect ROI files if not provided
    if roi_files is None:
        roi_files = detect_roi_files(sample_dir, sample_name)
    
    if not roi_files:
        logging.warning("No ROI files found")
        return {}
    
    # Extract coordinates for each ROI
    roi_ranges = {}
    for roi_file_path, roi_label in roi_files:
        x_min, x_max, y_min, y_max = calculate_roi_range(
            roi_file_path, total_coords, sample_name, roi_label
        )
        
        if all(coord is not None for coord in [x_min, x_max, y_min, y_max]):
            roi_ranges[f"{sample_name} - {roi_label}"] = {
                'x_min': x_min, 'x_max': x_max,
                'y_min': y_min, 'y_max': y_max
            }
    
    # Calculate overall image range
    overall_x_min = total_coords['X Coordinate'].min()
    overall_x_max = total_coords['X Coordinate'].max()
    overall_y_min = total_coords['Y Coordinate'].min()
    overall_y_max = total_coords['Y Coordinate'].max()
    
    roi_ranges[f"{sample_name} whole image"] = {
        'x_min': overall_x_min, 'x_max': overall_x_max,
        'y_min': overall_y_min, 'y_max': overall_y_max
    }
    
    logging.info(f"{sample_name} whole image coordinates:")
    logging.info(f"  X: {overall_x_min} - {overall_x_max}")
    logging.info(f"  Y: {overall_y_min} - {overall_y_max}")
    
    # Save results to file
    if output_path is None:
        output_path = sample_dir_path / f"{sample_name}_ranges.txt"
    
    save_roi_ranges(roi_ranges, output_path)
    
    return roi_ranges


def save_roi_ranges(roi_ranges: Dict, output_path: str):
    """
    Save ROI coordinate ranges to a text file.
    
    Args:
        roi_ranges (dict): Dictionary containing ROI coordinate ranges
        output_path (str): Output file path
    """
    with open(output_path, 'w', encoding='utf-8') as f:
        for roi_name, coords in roi_ranges.items():
            f.write(f"{roi_name} Rectangle coordinates:\n")
            f.write(f"  X: {coords['x_min']} - {coords['x_max']}\n")
            f.write(f"  Y: {coords['y_min']} - {coords['y_max']}\n\n")
    
    logging.info(f"Results saved to {output_path}")


def main():
    """Main entry point for ROI coordinate extraction."""
    parser = argparse.ArgumentParser(
        description="Extract ROI coordinates from Loupe Browser exported CSV files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python roi_extractor.py --sample E14.5
  python roi_extractor.py --sample P3 --base_dir /data/coordinates/
  python roi_extractor.py --sample custom_sample --base_dir /data/ --verbose

Expected file structure:
  base_dir/
  └── sample_name/
      ├── sample_name-all.csv      (total coordinates)
      ├── sample_name-CS1.csv      (ROI 1)
      ├── sample_name-CS2.csv      (ROI 2)
      └── sample_name-WT1.csv      (ROI 3)
        """
    )
    
    parser.add_argument('--sample', required=True,
                        help="Sample name (e.g., E14.5, E18.5, P3)")
    parser.add_argument('--base_dir', default='.',
                        help="Base directory containing sample folders (default: current directory)")
    parser.add_argument('--output', default=None,
                        help="Output file path (default: sample_dir/sample_ranges.txt)")
    parser.add_argument('--verbose', action='store_true',
                        help="Enable verbose logging")
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    try:
        # Construct sample directory path
        sample_dir = os.path.join(args.base_dir, args.sample)
        
        # Extract ROI coordinates
        roi_ranges = extract_roi_coordinates(
            args.sample, sample_dir, args.output
        )
        
        if roi_ranges:
            print(f"✓ Successfully extracted coordinates for {len(roi_ranges)} ROIs")
            print(f"✓ Results saved to coordinate range file")
        else:
            print("⚠ No ROI coordinates extracted")
            
    except Exception as e:
        logging.error(f"Extraction failed: {e}")
        print(f"✗ Failed to extract ROI coordinates: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())# Version: 1.0.7
