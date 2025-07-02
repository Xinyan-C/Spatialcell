#!/usr/bin/env python3
"""
ROI Coordinate Validation Tool for Spatialcell Pipeline

This tool validates ROI coordinates by visualizing them overlaid on source images,
helping users verify that ROI extraction was performed correctly.

Author: Xinyan
License: MIT
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import argparse
import os
import logging
from pathlib import Path
from typing import Dict, Optional, List
import warnings

# Handle different image reading libraries
try:
    from tifffile import imread
    TIFF_AVAILABLE = True
except ImportError:
    try:
        from PIL import Image
        TIFF_AVAILABLE = False
    except ImportError:
        raise ImportError("Please install either tifffile or Pillow for image reading")


def setup_logging(verbose=False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def read_image(image_path: str) -> np.ndarray:
    """
    Read image from file using available libraries.
    
    Args:
        image_path (str): Path to image file
        
    Returns:
        np.ndarray: Image array
    """
    if not Path(image_path).exists():
        raise FileNotFoundError(f"Image file not found: {image_path}")
    
    try:
        if TIFF_AVAILABLE:
            return imread(image_path)
        else:
            # Use PIL as fallback
            with Image.open(image_path) as img:
                return np.array(img)
    except Exception as e:
        raise ValueError(f"Failed to read image {image_path}: {e}")


def parse_coordinate_file(coords_file: str) -> Dict:
    """
    Parse coordinate file and extract ROI information.
    
    Args:
        coords_file (str): Path to coordinate range file
        
    Returns:
        dict: Dictionary containing ROI coordinates
    """
    if not Path(coords_file).exists():
        raise FileNotFoundError(f"Coordinate file not found: {coords_file}")
    
    rois = {}
    
    try:
        with open(coords_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith('#'):  # Skip empty lines and comments
                i += 1
                continue
            
            # Identify ROI name (contains coordinate indicators)
            coord_indicators = [
                'Rectangle coordinates:',
                '矩形框的坐标范围:',
                'coordinates:',
                'whole image'
            ]
            
            if any(indicator in line for indicator in coord_indicators):
                # Extract ROI name
                roi_name = line
                for indicator in coord_indicators:
                    if indicator in line:
                        roi_name = line.split(indicator)[0].strip()
                        break
                
                roi_coords = {}
                
                # Parse X and Y coordinates from following lines
                i += 1
                while i < len(lines) and lines[i].strip():  # Until empty line or EOF
                    next_line = lines[i].strip()
                    
                    if 'X:' in next_line:
                        x_part = next_line.split('X:')[1].strip()
                        x_min, x_max = map(float, x_part.split(' - '))
                        roi_coords['x_min'] = x_min
                        roi_coords['x_max'] = x_max
                        
                    elif 'Y:' in next_line:
                        y_part = next_line.split('Y:')[1].strip()
                        y_min, y_max = map(float, y_part.split(' - '))
                        roi_coords['y_min'] = y_min
                        roi_coords['y_max'] = y_max
                        
                    i += 1
                
                # Only add if we have all coordinates
                if len(roi_coords) == 4:
                    rois[roi_name] = roi_coords
                    logging.debug(f"Parsed ROI: {roi_name}")
                
            else:
                i += 1
                
    except Exception as e:
        raise ValueError(f"Error parsing coordinate file: {e}")
    
    logging.info(f"Successfully parsed {len(rois)} ROIs from coordinate file")
    return rois


def validate_roi_coordinates(rois: Dict, image_shape: tuple) -> List[str]:
    """
    Validate ROI coordinates against image dimensions.
    
    Args:
        rois (dict): Dictionary containing ROI coordinates
        image_shape (tuple): Image shape (height, width)
        
    Returns:
        list: List of validation warnings
    """
    warnings_list = []
    img_height, img_width = image_shape[:2]
    
    for roi_name, coords in rois.items():
        # Check if coordinates are within image bounds
        if coords['x_min'] < 0 or coords['x_max'] >= img_width:
            warnings_list.append(f"{roi_name}: X coordinates outside image bounds")
            
        if coords['y_min'] < 0 or coords['y_max'] >= img_height:
            warnings_list.append(f"{roi_name}: Y coordinates outside image bounds")
            
        # Check if min < max
        if coords['x_min'] >= coords['x_max']:
            warnings_list.append(f"{roi_name}: Invalid X range (min >= max)")
            
        if coords['y_min'] >= coords['y_max']:
            warnings_list.append(f"{roi_name}: Invalid Y range (min >= max)")
            
        # Check for reasonable ROI size
        roi_width = coords['x_max'] - coords['x_min']
        roi_height = coords['y_max'] - coords['y_min']
        
        if roi_width < 10 or roi_height < 10:
            warnings_list.append(f"{roi_name}: ROI too small ({roi_width}x{roi_height})")
    
    return warnings_list


def plot_image_with_rois(image: np.ndarray, rois: Dict, 
                        output_path: Optional[str] = None,
                        figsize: tuple = (12, 10),
                        dpi: int = 300) -> None:
    """
    Plot image with ROI overlays.
    
    Args:
        image (np.ndarray): Image array
        rois (dict): Dictionary containing ROI coordinates
        output_path (str, optional): Output file path
        figsize (tuple): Figure size
        dpi (int): DPI for saved figure
    """
    # Filter out whole image ROI for cleaner visualization
    plot_rois = {name: coords for name, coords in rois.items() 
                if 'whole image' not in name.lower()}
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Handle different image formats
    if len(image.shape) == 3:
        ax.imshow(image)
    else:
        ax.imshow(image, cmap='gray')
    
    # Define colors for different ROIs
    colors = ['red', 'green', 'blue', 'yellow', 'purple', 'orange', 'cyan', 'magenta']
    
    for i, (roi_name, coords) in enumerate(plot_rois.items()):
        x_min = coords['x_min']
        x_max = coords['x_max'] 
        y_min = coords['y_min']
        y_max = coords['y_max']
        
        # Create rectangle patch
        rect = Rectangle(
            (x_min, y_min), 
            x_max - x_min, 
            y_max - y_min,
            linewidth=2, 
            edgecolor=colors[i % len(colors)], 
            facecolor='none', 
            label=roi_name,
            alpha=0.8
        )
        ax.add_patch(rect)
        
        # Add text label
        ax.text(x_min + 10, y_min + 20, roi_name, 
               color=colors[i % len(colors)], fontsize=10, 
               bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.7))
    
    # Customize plot
    ax.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
    ax.set_title("Image with ROI Validation", fontsize=14, fontweight='bold')
    ax.axis('on')
    ax.grid(True, alpha=0.3)
    
    # Add image info
    img_info = f"Image size: {image.shape[1]}×{image.shape[0]} pixels"
    ax.text(0.02, 0.02, img_info, transform=ax.transAxes, 
           bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        logging.info(f"Validation plot saved to: {output_path}")
    else:
        plt.show()
    
    plt.close()


def generate_roi_report(rois: Dict, image_shape: tuple, 
                       sample_name: str, output_dir: str) -> None:
    """
    Generate a detailed ROI validation report.
    
    Args:
        rois (dict): Dictionary containing ROI coordinates
        image_shape (tuple): Image shape
        sample_name (str): Sample name
        output_dir (str): Output directory
    """
    report_path = Path(output_dir) / f"{sample_name}_roi_validation_report.txt"
    
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(f"ROI Validation Report for {sample_name}\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"Image Information:\n")
        f.write(f"  Dimensions: {image_shape[1]} × {image_shape[0]} pixels\n")
        f.write(f"  Total ROIs: {len(rois)}\n\n")
        
        # Validate coordinates
        validation_warnings = validate_roi_coordinates(rois, image_shape)
        
        if validation_warnings:
            f.write("Validation Warnings:\n")
            for warning in validation_warnings:
                f.write(f"  ⚠ {warning}\n")
            f.write("\n")
        else:
            f.write("✓ All ROI coordinates are valid\n\n")
        
        # ROI details
        f.write("ROI Details:\n")
        f.write("-" * 40 + "\n")
        
        for roi_name, coords in rois.items():
            width = coords['x_max'] - coords['x_min']
            height = coords['y_max'] - coords['y_min']
            area = width * height
            
            f.write(f"\n{roi_name}:\n")
            f.write(f"  X range: {coords['x_min']:.1f} - {coords['x_max']:.1f}\n")
            f.write(f"  Y range: {coords['y_min']:.1f} - {coords['y_max']:.1f}\n")
            f.write(f"  Size: {width:.1f} × {height:.1f} pixels\n")
            f.write(f"  Area: {area:.0f} pixels²\n")
    
    logging.info(f"Validation report saved to: {report_path}")


def main():
    """Main entry point for ROI validation."""
    parser = argparse.ArgumentParser(
        description="Validate ROI coordinates by overlaying them on source images",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python roi_validator.py --image_path tissue.tif --coords_file E14.5_ranges.txt --sample E14.5
  
  python roi_validator.py --image_path /data/images/tissue.tif \\
    --coords_file /data/coordinates/P3_ranges.txt --sample P3 \\
    --output_dir /data/validation/ --verbose
        """
    )
    
    parser.add_argument('--image_path', required=True,
                       help='Path to the source image file (TIFF, PNG, etc.)')
    parser.add_argument('--coords_file', required=True,
                       help='Path to coordinate range file')
    parser.add_argument('--sample', required=True,
                       help='Sample name for output file naming')
    parser.add_argument('--output_dir', default=None,
                       help='Output directory (default: same as coords_file directory)')
    parser.add_argument('--figsize', nargs=2, type=int, default=[12, 10],
                       help='Figure size in inches (width height)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='DPI for output image (default: 300)')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    try:
        # Set output directory
        if args.output_dir:
            output_dir = Path(args.output_dir)
        else:
            output_dir = Path(args.coords_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load image
        logging.info(f"Loading image: {args.image_path}")
        image = read_image(args.image_path)
        logging.info(f"Image loaded successfully: {image.shape}")
        
        # Parse coordinates
        logging.info(f"Parsing coordinates: {args.coords_file}")
        rois = parse_coordinate_file(args.coords_file)
        
        if not rois:
            logging.warning("No ROIs found in coordinate file")
            return 1
        
        # Validate coordinates
        validation_warnings = validate_roi_coordinates(rois, image.shape)
        if validation_warnings:
            logging.warning("Validation warnings found:")
            for warning in validation_warnings:
                logging.warning(f"  {warning}")
        
        # Generate validation plot
        output_image_path = output_dir / f"{args.sample}_roi_validation.png"
        plot_image_with_rois(
            image, rois, 
            output_path=str(output_image_path),
            figsize=tuple(args.figsize),
            dpi=args.dpi
        )
        
        # Generate validation report
        generate_roi_report(rois, image.shape, args.sample, str(output_dir))
        
        print("✓ ROI validation completed successfully!")
        print(f"✓ Validation plot: {output_image_path}")
        print(f"✓ Validation report: {output_dir}/{args.sample}_roi_validation_report.txt")
        
        if validation_warnings:
            print(f"⚠ {len(validation_warnings)} validation warnings found - check report")
        
    except Exception as e:
        logging.error(f"Validation failed: {e}")
        print(f"✗ ROI validation failed: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())# Version: 1.0.7
