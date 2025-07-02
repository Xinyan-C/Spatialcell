#!/usr/bin/env python3
"""
SVG to NPZ Converter for Spatialcell Pipeline

This script converts QuPath-exported SVG files containing nucleus detection results
into NPZ format for downstream spatial transcriptomics analysis.

Author: Xinyan

project: Spatialcell Pipeline

License: MIT
"""

import sys
import re
import xml.etree.ElementTree as ET
import numpy as np
import cv2
from scipy.sparse import csr_matrix, save_npz
import argparse
import logging
from pathlib import Path


def setup_logging(verbose=False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def parse_path_d(d_attr):
    """
    Extract all numbers from SVG path 'd' attribute using regular expressions,
    then pair them to return a list of points [(x1, y1), (x2, y2), ...].
    
    This method can handle cases where numbers lack clear separators.
    
    Args:
        d_attr (str): The 'd' attribute string from SVG path element
        
    Returns:
        list: List of coordinate tuples [(x1, y1), (x2, y2), ...]
    """
    # Regular expression to match integers or floats (allowing positive/negative signs)
    numbers = re.findall(r'[-+]?\d*\.\d+|[-+]?\d+', d_attr)
    
    try:
        numbers = [float(n) for n in numbers]
    except Exception as e:
        logging.error(f"Error converting numbers: {e}")
        return []
    
    if len(numbers) % 2 != 0:
        logging.warning("Extracted odd number of coordinates, possible format issue.")
        return []
    
    coords = []
    for i in range(0, len(numbers), 2):
        coords.append((numbers[i], numbers[i+1]))
    
    return coords


def svg_to_label_matrix(svg_path, image_height, image_width):
    """
    Convert QuPath-exported SVG file to label matrix.
    
    Prioritizes <path> elements for mask generation. If none found,
    tries <polygon> elements.
    
    Args:
        svg_path (str): Path to SVG file
        image_height (int): Original image height in pixels
        image_width (int): Original image width in pixels
        
    Returns:
        numpy.ndarray: 2D array (dtype=np.uint16) with pixel values 
                      corresponding to labels (1, 2, 3, ...)
                      
    Raises:
        SystemExit: If SVG parsing fails or no valid elements found
    """
    try:
        tree = ET.parse(svg_path)
    except Exception as e:
        logging.error(f"Error parsing SVG file: {e}")
        sys.exit(1)
    
    root = tree.getroot()
    
    # Auto-detect namespace
    ns = {}
    if root.tag.startswith("{"):
        uri, _ = root.tag[1:].split("}")
        ns = {'svg': uri}
        logging.info(f"Detected namespace: {ns}")
    else:
        logging.info("No namespace detected.")
    
    # Find <path> and <polygon> elements
    paths = root.findall('.//svg:path', ns) if ns else root.findall('.//path')
    polygons = root.findall('.//svg:polygon', ns) if ns else root.findall('.//polygon')
    
    logging.info(f"Found {len(paths)} path elements")
    logging.info(f"Found {len(polygons)} polygon elements")
    
    # Create empty mask (16-bit integer)
    mask = np.zeros((image_height, image_width), dtype=np.uint16)
    label = 1
    
    if len(paths) > 0:
        logging.info("Using <path> elements to generate mask")
        for path in paths:
            d_attr = path.attrib.get('d')
            if not d_attr:
                continue
            
            coords = parse_path_d(d_attr)
            if len(coords) < 3:
                continue  # Need at least three points to form a polygon
            
            pts = np.array(coords, dtype=np.int32).reshape((-1, 1, 2))
            cv2.fillPoly(mask, [pts], color=label)
            label += 1
            
    elif len(polygons) > 0:
        logging.info("Using <polygon> elements to generate mask")
        for poly in polygons:
            points_str = poly.attrib.get('points')
            if not points_str:
                continue
            
            points = []
            for point in points_str.strip().split():
                try:
                    x, y = map(float, point.split(','))
                    points.append([int(round(x)), int(round(y))])
                except Exception as e:
                    logging.error(f"Error parsing point '{point}': {e}")
                    continue
            
            if len(points) < 3:
                continue
            
            pts = np.array(points, dtype=np.int32).reshape((-1, 1, 2))
            cv2.fillPoly(mask, [pts], color=label)
            label += 1
    else:
        logging.error("No <path> or <polygon> elements found. Cannot determine SVG content.")
        sys.exit(1)
    
    return mask


def convert_svg_to_npz(svg_path, height, width, output_path, verbose=False):
    """
    Main conversion function.
    
    Args:
        svg_path (str): Path to input SVG file
        height (int): Image height in pixels
        width (int): Image width in pixels
        output_path (str): Path for output NPZ file
        verbose (bool): Enable verbose logging
        
    Returns:
        int: Number of detected objects
    """
    setup_logging(verbose)
    
    # Validate input files
    if not Path(svg_path).exists():
        raise FileNotFoundError(f"SVG file not found: {svg_path}")
    
    # Convert SVG to label matrix
    logging.info(f"Converting SVG file: {svg_path}")
    label_matrix = svg_to_label_matrix(svg_path, height, width)
    
    # Convert to sparse matrix and save as NPZ
    label_sparse = csr_matrix(label_matrix)
    save_npz(output_path, label_sparse)
    
    # Count unique objects
    num_objects = len(np.unique(label_matrix)) - 1  # Subtract 1 for background
    
    logging.info(f"Found {num_objects} nucleus objects")
    logging.info(f"Label matrix saved to {output_path}")
    
    return num_objects


def main():
    """Main entry point for command line usage."""
    parser = argparse.ArgumentParser(
        description="Convert QuPath-exported SVG files to NPZ format for spatial transcriptomics analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python svg_to_npz.py --svg nuclei.svg --height 2048 --width 2048 --output labels.npz
  python svg_to_npz.py --svg nuclei.svg --height 2048 --width 2048 --output labels.npz --verbose
        """
    )
    
    parser.add_argument('--svg', type=str, required=True,
                        help='Path to QuPath-exported SVG file')
    parser.add_argument('--height', type=int, required=True,
                        help='Original image height in pixels')
    parser.add_argument('--width', type=int, required=True,
                        help='Original image width in pixels')
    parser.add_argument('--output', type=str, required=True,
                        help='Output NPZ file path')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose logging')
    
    args = parser.parse_args()
    
    try:
        num_objects = convert_svg_to_npz(
            args.svg, args.height, args.width, args.output, args.verbose
        )
        print(f"✓ Successfully converted {num_objects} objects to {args.output}")
        
    except Exception as e:
        logging.error(f"Conversion failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()# Version: 1.0.7
