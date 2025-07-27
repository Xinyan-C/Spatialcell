#!/usr/bin/env python3
"""
TopAct Spatial Classification Visualizer for Spatialcell Pipeline

This module provides comprehensive visualization tools for spatial cell type classification
results from TopAct framework, supporting multiple output formats and ROI-based analysis.

Built upon TopAct framework: https://gitlab.com/kfbenjamin/topact.git

Author: Xinyan
License: Apache 2.0
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import load
import os
import sys
import re
from PIL import Image
import matplotlib.colors as mcolors
from natsort import natsorted
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import logging
import warnings

# Handle custom modules with fallback
try:
    from ..utils.color_schemes import get_cell_type_colors
    XINYAN_PALETTE_AVAILABLE = True
except ImportError:
    XINYAN_PALETTE_AVAILABLE = False
    warnings.warn("Color schemes module not found. Using matplotlib default colors.")

try:
    from ..utils import rename_utils  # Cell type renaming utilities
    RENAME_UTILS_AVAILABLE = True
except ImportError:
    RENAME_UTILS_AVAILABLE = False
    warnings.warn("rename_utils module not found. Cell type renaming disabled.")

# Avoid PIL decompression bomb warning for large images
Image.MAX_IMAGE_PIXELS = None


def setup_logging(verbose: bool = False) -> None:
    """
    Setup logging configuration.
    
    Args:
        verbose (bool): Enable verbose logging
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(levelname)s - %(message)s'
    )


def str2bool(v: Any) -> bool:
    """
    Convert string to boolean value.
    
    Args:
        v: Value to convert
        
    Returns:
        bool: Boolean value
        
    Raises:
        argparse.ArgumentTypeError: If value cannot be converted
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected')


def parse_roi_coordinates(roi_file: str) -> Dict[str, Tuple[float, float, float, float]]:
    """
    Parse ROI coordinates from coordinate file.
    
    Args:
        roi_file (str): Path to ROI coordinate file
        
    Returns:
        dict: ROI coordinates as {name: (x_min, x_max, y_min, y_max)}
        
    Raises:
        FileNotFoundError: If ROI file doesn't exist
        ValueError: If no ROIs found in file
    """
    if not Path(roi_file).exists():
        raise FileNotFoundError(f"ROI file not found: {roi_file}")
    
    with open(roi_file, 'r', encoding='utf-8') as f:
        content = f.read()
        
    # Parse coordinate patterns (handle multiple formats)
    pattern = r'([A-Za-z]+\d*(?:\.\d+)?\s*-\s*\w+\d*)\s*(?:Rectangle coordinates|矩形框的坐标范围):\s*X:\s*([\d\.\-]+)\s*-\s*([\d\.\-]+)\s*Y:\s*([\d\.\-]+)\s*-\s*([\d\.\-]+)'
    matches = re.findall(pattern, content)
    
    roi_dict = {}
    for match in matches:
        roi_name = match[0].strip()
        x_min, x_max, y_min, y_max = map(float, match[1:5])
        roi_dict[roi_name] = (x_min, x_max, y_min, y_max)
    
    if not roi_dict:
        raise ValueError(f"No ROIs found in {roi_file}")
    
    logging.info(f"Parsed {len(roi_dict)} ROI regions: {list(roi_dict.keys())}")
    return roi_dict


def hex_to_rgba(hex_color: str) -> np.ndarray:
    """
    Convert hex color to RGBA array.
    
    Args:
        hex_color (str): Hex color string (e.g., '#E31A1C')
        
    Returns:
        np.ndarray: RGBA color array
    """
    if hex_color.startswith('#'):
        hex_color = hex_color[1:]
    
    # Convert hex to RGB (0-1 range)
    rgb = tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    # Add alpha channel
    return np.array([rgb[0], rgb[1], rgb[2], 1.0])


def get_color_mapping(class_names: List[str], color_scheme: str = "primary") -> Dict[str, np.ndarray]:
    """
    Generate color mapping for cell types based on specified scheme.
    
    Args:
        class_names (list): List of cell type names
        color_scheme (str): Color scheme name
        
    Returns:
        dict: Mapping of class names to RGBA color arrays
    """
    try:
        # Import color schemes
        from ..utils.color_schemes import AVAILABLE_SCHEMES
        
        # Check if the specified color scheme exists
        if color_scheme in AVAILABLE_SCHEMES:
            scheme_colors = AVAILABLE_SCHEMES[color_scheme]
            
            # Create color mapping
            type_to_color = {}
            
            # First, map exact matches from the scheme
            for cell_type in class_names:
                if cell_type in scheme_colors:
                    hex_color = scheme_colors[cell_type]
                    type_to_color[cell_type] = hex_to_rgba(hex_color)
            
            # For unmapped cell types, generate additional colors
            unmapped_types = [ct for ct in class_names if ct not in type_to_color]
            if unmapped_types:
                logging.warning(f"Using fallback colors for {len(unmapped_types)} unmapped cell types: {unmapped_types}")
                # Use a colormap for unmapped types
                cmap = plt.cm.Set3
                for i, cell_type in enumerate(unmapped_types):
                    color = cmap(i / max(len(unmapped_types), 1))
                    type_to_color[cell_type] = np.array(color)
            
            logging.info(f"Color mapping generated using '{color_scheme}' scheme: {len(type_to_color)} colors assigned")
            return type_to_color
        
        else:
            logging.warning(f"Color scheme '{color_scheme}' not found in AVAILABLE_SCHEMES. Available schemes: {list(AVAILABLE_SCHEMES.keys())}")
            raise ValueError(f"Invalid color scheme: {color_scheme}")
            
    except ImportError as e:
        logging.warning(f"Failed to import color schemes: {e}. Using matplotlib fallback.")
    except Exception as e:
        logging.warning(f"Failed to load color scheme '{color_scheme}': {e}. Using matplotlib fallback.")
    
    # Fallback to matplotlib colormap
    logging.info(f"Using matplotlib fallback colormap for {len(class_names)} cell types")
    
    # Try to use the color_scheme as a matplotlib colormap
    try:
        cmap = plt.get_cmap(color_scheme)
        n_colors = min(len(class_names), cmap.N if hasattr(cmap, 'N') else 20)
        colors = [cmap(i / max(1, n_colors - 1)) for i in range(n_colors)]
        
        # Extend with default colors if needed
        if len(class_names) > len(colors):
            extra_colors = plt.cm.tab20(np.linspace(0, 1, len(class_names) - len(colors)))
            colors.extend(extra_colors)
            
        return dict(zip(class_names, [np.array(color) for color in colors[:len(class_names)]]))
        
    except Exception as e:
        logging.warning(f"Invalid matplotlib colormap '{color_scheme}': {e}. Using default tab20.")
        
        # Final fallback to tab20
        n_colors = min(len(class_names), 20)
        colors = plt.cm.tab20(np.linspace(0, 1, n_colors))
        if len(class_names) > 20:
            # Use tab20 + tab20b for more colors
            extra_colors = plt.cm.tab20b(np.linspace(0, 1, len(class_names) - 20))
            colors = np.concatenate([colors, extra_colors])
        
        return dict(zip(class_names, [np.array(color) for color in colors[:len(class_names)]]))


def apply_cell_type_renaming(class_names: List[str], sample: str, 
                           enable_renaming: bool = False) -> List[str]:
    """
    Apply cell type renaming if enabled and available.
    
    Args:
        class_names (list): Original class names
        sample (str): Sample identifier
        enable_renaming (bool): Whether to enable renaming
        
    Returns:
        list: Renamed or original class names
    """
    if not enable_renaming or not RENAME_UTILS_AVAILABLE:
        return class_names
    
    try:
        rename_dict = rename_utils.get_rename_dict(sample)
        if rename_dict:
            renamed_classes = [rename_dict.get(name, name) for name in class_names]
            logging.info(f"Applied cell type renaming for sample {sample}")
            logging.debug(f"Original: {class_names}")
            logging.debug(f"Renamed: {renamed_classes}")
            return renamed_classes
        else:
            logging.warning(f"No rename dictionary found for sample {sample}")
            
    except Exception as e:
        logging.warning(f"Failed to apply cell type renaming: {e}")
    
    return class_names


def load_classification_data(sd_path: str, outfile_path: str, clf_path: str) -> Tuple[Any, np.ndarray, List[str]]:
    """
    Load classification data and results.
    
    Args:
        sd_path (str): Path to CountGrid object
        outfile_path (str): Path to classification results
        clf_path (str): Path to classifier
        
    Returns:
        tuple: (CountGrid object, confidence matrix, class names)
        
    Raises:
        FileNotFoundError: If required files don't exist
    """
    # Validate file existence
    for path, name in [(sd_path, "CountGrid"), (outfile_path, "Classification results"), (clf_path, "Classifier")]:
        if not Path(path).exists():
            raise FileNotFoundError(f"{name} file not found: {path}")
    
    # Load data
    sd = load(sd_path)
    confidence_matrix = np.load(outfile_path)
    clf = load(clf_path)
    
    logging.info(f"Loaded classification data:")
    logging.info(f"  CountGrid: {len(sd.table)} data points")
    logging.info(f"  Confidence matrix shape: {confidence_matrix.shape}")
    logging.info(f"  Classifier classes: {len(clf.classes)}")
    
    return sd, confidence_matrix, list(clf.classes)


def process_classification_results(sd: Any, confidence_matrix: np.ndarray, 
                                 roi_coords: Tuple[float, float, float, float]) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """
    Process classification results for a specific ROI.
    
    Args:
        sd: CountGrid object
        confidence_matrix (np.ndarray): Classification confidence scores
        roi_coords (tuple): ROI coordinates (x_min, x_max, y_min, y_max)
        
    Returns:
        tuple: (coordinates DataFrame, classifications, max probabilities)
    """
    # Extract unique spatial coordinates
    unique_coords = sd.table[['x', 'y']].drop_duplicates().reset_index(drop=True)
    
    # Use maximum scale (last scale dimension)
    num_scales = confidence_matrix.shape[1]
    scale_idx = num_scales - 1
    probs = confidence_matrix[:, scale_idx, :]
    
    logging.info(f"Using scale index {scale_idx} of {num_scales} available scales")
    
    # Calculate classifications and max probabilities
    classifications = np.argmax(probs, axis=1)
    max_probs = np.max(probs, axis=1)
    
    # Mark unclassified points (all probabilities = -1)
    unclassified_mask = np.all(probs == -1, axis=1)
    classifications[unclassified_mask] = -1
    
    # Filter coordinates within ROI
    x_min, x_max, y_min, y_max = roi_coords
    mask = (
        (unique_coords['x'] >= x_min) & (unique_coords['x'] <= x_max) &
        (unique_coords['y'] >= y_min) & (unique_coords['y'] <= y_max)
    )
    
    coords_roi = unique_coords[mask].reset_index(drop=True)
    classifications_roi = classifications[mask]
    max_probs_roi = max_probs[mask]
    
    logging.info(f"ROI filtering: {len(coords_roi)} points within bounds")
    
    return coords_roi, classifications_roi, max_probs_roi


def create_scatter_plot(coords: pd.DataFrame, classifications: np.ndarray, 
                       class_names: List[str], type_to_color: Dict[str, Any],
                       roi_name: str, point_size: int = 10, point_shape: str = 's') -> plt.Figure:
    """
    Create scatter plot visualization of spatial classification.
    
    Args:
        coords (pd.DataFrame): Spatial coordinates
        classifications (np.ndarray): Classification labels
        class_names (list): Cell type names
        type_to_color (dict): Color mapping
        roi_name (str): ROI name for title
        point_size (int): Size of scatter points
        point_shape (str): Shape of scatter points
        
    Returns:
        plt.Figure: Generated figure
    """
    fig, ax = plt.subplots(figsize=(10, 14))
    ax.invert_yaxis()
    ax.set_aspect('equal', adjustable='box')
    
    # Get existing classes in this ROI
    existing_classes = np.unique(classifications)
    existing_classes = existing_classes[existing_classes != -1]
    
    # Create legend entries for existing classes only
    plot_class_names = [name for i, name in enumerate(class_names) if i in existing_classes]
    
    # Plot each cell type
    for i, name in enumerate(class_names):
        if i in existing_classes:
            mask = classifications == i
            color = type_to_color[name]
            ax.scatter(coords['x'][mask], coords['y'][mask],
                      c=[color], label=name, s=point_size, 
                      marker=point_shape, alpha=1.0)
    
    # Configure plot
    ax.grid(False)
    
    # Create legend
    handles = [plt.Line2D([0], [0], marker=point_shape, linestyle='None',
                         color=type_to_color[name], label=name, markersize=10)
              for name in plot_class_names]
    ax.legend(handles=handles, title='Cell type', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    logging.info(f"Created scatter plot for {roi_name} with {len(plot_class_names)} cell types")
    
    return fig


def create_overlay_plot(coords: pd.DataFrame, classifications: np.ndarray,
                       class_names: List[str], type_to_color: Dict[str, Any],
                       background_image: np.ndarray, roi_coords: Tuple[float, float, float, float],
                       roi_name: str, point_size: int = 10, point_shape: str = 's') -> plt.Figure:
    """
    Create overlay plot with background image and classification results.
    
    Args:
        coords (pd.DataFrame): Spatial coordinates
        classifications (np.ndarray): Classification labels
        class_names (list): Cell type names
        type_to_color (dict): Color mapping
        background_image (np.ndarray): Background image
        roi_coords (tuple): ROI coordinates for cropping
        roi_name (str): ROI name
        point_size (int): Size of scatter points
        point_shape (str): Shape of scatter points
        
    Returns:
        plt.Figure: Generated figure
    """
    x_min, x_max, y_min, y_max = roi_coords
    
    fig, ax = plt.subplots(figsize=(10, 14))
    ax.invert_yaxis()
    ax.set_aspect('equal', adjustable='box')
    
    # Display background image
    ax.imshow(background_image, extent=[x_min, x_max, y_max, y_min], alpha=1.0)
    
    # Get existing classes
    existing_classes = np.unique(classifications)
    existing_classes = existing_classes[existing_classes != -1]
    
    plot_class_names = [name for i, name in enumerate(class_names) if i in existing_classes]
    
    # Plot classifications
    for i, name in enumerate(class_names):
        if i in existing_classes:
            mask = classifications == i
            color = type_to_color[name]
            ax.scatter(coords['x'][mask], coords['y'][mask],
                      c=[color], label=name, s=point_size,
                      marker=point_shape, alpha=1.0)
    
    # Configure plot
    ax.grid(False)
    
    # Create legend
    handles = [plt.Line2D([0], [0], marker=point_shape, linestyle='None',
                         color=type_to_color[name], label=name, markersize=10)
              for name in plot_class_names]
    ax.legend(handles=handles, title='Cell type', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    logging.info(f"Created overlay plot for {roi_name}")
    
    return fig


def create_side_by_side_plot(coords: pd.DataFrame, classifications: np.ndarray,
                           class_names: List[str], type_to_color: Dict[str, Any],
                           background_image: np.ndarray, roi_coords: Tuple[float, float, float, float],
                           roi_name: str, point_size: int = 10, point_shape: str = 's') -> plt.Figure:
    """
    Create side-by-side plot with background image and classification results.
    
    Args:
        coords (pd.DataFrame): Spatial coordinates
        classifications (np.ndarray): Classification labels
        class_names (list): Cell type names
        type_to_color (dict): Color mapping
        background_image (np.ndarray): Background image
        roi_coords (tuple): ROI coordinates
        roi_name (str): ROI name
        point_size (int): Size of scatter points
        point_shape (str): Shape of scatter points
        
    Returns:
        plt.Figure: Generated figure
    """
    x_min, x_max, y_min, y_max = roi_coords
    
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(20, 14), sharey=True)
    
    # Left panel: Background image
    ax_left.imshow(background_image, extent=[x_min, x_max, y_max, y_min])
    ax_left.invert_yaxis()
    ax_left.set_aspect('equal', adjustable='box')
    ax_left.grid(False)
    
    # Right panel: Classification results
    ax_right.invert_yaxis()
    ax_right.set_aspect('equal', adjustable='box')
    
    # Get existing classes
    existing_classes = np.unique(classifications)
    existing_classes = existing_classes[existing_classes != -1]
    
    plot_class_names = [name for i, name in enumerate(class_names) if i in existing_classes]
    
    # Plot classifications
    for i, name in enumerate(class_names):
        if i in existing_classes:
            mask = classifications == i
            color = type_to_color[name]
            ax_right.scatter(coords['x'][mask], coords['y'][mask],
                           c=[color], label=name, s=point_size,
                           marker=point_shape, alpha=1.0)
    
    ax_right.grid(False)
    
    # Create legend
    handles = [plt.Line2D([0], [0], marker=point_shape, linestyle='None',
                         color=type_to_color[name], label=name, markersize=10)
              for name in plot_class_names]
    ax_right.legend(handles=handles, title='Cell type', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    logging.info(f"Created side-by-side plot for {roi_name}")
    
    return fig



def process_roi_visualization(roi_name: str, roi_coords: Tuple[float, float, float, float],
                            args: argparse.Namespace, class_names: List[str],
                            renamed_class_names: List[str], type_to_color: Dict[str, Any]) -> None:
    """
    Process visualization for a single ROI.
    
    Args:
        roi_name (str): ROI identifier
        roi_coords (tuple): ROI coordinates
        args: Command line arguments
        class_names (list): Original class names
        renamed_class_names (list): Renamed class names
        type_to_color (dict): Color mapping
    """
    logging.info(f"Processing ROI: {roi_name}")
    
    # Construct file paths
    safe_roi_name = roi_name.replace(' ', '_')
    sd_filename = f"sd_{args.sample}_{safe_roi_name}.joblib"
    outfile_filename = f"outfile_{args.sample}_{safe_roi_name}.npy"
    
    sd_path = os.path.join(args.sd_dir, sd_filename)
    outfile_path = os.path.join(args.outfile_dir, outfile_filename)
    
    # Check if files exist
    try:
        sd, confidence_matrix, _ = load_classification_data(sd_path, outfile_path, args.clf_path)
    except FileNotFoundError as e:
        logging.warning(f"Skipping {roi_name}: {e}")
        return
    
    # Process classification results
    coords_roi, classifications_roi, max_probs_roi = process_classification_results(
        sd, confidence_matrix, roi_coords
    )
    
    if len(coords_roi) == 0:
        logging.warning(f"No data points found in ROI {roi_name}")
        return
    
    # Create scatter plot
    fig1 = create_scatter_plot(
        coords_roi, classifications_roi, renamed_class_names, type_to_color,
        roi_name, args.point_size, args.point_shape
    )
    
    scatter_path = os.path.join(args.output_dir, f"Spatial_Classification_{safe_roi_name}.pdf")
    fig1.savefig(scatter_path, dpi=600, bbox_inches='tight')
    plt.close(fig1)
    logging.info(f"Scatter plot saved: {scatter_path}")
    
    # Process background image if provided
    if args.background_image and Path(args.background_image).exists():
        try:
            img = Image.open(args.background_image)
            x_min, x_max, y_min, y_max = roi_coords
            img_roi = img.crop((x_min, y_min, x_max, y_max))
            img_array = np.array(img_roi)
            
            # Create overlay plot
            fig2 = create_overlay_plot(
                coords_roi, classifications_roi, renamed_class_names, type_to_color,
                img_array, roi_coords, roi_name, args.point_size, args.point_shape
            )
            
            overlay_path = os.path.join(args.output_dir, f"Spatial_Classification_{safe_roi_name}_overlay.pdf")
            fig2.savefig(overlay_path, dpi=600, bbox_inches='tight')
            plt.close(fig2)
            logging.info(f"Overlay plot saved: {overlay_path}")
            
            # Create side-by-side plot
            fig3 = create_side_by_side_plot(
                coords_roi, classifications_roi, renamed_class_names, type_to_color,
                img_array, roi_coords, roi_name, args.point_size, args.point_shape
            )
            
            side_by_side_path = os.path.join(args.output_dir, f"Spatial_Classification_{safe_roi_name}_side_by_side.pdf")
            fig3.savefig(side_by_side_path, dpi=600, bbox_inches='tight')
            plt.close(fig3)
            logging.info(f"Side-by-side plot saved: {side_by_side_path}")
            
        except Exception as e:
            logging.warning(f"Failed to process background image for {roi_name}: {e}")
    



def save_figures_with_auto_mkdir(figures: List[plt.Figure], 
                                output_dir: str, 
                                filename_prefix: str = "roi",
                                file_format: str = "pdf",
                                dpi: int = 600) -> List[str]:
    """
    Save figures with automatic directory creation.
    
    Args:
        figures: List of matplotlib figures
        output_dir: Output directory path
        filename_prefix: Prefix for filenames
        file_format: File format (pdf, png, svg, etc.)
        dpi: Resolution for raster formats
        
    Returns:
        List of saved file paths
    """
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    saved_files = []
    
    for i, fig in enumerate(figures):
        filename = f"{filename_prefix}_{i+1}.{file_format}"
        filepath = output_path / filename
        
        try:
            fig.savefig(str(filepath), dpi=dpi, bbox_inches='tight')
            saved_files.append(str(filepath))
            logging.info(f"Saved figure: {filepath}")
        except Exception as e:
            logging.error(f"Failed to save figure {filename}: {e}")
    
    return saved_files


def visualize_for_notebook(sample: str, 
                          sd_dir: str, 
                          outfile_dir: str, 
                          clf_dir: str, 
                          roi_file: str,
                          background_image: Optional[str] = None,
                          point_size: int = 10,
                          point_shape: str = 's',
                          color_scheme: str = "primary",
                          rename_cell_types: bool = False,
                          show_legend: bool = True) -> List[plt.Figure]:
    """
    Create visualizations optimized for Jupyter notebook display.
    
    Args:
        sample (str): Sample identifier
        sd_dir (str): Directory containing CountGrid files
        outfile_dir (str): Directory containing classification results
        clf_dir (str): Directory containing classifier files
        roi_file (str): Path to ROI coordinates file
        background_image (str, optional): Path to background tissue image
        point_size (int): Size of scatter points
        point_shape (str): Shape of scatter points
        color_scheme (str): Color scheme to use
        rename_cell_types (bool): Enable cell type renaming
        show_legend (bool): Whether to show legend
        
    Returns:
        List[plt.Figure]: List of matplotlib figures for display
    """
    setup_logging(False)  # Disable verbose logging for notebook
    
    try:
        # Parse ROI coordinates
        roi_dict = parse_roi_coordinates(roi_file)
        
        # Load classifier to get class names
        clf_path = os.path.join(clf_dir, f"clf_{sample}.joblib")
        if not Path(clf_path).exists():
            raise FileNotFoundError(f"Classifier not found: {clf_path}")
        
        clf = load(clf_path)
        class_names = list(clf.classes)
        
        # Apply cell type renaming if enabled
        renamed_class_names = apply_cell_type_renaming(class_names, sample, rename_cell_types)
        
        # Generate color mapping with fallback
        try:
            type_to_color = get_color_mapping(renamed_class_names, color_scheme)
        except Exception as e:
            logging.warning(f"Failed to load color scheme '{color_scheme}': {e}")
            logging.warning("Falling back to 'primary' color scheme")
            type_to_color = get_color_mapping(renamed_class_names, "primary")
        
        figures = []
        
        # Process each ROI
        for roi_name, roi_coords in roi_dict.items():
            try:
                logging.info(f"Processing ROI: {roi_name}")
                
                # Construct file paths
                safe_roi_name = roi_name.replace(' ', '_')
                sd_filename = f"sd_{sample}_{safe_roi_name}.joblib"
                outfile_filename = f"outfile_{sample}_{safe_roi_name}.npy"
                
                sd_path = os.path.join(sd_dir, sd_filename)
                outfile_path = os.path.join(outfile_dir, outfile_filename)
                
                # Load classification data
                sd, confidence_matrix, _ = load_classification_data(sd_path, outfile_path, clf_path)
                
                # Process classification results
                coords_roi, classifications_roi, max_probs_roi = process_classification_results(
                    sd, confidence_matrix, roi_coords
                )
                
                if len(coords_roi) == 0:
                    logging.warning(f"No data points found in ROI {roi_name}")
                    continue
                
                # Create scatter plot
                fig = create_scatter_plot(
                    coords_roi, classifications_roi, renamed_class_names, type_to_color,
                    roi_name, point_size, point_shape
                )
                
                # Adjust legend if requested
                if not show_legend:
                    if fig.axes:
                        legend = fig.axes[0].get_legend()
                        if legend:
                            legend.remove()
                
                figures.append(fig)
                
                # Create background overlay if image provided
                if background_image and Path(background_image).exists():
                    try:
                        img = Image.open(background_image)
                        x_min, x_max, y_min, y_max = roi_coords
                        img_roi = img.crop((x_min, y_min, x_max, y_max))
                        img_array = np.array(img_roi)
                        
                        # Create overlay plot
                        fig_overlay = create_overlay_plot(
                            coords_roi, classifications_roi, renamed_class_names, type_to_color,
                            img_array, roi_coords, roi_name, point_size, point_shape
                        )
                        
                        if not show_legend:
                            if fig_overlay.axes:
                                legend = fig_overlay.axes[0].get_legend()
                                if legend:
                                    legend.remove()
                        
                        figures.append(fig_overlay)
                        
                    except Exception as e:
                        logging.warning(f"Failed to create overlay for {roi_name}: {e}")
                
            except Exception as e:
                logging.error(f"Failed to process ROI {roi_name}: {e}")
                continue
        
        logging.info(f"Generated {len(figures)} visualization figures")
        return figures
        
    except Exception as e:
        logging.error(f"Visualization failed: {e}")
        raise


def display_roi_visualizations(sample: str,
                              sd_dir: str,
                              outfile_dir: str, 
                              clf_dir: str,
                              roi_file: str,
                              **kwargs) -> None:
    """
    Convenience function to display all ROI visualizations in notebook.
    
    Args:
        sample (str): Sample identifier
        sd_dir (str): Directory containing CountGrid files
        outfile_dir (str): Directory containing classification results
        clf_dir (str): Directory containing classifier files
        roi_file (str): Path to ROI coordinates file
        **kwargs: Additional arguments for visualize_for_notebook
    """
    figures = visualize_for_notebook(
        sample=sample,
        sd_dir=sd_dir,
        outfile_dir=outfile_dir,
        clf_dir=clf_dir,
        roi_file=roi_file,
        **kwargs
    )
    
    for i, fig in enumerate(figures):
        plt.figure(fig.number)
        plt.show()

def main():
    """Main entry point for spatial classification visualization."""
    parser = argparse.ArgumentParser(
        description="Visualize spatial cell type classification results from TopAct framework",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic visualization
  python annotation_visualizer.py --sample E14.5 \\
    --sd_dir /data/results/ --outfile_dir /data/results/ \\
    --clf_dir /data/classifiers/ --roi_file /data/E14.5_ranges.txt

  # With background image and custom colors
  python annotation_visualizer.py --sample P3 \\
    --sd_dir /data/results/ --outfile_dir /data/results/ \\
    --clf_dir /data/classifiers/ --roi_file /data/P3_ranges.txt \\
    --background_image /data/tissue.tif --color_scheme Set3 \\
    --rename_cell_types True --verbose

This tool generates comprehensive visualizations including:
- Scatter plots of spatial classification
- Background image overlays
- Side-by-side comparisons

        """
    )
    
    # Required parameters
    parser.add_argument("--sample", type=str, required=True,
                        help="Sample identifier (e.g., E14.5, E18.5, P3)")
    parser.add_argument("--sd_dir", type=str, required=True,
                        help="Directory containing CountGrid (sd) files")
    parser.add_argument("--outfile_dir", type=str, required=True,
                        help="Directory containing classification result (outfile) files")
    parser.add_argument("--clf_dir", type=str, required=True,
                        help="Directory containing trained classifier files")
    parser.add_argument("--roi_file", type=str, required=True,
                        help="Path to ROI coordinates file")
    
    # Optional parameters
    parser.add_argument("--output_dir", type=str, default=None,
                        help="Output directory for visualizations (default: current directory)")
    parser.add_argument("--background_image", type=str, default=None,
                        help="Path to background tissue image (optional)")
    parser.add_argument("--point_size", type=int, default=10,
                        help="Size of scatter plot points (default: 10)")
    parser.add_argument("--point_shape", type=str, default='s',
                        help="Shape of scatter plot points (default: 's')")
    parser.add_argument("--color_scheme", type=str, default="primary",
                        help="Color scheme (primary, scientific, functional, modern, warm, golden, article_reproducibility) (default: primary)")
    parser.add_argument("--rename_cell_types", type=str2bool, default=False,
                        help="Enable cell type renaming (default: False)")
    parser.add_argument("--verbose", action='store_true',
                        help="Enable verbose logging")
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    try:
        # Set output directory
        if args.output_dir:
            args.output_dir = Path(args.output_dir)
            args.output_dir.mkdir(parents=True, exist_ok=True)
        else:
            args.output_dir = Path.cwd()
        
        # Parse ROI coordinates
        roi_dict = parse_roi_coordinates(args.roi_file)
        
        # Load classifier to get class names
        args.clf_path = os.path.join(args.clf_dir, f"clf_{args.sample}.joblib")
        if not Path(args.clf_path).exists():
            raise FileNotFoundError(f"Classifier not found: {args.clf_path}")
        
        clf = load(args.clf_path)
        class_names = list(clf.classes)
        logging.info(f"Loaded classifier with {len(class_names)} cell types")
        
        # Apply cell type renaming if enabled
        renamed_class_names = apply_cell_type_renaming(class_names, args.sample, args.rename_cell_types)
        
        # Sort class names naturally and handle duplicates
        class_name_pairs = [(i, name) for i, name in enumerate(renamed_class_names)]
        class_name_pairs = natsorted(class_name_pairs, key=lambda x: x[1])
        
        # Remove duplicates while preserving all original indices
        unique_class_names = []
        name_to_indices = {}
        for class_idx, name in class_name_pairs:
            if name not in name_to_indices:
                name_to_indices[name] = []
                unique_class_names.append(name)
            name_to_indices[name].append(class_idx)
        
        # Generate color mapping
        type_to_color = get_color_mapping(unique_class_names, args.color_scheme)
        logging.info(f"Generated color mapping for {len(unique_class_names)} unique cell types")
        
        # Process each ROI
        for roi_name, roi_coords in roi_dict.items():
            try:
                process_roi_visualization(roi_name, roi_coords, args, class_names, 
                                        renamed_class_names, type_to_color)
            except Exception as e:
                logging.error(f"Failed to process ROI {roi_name}: {e}")
                continue
        
        print("Spatial classification visualization completed successfully!")
        print(f"Results saved to: {args.output_dir}")
        
    except Exception as e:
        logging.error(f"Visualization failed: {e}")
        print(f"Visualization failed: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())