"""
Utilities Module for Spatialcell Pipeline

This module provides utility tools for ROI coordinate extraction, validation,
spatial analysis, and color scheme management for spatial transcriptomics analysis.
"""

# ROI processing utilities
from .roi_extractor import extract_roi_coordinates, calculate_roi_range, detect_roi_files
from .roi_validator import parse_coordinate_file, plot_image_with_rois, validate_roi_coordinates

# Enhanced spatial analysis (TopAct extensions)
from .spatial_analysis import (
    EnhancedCountGrid, FloatExpressionGrid, SpatialClassificationWorker,
    combine_coordinates, split_coordinates, get_x_coordinate, get_y_coordinate,
    extract_classifications_from_confidence, extract_classification_image
)

# Color schemes for visualization
from .color_schemes import (
    get_cell_type_colors, get_color_families, list_available_schemes,
    get_scheme_description, validate_color_scheme, PRIMARY_CELL_TYPE_COLORS,
    AVAILABLE_SCHEMES
)

__version__ = "0.1.0"
__author__ = "Xinyan"

# Define what gets imported with "from utils import *"
__all__ = [
    # ROI utilities
    'extract_roi_coordinates', 'calculate_roi_range', 'detect_roi_files',
    'parse_coordinate_file', 'plot_image_with_rois', 'validate_roi_coordinates',
    
    # Spatial analysis
    'EnhancedCountGrid', 'FloatExpressionGrid', 'SpatialClassificationWorker',
    'combine_coordinates', 'split_coordinates', 'get_x_coordinate', 'get_y_coordinate',
    'extract_classifications_from_confidence', 'extract_classification_image',
    
    # Color schemes
    'get_cell_type_colors', 'get_color_families', 'list_available_schemes',
    'get_scheme_description', 'validate_color_scheme', 'PRIMARY_CELL_TYPE_COLORS',
    'AVAILABLE_SCHEMES'
]# Version: 1.0.7
