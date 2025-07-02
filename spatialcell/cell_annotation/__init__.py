"""
Cell Annotation Module for Spatialcell Pipeline

This module provides comprehensive tools for spatial cell type annotation using TopAct framework,
including classifier training, spatial processing, and visualization capabilities.

Built upon TopAct framework: https://github.com/carmonalab/TopAct
"""

# Classifier training utilities
from .classifier_trainer import (
    train_topact_classifier,
    load_reference_data,
    subset_by_time_point,
    train_multiple_time_points,
    save_classifier
)

# Spatial annotation processing
from .annotation_processor import (
    process_sample_annotation,
    parse_roi_coordinates,
    load_hd_spatial_data,
    process_hd_data_with_cell_labels,
    build_spatial_countgrid,
    perform_multiscale_classification
)

# Visualization tools  
from .annotation_visualizer import (
    parse_roi_coordinates as parse_roi_coords_viz,
    load_classification_data,
    process_classification_results,
    create_scatter_plot,
    create_overlay_plot,
    save_classification_csv,
    apply_cell_type_renaming
)

__version__ = "0.1.0"
__author__ = "Xinyan"

# Define what gets imported with "from cell_annotation import *"
__all__ = [
    # Classifier training
    'train_topact_classifier',
    'load_reference_data', 
    'subset_by_time_point',
    'train_multiple_time_points',
    'save_classifier',
    
    # Spatial processing
    'process_sample_annotation',
    'parse_roi_coordinates',
    'load_hd_spatial_data',
    'process_hd_data_with_cell_labels',
    'build_spatial_countgrid', 
    'perform_multiscale_classification',
    
    # Visualization
    'parse_roi_coords_viz',
    'load_classification_data',
    'process_classification_results',
    'create_scatter_plot',
    'create_overlay_plot', 
    'save_classification_csv',
    'apply_cell_type_renaming'
]# Version: 1.0.7
