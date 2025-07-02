"""
Spatial Segmentation Module for Spatialcell Pipeline

This module integrates bin2cell functionality for comprehensive spatial transcriptomics analysis,
including cell segmentation, label expansion, and spatial visualization.

Built upon the bin2cell package: https://github.com/BayraktarLab/bin2cell
"""

from .spatial_processor import (
    process_spatial_data,
    read_regions,
    extract_roi_coordinates,
    visualize_labels_in_regions,
    visualize_after_insert_labels,
    visualize_after_expansion,
    visualize_gex_labels,
    visualize_joint_labels,
    generate_roi_data
)

__version__ = "0.1.0"
__author__ = "Xinyan"

# Define what gets imported with "from spatial_segmentation import *"
__all__ = [
    'process_spatial_data',
    'read_regions',
    'extract_roi_coordinates',
    'visualize_labels_in_regions',
    'visualize_after_insert_labels', 
    'visualize_after_expansion',
    'visualize_gex_labels',
    'visualize_joint_labels',
    'generate_roi_data'
]# Version: 1.0.7
