"""
Workflows Module for Spatialcell Pipeline

This module provides high-level workflow orchestration for the complete Spatialcell
analysis pipeline, including single sample processing, batch operations, and
configuration management.
"""

from .complete_pipeline import run_complete_pipeline
from .batch_processor import (
    run_batch_processing, 
    BatchConfig,
    load_sample_configs,
    create_batch_configs_from_template,
    process_single_sample
)
from .pipeline_config import (
    PipelineConfig,
    validate_pipeline_inputs,
    create_example_config
)

__version__ = "0.1.0"
__author__ = "Xinyan"

# Define what gets imported with "from workflows import *"
__all__ = [
    # Complete pipeline
    'run_complete_pipeline',
    
    # Batch processing
    'run_batch_processing',
    'BatchConfig', 
    'load_sample_configs',
    'create_batch_configs_from_template',
    'process_single_sample',
    
    # Configuration management
    'PipelineConfig',
    'validate_pipeline_inputs',
    'create_example_config'
]# Version: 1.0.7
