#!/usr/bin/env python3
"""
Pipeline Configuration Management for Spatialcell

This module provides comprehensive configuration management for the Spatialcell pipeline,
supporting YAML and JSON configuration files with validation and default value handling.

Author: Xinyan
License: MIT
"""

import os
import json
import yaml
from pathlib import Path
from typing import Dict, Any, Optional, List, Union
from dataclasses import dataclass, field, asdict
import logging


@dataclass
class PipelineConfig:
    """
    Comprehensive configuration class for Spatialcell pipeline.
    
    Contains all parameters needed for end-to-end pipeline execution,
    with sensible defaults and validation capabilities.
    """
    
    # =========================================================================
    # Basic Sample Information
    # =========================================================================
    sample_name: str = "sample"
    output_dir: str = "./spatialcell_output"
    
    # =========================================================================
    # Pipeline Step Controls
    # =========================================================================
    enable_roi_extraction: bool = True
    enable_svg_conversion: bool = True
    enable_spatial_segmentation: bool = True
    enable_cell_annotation: bool = True
    enable_visualization: bool = True
    
    # =========================================================================
    # ROI Extraction Parameters
    # =========================================================================
    loupe_export_dir: Optional[str] = None
    roi_coordinate_file: Optional[str] = None
    
    # =========================================================================
    # SVG Conversion Parameters
    # =========================================================================
    svg_file: Optional[str] = None
    image_height: int = 2048
    image_width: int = 2048
    npz_labels_file: Optional[str] = None
    
    # =========================================================================
    # Spatial Segmentation Parameters
    # =========================================================================
    visium_data_path: Optional[str] = None
    source_image_path: Optional[str] = None
    
    # StarDist parameters
    stardist_prob_thresh: float = 0.05
    stardist_nms_thresh: float = 0.5
    
    # Label expansion parameters
    labels_key: str = "labels_joint"
    expansion_algorithm: str = "max_bin_distance"
    max_bin_distance: int = 2
    volume_ratio: float = 4.0
    nearest_neighbors: int = 4
    subset_pca: bool = True
    
    # =========================================================================
    # Cell Annotation Parameters
    # =========================================================================
    classifier_path: Optional[str] = None
    hd_expression_path: Optional[str] = None
    hd_positions_path: Optional[str] = None
    cell_labels_column: str = "labels_qupath_expanded"
    
    # Multi-scale classification parameters
    min_scale: float = 3.0
    max_scale: float = 9.0
    
    # =========================================================================
    # Visualization Parameters
    # =========================================================================
    background_image_path: Optional[str] = None
    color_scheme: str = "primary"
    point_size: int = 10
    point_shape: str = "s"
    rename_cell_types: bool = False
    
    # =========================================================================
    # Computational Resources
    # =========================================================================
    num_processes: int = 80
    num_threads: int = 80
    memory_limit_gb: int = 600
    
    # =========================================================================
    # General Settings
    # =========================================================================
    verbose: bool = False
    random_seed: Optional[int] = None
    
    # =========================================================================
    # Derived Properties (set automatically)
    # =========================================================================
    _derived_paths: Dict[str, str] = field(default_factory=dict, init=False)
    
    def __post_init__(self):
        """Post-initialization setup and path derivation."""
        self._setup_derived_paths()
        self._validate_basic_config()
    
    def _setup_derived_paths(self) -> None:
        """Setup derived file paths based on sample name and output directory."""
        base_output = Path(self.output_dir)
        
        # Set default paths if not provided
        if not self.roi_coordinate_file:
            self.roi_coordinate_file = str(base_output / f"{self.sample_name}_ranges.txt")
        
        if not self.npz_labels_file:
            self.npz_labels_file = str(base_output / "preprocessing" / f"{self.sample_name}_labels.npz")
        
        # Store derived paths for reference
        self._derived_paths = {
            "base_output": str(base_output),
            "preprocessing_dir": str(base_output / "preprocessing"),
            "segmentation_dir": str(base_output / "spatial_segmentation"),
            "annotation_dir": str(base_output / "cell_annotation"),
            "visualization_dir": str(base_output / "visualizations"),
            "logs_dir": str(base_output / "logs"),
            "reports_dir": str(base_output / "reports")
        }
    
    def _validate_basic_config(self) -> None:
        """Validate basic configuration parameters."""
        if not self.sample_name or not self.sample_name.strip():
            raise ValueError("Sample name cannot be empty")
        
        if self.image_height <= 0 or self.image_width <= 0:
            raise ValueError("Image dimensions must be positive")
        
        if self.min_scale >= self.max_scale:
            raise ValueError("min_scale must be less than max_scale")
        
        if self.num_processes <= 0 or self.num_threads <= 0:
            raise ValueError("Number of processes and threads must be positive")
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'PipelineConfig':
        """
        Create configuration from dictionary.
        
        Args:
            config_dict (dict): Configuration parameters
            
        Returns:
            PipelineConfig: Configured pipeline object
        """
        # Filter out unknown parameters
        valid_fields = {f.name for f in cls.__dataclass_fields__.values() if f.init}
        filtered_dict = {k: v for k, v in config_dict.items() if k in valid_fields}
        
        return cls(**filtered_dict)
    
    @classmethod
    def from_file(cls, config_path: str) -> 'PipelineConfig':
        """
        Load configuration from YAML or JSON file.
        
        Args:
            config_path (str): Path to configuration file
            
        Returns:
            PipelineConfig: Loaded configuration
            
        Raises:
            FileNotFoundError: If config file doesn't exist
            ValueError: If file format is unsupported or invalid
        """
        config_path = Path(config_path)
        
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                if config_path.suffix.lower() in ['.yaml', '.yml']:
                    config_dict = yaml.safe_load(f)
                elif config_path.suffix.lower() == '.json':
                    config_dict = json.load(f)
                else:
                    raise ValueError(f"Unsupported config file format: {config_path.suffix}")
            
            return cls.from_dict(config_dict)
            
        except (yaml.YAMLError, json.JSONDecodeError) as e:
            raise ValueError(f"Invalid configuration file format: {e}")
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert configuration to dictionary.
        
        Returns:
            dict: Configuration as dictionary
        """
        config_dict = asdict(self)
        # Remove private fields
        config_dict = {k: v for k, v in config_dict.items() if not k.startswith('_')}
        return config_dict
    
    def to_yaml(self) -> str:
        """
        Convert configuration to YAML string.
        
        Returns:
            str: Configuration as YAML
        """
        return yaml.dump(self.to_dict(), default_flow_style=False, sort_keys=False)
    
    def to_json(self) -> str:
        """
        Convert configuration to JSON string.
        
        Returns:
            str: Configuration as JSON
        """
        return json.dumps(self.to_dict(), indent=2)
    
    def save(self, output_path: str, format: str = 'yaml') -> None:
        """
        Save configuration to file.
        
        Args:
            output_path (str): Output file path
            format (str): Output format ('yaml' or 'json')
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            if format.lower() == 'yaml':
                f.write(self.to_yaml())
            elif format.lower() == 'json':
                f.write(self.to_json())
            else:
                raise ValueError(f"Unsupported format: {format}")
        
        logging.info(f"Configuration saved to: {output_path}")
    
    def get_derived_path(self, path_key: str) -> str:
        """
        Get derived output path.
        
        Args:
            path_key (str): Path key ('base_output', 'segmentation_dir', etc.)
            
        Returns:
            str: Derived path
            
        Raises:
            KeyError: If path key not found
        """
        if path_key not in self._derived_paths:
            available_keys = list(self._derived_paths.keys())
            raise KeyError(f"Path key '{path_key}' not found. Available: {available_keys}")
        
        return self._derived_paths[path_key]
    
    def create_output_directories(self) -> None:
        """Create all output directories."""
        for dir_path in self._derived_paths.values():
            Path(dir_path).mkdir(parents=True, exist_ok=True)
        
        logging.info("Created output directory structure")
    
    def validate_file_paths(self) -> List[str]:
        """
        Validate that required input files exist.
        
        Returns:
            list: List of validation error messages
        """
        errors = []
        
        # Check required files based on enabled steps
        if self.enable_roi_extraction and self.loupe_export_dir:
            if not Path(self.loupe_export_dir).exists():
                errors.append(f"Loupe export directory not found: {self.loupe_export_dir}")
        
        if self.enable_svg_conversion and self.svg_file:
            if not Path(self.svg_file).exists():
                errors.append(f"SVG file not found: {self.svg_file}")
        
        if self.enable_spatial_segmentation:
            if self.visium_data_path and not Path(self.visium_data_path).exists():
                errors.append(f"Visium data path not found: {self.visium_data_path}")
            
            if self.source_image_path and not Path(self.source_image_path).exists():
                errors.append(f"Source image not found: {self.source_image_path}")
        
        if self.enable_cell_annotation:
            if self.classifier_path and not Path(self.classifier_path).exists():
                errors.append(f"Classifier file not found: {self.classifier_path}")
            
            if self.hd_expression_path and not Path(self.hd_expression_path).exists():
                errors.append(f"HD expression data not found: {self.hd_expression_path}")
            
            if self.hd_positions_path and not Path(self.hd_positions_path).exists():
                errors.append(f"HD positions file not found: {self.hd_positions_path}")
        
        if self.enable_visualization and self.background_image_path:
            if not Path(self.background_image_path).exists():
                errors.append(f"Background image not found: {self.background_image_path}")
        
        return errors


def validate_pipeline_inputs(config: PipelineConfig) -> None:
    """
    Comprehensive validation of pipeline inputs.
    
    Args:
        config (PipelineConfig): Pipeline configuration to validate
        
    Raises:
        ValueError: If validation fails
    """
    logging.info("Validating pipeline configuration...")
    
    # Validate file paths
    file_errors = config.validate_file_paths()
    if file_errors:
        error_msg = "Input validation failed:\n" + "\n".join(f"  - {error}" for error in file_errors)
        raise ValueError(error_msg)
    
    # Validate parameter ranges
    if config.stardist_prob_thresh < 0 or config.stardist_prob_thresh > 1:
        raise ValueError("StarDist probability threshold must be between 0 and 1")
    
    if config.stardist_nms_thresh < 0 or config.stardist_nms_thresh > 1:
        raise ValueError("StarDist NMS threshold must be between 0 and 1")
    
    if config.volume_ratio <= 0:
        raise ValueError("Volume ratio must be positive")
    
    if config.nearest_neighbors <= 0:
        raise ValueError("Number of nearest neighbors must be positive")
    
    if config.point_size <= 0:
        raise ValueError("Point size must be positive")
    
    # Validate color scheme
    valid_color_schemes = ['primary', 'scientific', 'functional', 'modern', 'warm', 'golden']
    if config.color_scheme not in valid_color_schemes:
        raise ValueError(f"Invalid color scheme. Valid options: {valid_color_schemes}")
    
    # Validate point shape
    valid_point_shapes = ['o', 's', '^', 'v', '<', '>', 'D', 'd', 'p', '*', 'h', 'H', '+', 'x']
    if config.point_shape not in valid_point_shapes:
        raise ValueError(f"Invalid point shape. Valid options: {valid_point_shapes}")
    
    # Validate expansion algorithm
    valid_algorithms = ['max_bin_distance', 'volume_ratio']
    if config.expansion_algorithm not in valid_algorithms:
        raise ValueError(f"Invalid expansion algorithm. Valid options: {valid_algorithms}")
    
    # Validate labels key
    valid_labels = ['labels_qupath', 'labels_qupath_expanded', 'labels_gex', 'labels_joint']
    if config.labels_key not in valid_labels:
        raise ValueError(f"Invalid labels key. Valid options: {valid_labels}")
    
    logging.info("Configuration validation passed")


def create_example_config(output_path: str, sample_name: str = "example_sample") -> None:
    """
    Create an example configuration file with documentation.
    
    Args:
        output_path (str): Path to save example config
        sample_name (str): Sample name for example
    """
    config = PipelineConfig(
        sample_name=sample_name,
        output_dir=f"./spatialcell_output_{sample_name}",
        
        # Enable all pipeline steps
        enable_roi_extraction=True,
        enable_svg_conversion=True,
        enable_spatial_segmentation=True,
        enable_cell_annotation=True,
        enable_visualization=True,
        
        # ROI extraction (adjust paths as needed)
        loupe_export_dir=f"./data/{sample_name}/loupe_exports/",
        
        # SVG conversion
        svg_file=f"./data/{sample_name}/qupath_exports/nuclei.svg",
        image_height=4096,
        image_width=4096,
        
        # Spatial segmentation
        visium_data_path=f"./data/{sample_name}/visium/",
        source_image_path=f"./data/{sample_name}/tissue_hires_image.tif",
        
        # Cell annotation
        classifier_path=f"./models/clf_{sample_name}.joblib",
        hd_expression_path=f"./data/{sample_name}/hd_data/matrix.h5",
        hd_positions_path=f"./data/{sample_name}/hd_data/positions.parquet",
        
        # Visualization
        background_image_path=f"./data/{sample_name}/tissue_fullres_image.tif",
        color_scheme="primary",
        rename_cell_types=True,
        
        # Computational resources
        num_processes=80,
        memory_limit_gb=600,
        
        # General settings
        verbose=True
    )
    
    # Save with documentation comments
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        f.write(f"""# Spatialcell Pipeline Configuration
# Generated example for sample: {sample_name}

# Basic sample information
sample_name: {config.sample_name}
output_dir: {config.output_dir}

# Pipeline step controls (set to false to skip steps)
enable_roi_extraction: {config.enable_roi_extraction}
enable_svg_conversion: {config.enable_svg_conversion}
enable_spatial_segmentation: {config.enable_spatial_segmentation}
enable_cell_annotation: {config.enable_cell_annotation}
enable_visualization: {config.enable_visualization}

# ROI extraction parameters
loupe_export_dir: {config.loupe_export_dir}

# SVG conversion parameters
svg_file: {config.svg_file}
image_height: {config.image_height}
image_width: {config.image_width}

# Spatial segmentation parameters
visium_data_path: {config.visium_data_path}
source_image_path: {config.source_image_path}

# StarDist parameters
stardist_prob_thresh: {config.stardist_prob_thresh}  # 0.01-0.3 range
stardist_nms_thresh: {config.stardist_nms_thresh}    # 0.3-0.7 range

# Label expansion parameters
labels_key: {config.labels_key}                      # labels_qupath_expanded, labels_joint, etc.
expansion_algorithm: {config.expansion_algorithm}    # max_bin_distance or volume_ratio
max_bin_distance: {config.max_bin_distance}
volume_ratio: {config.volume_ratio}
nearest_neighbors: {config.nearest_neighbors}
subset_pca: {config.subset_pca}

# Cell annotation parameters
classifier_path: {config.classifier_path}
hd_expression_path: {config.hd_expression_path}
hd_positions_path: {config.hd_positions_path}
cell_labels_column: {config.cell_labels_column}

# Multi-scale classification parameters (in micrometers)
min_scale: {config.min_scale}
max_scale: {config.max_scale}

# Visualization parameters
background_image_path: {config.background_image_path}
color_scheme: {config.color_scheme}                   # primary, scientific, functional, modern, warm, golden
point_size: {config.point_size}
point_shape: {config.point_shape}
rename_cell_types: {config.rename_cell_types}

# Computational resources
num_processes: {config.num_processes}
num_threads: {config.num_threads}
memory_limit_gb: {config.memory_limit_gb}

# General settings
verbose: {config.verbose}
""")
    
    logging.info(f"Example configuration saved to: {output_path}")


if __name__ == "__main__":
    # Example usage and testing
    logging.basicConfig(level=logging.INFO)
    
    print("Pipeline Configuration Management")
    print("=" * 40)
    
    # Create example configuration
    example_config = PipelineConfig(
        sample_name="test_sample",
        output_dir="./test_output",
        verbose=True
    )
    
    print("Example configuration:")
    print(example_config.to_yaml())
    
    # Test validation
    try:
        validate_pipeline_inputs(example_config)
        print("Configuration validation: PASSED")
    except ValueError as e:
        print(f"Configuration validation: FAILED - {e}")# Version: 1.0.7
