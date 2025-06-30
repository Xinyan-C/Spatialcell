#!/usr/bin/env python3
"""
SpatialCell Basic Tutorial
==========================

This tutorial demonstrates how to run a complete spatial transcriptomics
analysis using the SpatialCell pipeline.

Author: Xinyan
Email: keepandon@gmail.com
"""

import os
import yaml
from pathlib import Path
import sys

# Add spatialcell to path if running from source
sys.path.append(str(Path(__file__).parent.parent))

from spatialcell.workflows.main import SpatialCellPipeline
from spatialcell.utils.config_manager import load_config


def load_configuration(config_path):
    """Load analysis configuration from YAML file"""
    print(f"Loading configuration from: {config_path}")
    config = load_config(config_path)
    return config


def run_basic_analysis():
    """Run a basic SpatialCell analysis"""
    
    # 1. Load configuration
    config_path = "examples/config_example.yml"
    config = load_configuration(config_path)
    
    # 2. Initialize pipeline
    print("Initializing SpatialCell pipeline...")
    pipeline = SpatialCellPipeline(config)
    
    # 3. Run preprocessing
    print("Step 1: Data preprocessing...")
    pipeline.run_preprocessing()
    
    # 4. Run segmentation
    print("Step 2: Spatial segmentation...")
    pipeline.run_segmentation()
    
    # 5. Run classification
    print("Step 3: Cell classification...")
    pipeline.run_classification()
    
    # 6. Generate visualizations
    print("Step 4: Creating visualizations...")
    pipeline.create_visualizations()
    
    print("✅ Analysis completed successfully!")
    print(f"Results saved to: {config['output_dir']}")


def run_step_by_step_analysis():
    """Run analysis with manual control over each step"""
    
    config_path = "examples/config_example.yml"
    config = load_configuration(config_path)
    
    # Initialize pipeline
    pipeline = SpatialCellPipeline(config)
    
    # Step 1: QuPath processing (if needed)
    print("Step 1: QuPath nucleus detection...")
    # Assuming QuPath script has been run separately
    # pipeline.run_qupath_detection()
    
    # Step 2: SVG to NPZ conversion
    print("Step 2: Converting SVG to NPZ...")
    pipeline.convert_svg_to_npz()
    
    # Step 3: Bin2cell processing
    print("Step 3: Running Bin2cell analysis...")
    pipeline.run_bin2cell()
    
    # Step 4: Classifier training (if needed)
    print("Step 4: Training/loading classifier...")
    pipeline.setup_classifier()
    
    # Step 5: TopAct classification
    print("Step 5: Running TopAct classification...")
    pipeline.run_topact()
    
    # Step 6: Visualization
    print("Step 6: Generating plots...")
    pipeline.create_plots()
    
    print("✅ Step-by-step analysis completed!")


def customize_analysis_example():
    """Example of customizing analysis parameters"""
    
    # Load base configuration
    config = load_config("examples/config_example.yml")
    
    # Customize parameters
    config['segmentation']['algorithm'] = 'volume_ratio'
    config['segmentation']['volume_ratio'] = 6.0
    config['classification']['min_scale'] = 2.0
    config['classification']['max_scale'] = 12.0
    config['visualization']['color_scheme'] = 'tab20'
    
    # Run with custom parameters
    pipeline = SpatialCellPipeline(config)
    pipeline.run_full_analysis()
    
    print("✅ Custom analysis completed!")


if __name__ == "__main__":
    print("SpatialCell Basic Tutorial")
    print("=" * 50)
    
    # Choose analysis mode
    analysis_mode = input("""
Choose analysis mode:
1. Basic analysis (recommended for beginners)
2. Step-by-step analysis (for advanced users)
3. Custom parameters example
Enter choice (1-3): """).strip()
    
    try:
        if analysis_mode == "1":
            run_basic_analysis()
        elif analysis_mode == "2":
            run_step_by_step_analysis()
        elif analysis_mode == "3":
            customize_analysis_example()
        else:
            print("Invalid choice. Running basic analysis...")
            run_basic_analysis()
            
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        print("Please check your configuration and input files.")
        sys.exit(1)