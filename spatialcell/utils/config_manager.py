"""
Configuration management for SpatialCell
"""

import yaml
import os

def load_config(config_path):
    """Load configuration from YAML file"""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config

def validate_config(config):
    """Basic configuration validation"""
    required_keys = ['sample_info', 'input_paths', 'output_dir']
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")
    
    return True
# Version: 1.0.7
