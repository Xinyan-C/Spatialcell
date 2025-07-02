#!/usr/bin/env python3
"""
Flexible cell type renaming utilities for SpatialCell

This module provides utilities for converting between abbreviated and full names
of cell types, allowing users to define custom mapping rules for better 
visualization and analysis.

Author: Xinyan
License: MIT
"""

import json
import yaml
from typing import Dict, List, Optional, Union
import warnings


class CellTypeRenamer:
    """
    Flexible cell type renaming utility.
    
    Supports custom mapping rules for converting between different naming conventions,
    such as abbreviated codes to full names or vice versa.
    """
    
    def __init__(self, rename_dict: Optional[Dict[str, str]] = None):
        """
        Initialize the renamer with custom mapping rules.
        
        Args:
            rename_dict (Dict[str, str], optional): Mapping from original to new names
        """
        self.rename_dict = rename_dict or {}
    
    def load_from_dict(self, rename_dict: Dict[str, str]) -> None:
        """Load renaming rules from a dictionary."""
        self.rename_dict.update(rename_dict)
    
    def load_from_file(self, file_path: str) -> None:
        """
        Load renaming rules from a file (JSON or YAML).
        
        Args:
            file_path (str): Path to configuration file
        """
        try:
            if file_path.endswith('.json'):
                with open(file_path, 'r') as f:
                    rules = json.load(f)
            elif file_path.endswith(('.yml', '.yaml')):
                with open(file_path, 'r') as f:
                    rules = yaml.safe_load(f)
            else:
                raise ValueError("File must be JSON or YAML format")
            
            self.rename_dict.update(rules)
            
        except Exception as e:
            warnings.warn(f"Failed to load renaming rules from {file_path}: {e}")
    
    def rename_single(self, name: str) -> str:
        """Rename a single cell type."""
        return self.rename_dict.get(name, name)
    
    def rename_list(self, names: List[str]) -> List[str]:
        """Rename a list of cell types."""
        return [self.rename_single(name) for name in names]
    
    def reverse_mapping(self) -> Dict[str, str]:
        """Get the reverse mapping (new_name -> original_name)."""
        return {v: k for k, v in self.rename_dict.items()}
    
    def get_available_mappings(self) -> Dict[str, str]:
        """Get all available mappings."""
        return self.rename_dict.copy()


# Common renaming utilities
def create_abbreviation_mapper(abbreviations: Dict[str, str]) -> CellTypeRenamer:
    """
    Create a renamer for converting abbreviations to full names.
    
    Args:
        abbreviations (Dict[str, str]): Mapping from abbreviation to full name
        
    Returns:
        CellTypeRenamer: Configured renamer instance
        
    Example:
        abbrev_map = {
            "SuSC": "Suture Stem Cell",
            "Pre-OB": "Pre-osteoblast",
            "OB": "Osteoblast"
        }
        renamer = create_abbreviation_mapper(abbrev_map)
    """
    return CellTypeRenamer(abbreviations)


def create_standardization_mapper(standards: Dict[str, str]) -> CellTypeRenamer:
    """
    Create a renamer for standardizing cell type names.
    
    Args:
        standards (Dict[str, str]): Mapping from variant to standard name
        
    Returns:
        CellTypeRenamer: Configured renamer instance
        
    Example:
        standard_map = {
            "B cells": "B cell",
            "B-cells": "B cell", 
            "B_cell": "B cell"
        }
        renamer = create_standardization_mapper(standard_map)
    """
    return CellTypeRenamer(standards)


# Pre-defined common mappings (examples)
COMMON_ABBREVIATIONS = {
    # Bone/cartilage lineage
    "SuSC": "Suture Stem Cell",
    "Pre-OB": "Pre-osteoblast", 
    "OB": "Osteoblast",
    "OC": "Osteoclast",
    "CC": "Chondrocyte",
    
    # Immune cells
    "MP": "Macrophage",
    "DC": "Dendritic cell",
    "BC": "B cell",
    "TC": "T cell",
    "NK": "Natural killer cell",
    
    # Stromal cells
    "FB": "Fibroblast",
    "EC": "Endothelial cell",
    "PC": "Pericyte",
    "SM": "Smooth muscle",
    
    # Neural
    "NPC": "Neural progenitor cell",
    "NEU": "Neuron",
    "AST": "Astrocyte",
    "OLI": "Oligodendrocyte"
}

COMMON_STANDARDIZATIONS = {
    # Pluralization standardization
    "B cells": "B cell",
    "T cells": "T cell", 
    "Neurons": "Neuron",
    "Osteoblasts": "Osteoblast",
    "Osteoclasts": "Osteoclast",
    "Chondrocytes": "Chondrocyte",
    "Macrophages": "Macrophage",
    "Fibroblasts": "Fibroblast",
    
    # Capitalization standardization  
    "b cell": "B cell",
    "t cell": "T cell",
    "macrophage": "Macrophage",
    "osteoblast": "Osteoblast",
    
    # Separator standardization
    "B_cell": "B cell",
    "T_cell": "T cell",
    "Pre_osteoblast": "Pre-osteoblast"
}


# Global renamer instances (for backward compatibility)
_default_renamer = None


def get_rename_dict(sample_name: str = None) -> Dict[str, str]:
    """
    Get renaming dictionary (backward compatibility function).
    
    Args:
        sample_name (str, optional): Sample name (ignored in new version)
        
    Returns:
        Dict[str, str]: Common standardization mappings
        
    Note:
        This function maintains backward compatibility with the old interface.
        For new code, use CellTypeRenamer class directly.
    """
    return COMMON_STANDARDIZATIONS.copy()


def rename_cell_types(data_frame, sample_name: str = None, 
                     cell_type_column: str = 'cell_type',
                     custom_mapping: Optional[Dict[str, str]] = None):
    """
    Rename cell types in a DataFrame.
    
    Args:
        data_frame: DataFrame containing cell type information
        sample_name (str, optional): Sample name (for backward compatibility)
        cell_type_column (str): Column name containing cell types
        custom_mapping (Dict[str, str], optional): Custom renaming rules
        
    Returns:
        DataFrame with renamed cell types
    """
    import pandas as pd
    
    mapping = custom_mapping or COMMON_STANDARDIZATIONS
    renamer = CellTypeRenamer(mapping)
    
    df_copy = data_frame.copy()
    df_copy[cell_type_column] = df_copy[cell_type_column].apply(renamer.rename_single)
    
    return df_copy


# Convenience functions
def expand_abbreviations(cell_types: Union[str, List[str]], 
                        custom_mapping: Optional[Dict[str, str]] = None) -> Union[str, List[str]]:
    """
    Expand abbreviated cell type names to full names.
    
    Args:
        cell_types: Single cell type or list of cell types
        custom_mapping: Custom abbreviation mapping
        
    Returns:
        Expanded cell type name(s)
    """
    mapping = custom_mapping or COMMON_ABBREVIATIONS
    renamer = CellTypeRenamer(mapping)
    
    if isinstance(cell_types, str):
        return renamer.rename_single(cell_types)
    else:
        return renamer.rename_list(cell_types)


def standardize_names(cell_types: Union[str, List[str]],
                     custom_mapping: Optional[Dict[str, str]] = None) -> Union[str, List[str]]:
    """
    Standardize cell type names to consistent format.
    
    Args:
        cell_types: Single cell type or list of cell types
        custom_mapping: Custom standardization mapping
        
    Returns:
        Standardized cell type name(s)
    """
    mapping = custom_mapping or COMMON_STANDARDIZATIONS  
    renamer = CellTypeRenamer(mapping)
    
    if isinstance(cell_types, str):
        return renamer.rename_single(cell_types)
    else:
        return renamer.rename_list(cell_types)


def create_config_template(output_path: str) -> None:
    """
    Create a template configuration file for custom renaming rules.
    
    Args:
        output_path (str): Path where to save the template
    """
    template = {
        "abbreviations": {
            "SuSC": "Suture Stem Cell",
            "Pre-OB": "Pre-osteoblast",
            "OB": "Osteoblast"
        },
        "standardizations": {
            "B cells": "B cell",
            "Osteoblasts": "Osteoblast"
        },
        "custom_mappings": {
            "Your_Custom_Name": "Standard Name"
        }
    }
    
    if output_path.endswith('.json'):
        with open(output_path, 'w') as f:
            json.dump(template, f, indent=2)
    else:
        with open(output_path, 'w') as f:
            yaml.dump(template, f, default_flow_style=False)
    
    print(f"Template configuration saved to: {output_path}")


# Module information
__version__ = "1.0.0"
__author__ = "Xinyan"
__all__ = [
    'CellTypeRenamer',
    'create_abbreviation_mapper',
    'create_standardization_mapper', 
    'get_rename_dict',
    'rename_cell_types',
    'expand_abbreviations',
    'standardize_names',
    'create_config_template',
    'COMMON_ABBREVIATIONS',
    'COMMON_STANDARDIZATIONS'
]


if __name__ == "__main__":
    # Example usage
    print("SpatialCell Cell Type Renaming Utilities")
    print(f"Version: {__version__}")
    
    # Example 1: Basic usage
    renamer = CellTypeRenamer({"OB": "Osteoblast", "OC": "Osteoclast"})
    print(f"\nExample 1: {renamer.rename_single('OB')}")
    
    # Example 2: Abbreviation expansion
    expanded = expand_abbreviations(["SuSC", "Pre-OB", "OB"])
    print(f"Example 2: {expanded}")
    
    # Example 3: Name standardization
    standardized = standardize_names(["B cells", "Osteoblasts"])
    print(f"Example 3: {standardized}")
# Version: 1.0.7
