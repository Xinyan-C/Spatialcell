#!/usr/bin/env python3
"""
Comprehensive Color Schemes for Spatialcell Pipeline

This module provides scientifically-designed color palettes optimized for spatial 
transcriptomics visualization, with special consideration for cell type relationships,
biological function groupings, and publication quality.

Author: Xinyan
License: MIT
"""

import random
import itertools
import logging
from typing import List, Dict, Optional, Set, Any
import warnings


# =============================================================================
# Primary Cell Type Color Mapping (Current Default)
# =============================================================================

PRIMARY_CELL_TYPE_COLORS = {
    # Bone lineage - Warm colors representing osteogenic progression
    "Axin2⁺-SuSC": "#E31A1C",           # Bright red - Axin2+ suture stem cells
    "Gli1⁺-SuSC": "#33A02C",            # Forest green - Gli1+ suture stem cells  
    "Pre-osteoblast": "#fcc43c",         # Golden yellow - Pre-osteoblasts
    "Proliferating preosteoblast": "#f4bccc",  # Light pink - Proliferating pre-osteoblasts
    "Mmp13⁺-Pre-osteoblast": "#246c64", # Teal - Mmp13+ pre-osteoblasts
    "Immature osteoblast": "#377EB8",    # Royal blue - Immature osteoblasts
    "Mature osteoblast": "#7AC5C1",      # Aqua - Mature osteoblasts

    # Cartilage and mesenchyme - Purple/orange family
    "Suture mesenchyme": "#CAB2D6",      # Light lavender - Suture mesenchyme
    "Ligament-like mesenchyme": "#ec5414", # Orange-red - Ligament-like mesenchyme
    "Ectocranial SM": "#FF7F00",         # Orange - Ectocranial smooth muscle

    # Immune lineage - Varied colors for distinction
    "M1 Macrophage": "#e4bcd4",          # Light pink - M1 macrophages
    "M2 Macrophage": "#7570B3",          # Purple - M2 macrophages
    "Osteoclast": "#E7298A",             # Magenta - Osteoclasts
    "Mast cell": "#666666",              # Gray - Mast cells
    "Naive B cell": "#4c7c4c",           # Dark green - Naive B cells
    "Neutrophil": "#377EB8",             # Blue - Neutrophils
    "Proliferating B cell": "#ec7c5c",   # Coral - Proliferating B cells

    # Meningeal fibroblasts - Blue/green family
    "Arachnoid fibroblast": "#d1a3a4",   # Rose - Arachnoid fibroblasts
    "Inner dural fibroblast": "#8f61c1", # Purple - Inner dural fibroblasts
    "Outer dural fibroblast": "#6A3D9A", # Dark purple - Outer dural fibroblasts
    "Pro-inflammatory dural fibroblast": "#b4d434", # Lime green - Pro-inflammatory dural fibroblasts
    "Proliferating meningeal fibroblast": "#f4bccc", # Light pink - Proliferating meningeal fibroblasts

    # Neural lineage - Cool colors
    "Glial cell": "#bc80bd",             # Lavender - Glial cells
    "Neural progenitor cell": "#f4b484", # Peach - Neural progenitor cells
    "Neuron": "#80B1D3",                 # Light blue - Neurons

    # Other stromal cells - Earth tones
    "Chondrocyte": "#acd4ec",            # Light blue - Chondrocytes
    "Capillary endothelial cell": "#B15928", # Brown - Capillary endothelial cells
    "Erythroid progenitor cell": "#fbefcc",  # Cream - Erythroid progenitor cells
    "Pericyte": "#A6761D"                # Dark gold - Pericytes
}


# =============================================================================
# Alternative Color Schemes for Different Use Cases
# =============================================================================

# Scientific Journal Standard - High contrast, colorblind-friendly
SCIENTIFIC_JOURNAL_SCHEME = {
    'Gli1⁺-SuSC': '#E31A1C',           # Bright red - Stem cell 1
    'Axin2⁺-SuSC': '#FF7F00',          # Orange - Stem cell 2  
    'Pre-osteoblast': '#33A02C',        # Green - Pre-osteoblasts
    'Suture mesenchyme': '#F781BF',     # Pink - Suture mesenchyme
    'Inner dural fibroblast': '#1F78B4', # Blue - Inner dural fibroblasts
    'Outer dural fibroblast': '#6A3D9A', # Purple - Outer dural fibroblasts
    'Immature osteoblast': '#A6CEE3',   # Light blue - Immature osteoblasts
    'Mature osteoblast': '#B2DF8A',     # Light green - Mature osteoblasts
    'Ectocranial SM': '#FDBF6F',        # Light orange - Ectocranial smooth muscle
    'Ligament-like mesenchyme': '#CAB2D6', # Light purple - Ligament-like mesenchyme
    'Proliferating meningeal fibroblast': '#FFFF99', # Yellow - Proliferating meningeal fibroblasts
    'Osteoclast': '#FB9A99',            # Light red - Osteoclasts
    'Chondrocyte': '#B15928',           # Brown - Chondrocytes
    'M2 Macrophage': '#666666',         # Gray - M2 macrophages
    'Capillary endothelial cell': '#377EB8', # Deep blue - Capillary endothelial cells
    'Pericyte': '#4DAF4A',              # Medium green - Pericytes
    'Mast cell': '#984EA3',             # Deep purple - Mast cells
    'Neutrophil': '#FF1493',            # Deep pink - Neutrophils
    'Neuron': '#00CED1',                # Cyan - Neurons
    'B cell': '#8B4513'                 # Dark brown - B cells
}

# Functional Grouping - Colors grouped by biological function
FUNCTIONAL_GROUPING_SCHEME = {
    # Stem cell series - Red family
    'Gli1⁺-SuSC': '#E74C3C',           # Deep red
    'Axin2⁺-SuSC': '#EC7063',          # Medium red
    'Suture mesenchyme': '#F1948A',     # Light red
    
    # Osteoblast series - Blue-green family  
    'Pre-osteoblast': '#3498DB',        # Blue
    'Immature osteoblast': '#5DADE2',   # Medium blue
    'Mature osteoblast': '#85C1E9',     # Light blue
    
    # Fibroblast series - Green family
    'Inner dural fibroblast': '#27AE60', # Deep green
    'Outer dural fibroblast': '#58D68D', # Medium green
    'Proliferating meningeal fibroblast': '#82E0AA', # Light green
    'Ligament-like mesenchyme': '#ABEBC6', # Very light green
    
    # Vascular system - Purple family
    'Capillary endothelial cell': '#8E44AD', # Deep purple
    'Pericyte': '#BB8FCE',              # Light purple
    
    # Immune cells - Orange family
    'M2 Macrophage': '#E67E22',         # Deep orange
    'Mast cell': '#F39C12',             # Medium orange
    'Neutrophil': '#F7DC6F',            # Light orange
    'B cell': '#FAD5A5',                # Very light orange
    
    # Special cell types - Unique colors
    'Osteoclast': '#E91E63',            # Rose red - Osteoclasts
    'Chondrocyte': '#9C27B0',           # Purple-red - Chondrocytes
    'Ectocranial SM': '#FF5722',        # Deep orange-red - Smooth muscle
    'Neuron': '#00BCD4'                 # Cyan-blue - Neurons
}

# Modern Gradient - Contemporary design aesthetics
MODERN_GRADIENT_SCHEME = {
    'Gli1⁺-SuSC': '#667eea',           # Blue-purple
    'Axin2⁺-SuSC': '#764ba2',          # Deep purple
    'Pre-osteoblast': '#f093fb',        # Pink-purple
    'Suture mesenchyme': '#f5576c',     # Pink-red
    'Inner dural fibroblast': '#4facfe', # Bright blue
    'Outer dural fibroblast': '#00f2fe', # Cyan
    'Immature osteoblast': '#43e97b',   # Emerald green
    'Mature osteoblast': '#38f9d7',     # Mint green
    'Ectocranial SM': '#ffecd2',        # Warm cream
    'Ligament-like mesenchyme': '#fcb69f', # Peach
    'Proliferating meningeal fibroblast': '#a8edea', # Light cyan
    'Osteoclast': '#ffd89b',            # Golden yellow
    'Chondrocyte': '#19547b',           # Deep blue
    'M2 Macrophage': '#ff6b6b',         # Coral red
    'Capillary endothelial cell': '#4ecdc4', # Mint blue
    'Pericyte': '#45b7d1',              # Sky blue
    'Mast cell': '#96ceb4',             # Mint green
    'Neutrophil': '#ffeaa7',            # Light yellow
    'Neuron': '#dda0dd',                # Plum
    'B cell': '#98d8c8'                 # Sea green
}

# Warm Tone - Comfortable for extended viewing
WARM_TONE_SCHEME = {
    'Gli1⁺-SuSC': '#FF8A80',           # Coral pink
    'Axin2⁺-SuSC': '#FFAB91',          # Peach
    'Pre-osteoblast': '#FFCC02',        # Golden yellow
    'Suture mesenchyme': '#FFF176',     # Light yellow
    'Inner dural fibroblast': '#AED581', # Light green
    'Outer dural fibroblast': '#81C784', # Medium green
    'Immature osteoblast': '#4DB6AC',   # Tiffany blue
    'Mature osteoblast': '#64B5F6',     # Sky blue
    'Ectocranial SM': '#9575CD',        # Lavender
    'Ligament-like mesenchyme': '#F06292', # Pink
    'Proliferating meningeal fibroblast': '#A1887F', # Warm gray
    'Osteoclast': '#FF7043',            # Orange-red
    'Chondrocyte': '#8D6E63',           # Brown-gray
    'M2 Macrophage': '#FFB74D',         # Warm orange
    'Capillary endothelial cell': '#90A4AE', # Blue-gray
    'Pericyte': '#BCAAA4',              # Warm gray
    'Mast cell': '#CE93D8',             # Light purple
    'Neutrophil': '#80CBC4',            # Mint
    'Neuron': '#E1BEE7',                # Light purple
    'B cell': '#C5E1A5'                 # Light green
}

# Golden Ratio - Mathematical beauty based on Fibonacci sequence
GOLDEN_RATIO_SCHEME = {
    'Gli1⁺-SuSC': '#C41E3A',           # Deep red
    'Axin2⁺-SuSC': '#E6B800',          # Golden yellow
    'Pre-osteoblast': '#228B22',        # Forest green
    'Suture mesenchyme': '#FF6347',     # Tomato red
    'Inner dural fibroblast': '#4169E1', # Royal blue
    'Outer dural fibroblast': '#9932CC', # Dark orchid
    'Immature osteoblast': '#FF1493',   # Deep pink
    'Mature osteoblast': '#00CED1',     # Dark turquoise
    'Ectocranial SM': '#FF8C00',        # Dark orange
    'Ligament-like mesenchyme': '#8B008B', # Dark magenta
    'Proliferating meningeal fibroblast': '#DC143C', # Crimson
    'Osteoclast': '#32CD32',            # Lime green
    'Chondrocyte': '#FF69B4',           # Hot pink
    'M2 Macrophage': '#1E90FF',         # Dodger blue
    'Capillary endothelial cell': '#FFD700', # Gold
    'Pericyte': '#9370DB',              # Medium purple
    'Mast cell': '#FF4500',             # Orange red
    'Neutrophil': '#00FA9A',            # Medium spring green
    'Neuron': '#8A2BE2',                # Blue violet
    'B cell': '#FF6B6B'                 # Light red
}


# =============================================================================
# Comprehensive Color Palette Collections
# =============================================================================

# Extended color palette organized by color families
EXTENDED_COLOR_PALETTE = {
    'green_family': [
        "#9cbc1c", "#c4dc4c", "#b4d434", "#84a42c", "#74944c", "#6c8c54", 
        "#6c8c6c", "#6c945c", "#6c8454", "#5c8454", "#4c8444", "#4c7c4c", 
        "#4c6c4c", "#446c34", "#2c6c3c", "#5c745c", "#4c644c", "#3c543c"
    ],
    'green_sage': [
        "#acbc7c", "#c4d49c", "#c4d494", "#acbc8c", "#acb48c", "#acb494", 
        "#94a47c", "#bcc494", "#d4d4ac", "#bcbc9c", "#b4b49c", "#acac84", 
        "#aca474", "#ac9c74", "#9c9c84", "#949474", "#94945c", "#948c4c", 
        "#848c4c", "#847c54", "#7c744c", "#5c5434"
    ],
    'green_mint': [
        "#9cbcac", "#d4ece4", "#d4e4d4", "#b4d4cc", "#ccd4c4", "#bcd4bc", 
        "#bcccb4", "#b4bcac", "#a4bcb4", "#9cb494", "#84a484", "#84a494", 
        "#849c8c", "#748c74", "#7c8c6c", "#7c846c", "#84948c", "#6c7c74", 
        "#5c746c", "#545c4c", "#444c3c"
    ],
    'red_vibrant': [
        "#e40414", "#ec5414", "#ec4c2c", "#d42c24", "#d43c1c", "#c42c1c", 
        "#bc1c34", "#cc141c", "#a42424", "#ac1c24", "#8c1c24", "#842424", 
        "#7c1c1c", "#641414", "#4c1c1c"
    ],
    'red_coral': [
        "#ec6c3c", "#f4c4a4", "#fcc4b4", "#f4a494", "#f4b484", "#f49c74", 
        "#f48c64", "#ec845c", "#ec7c5c", "#e47464", "#d47c64", "#dc6c4c", 
        "#ec6c44", "#cc543c", "#d44434", "#bc5c4c", "#bc5444", "#b45c44", 
        "#bc3c24", "#b43c2c", "#a44434", "#944424", "#8c3c2c", "#9c2c24", 
        "#94342c", "#843424", "#642c1c"
    ],
    'pink_family': [
        "#f494a4", "#fcd4e4", "#e4ccd4", "#ecc4cc", "#f4bccc", "#ecb4c4", 
        "#e4acbc", "#e4a4ac", "#dcacac", "#dca4b4", "#cc949c", "#cc8c94", 
        "#dc6c84", "#dc6c7c", "#cc7474", "#bc7c8c", "#b47474", "#a45c64", 
        "#c45c6c", "#c45464", "#c45c5c", "#b44454", "#9c6c6c", "#ac5c5c", 
        "#9c4c54", "#74545c", "#84444c", "#7c444c", "#6c3434", "#642c2c"
    ],
    'purple_family': [
        "#bca4cc", "#f4cce4", "#dcc4e4", "#e4bcd4", "#cc94bc", "#bc7cb4", 
        "#a47cb4", "#a4549c", "#7c5484", "#7c547c", "#8c1c74", "#643c74", 
        "#6c246c", "#442454"
    ],
    'yellow_gold': [
        "#fcc43c", "#fcf49c", "#fcec6c", "#ecd454", "#fcdc5c", "#f4cc64", 
        "#e3ab4b", "#dc9434", "#dc8c3c", "#e48c3c", "#d4a434", "#dc9c34", 
        "#d48434", "#bc8454", "#bc8c2c", "#cc5c24", "#c47c14", "#bc6c34", 
        "#9c6424", "#9c5424"
    ],
    'blue_family': [
        "#2c5ca4", "#acd4ec", "#a4bcdc", "#8caccc", "#8cacdc", "#6c9cc4", 
        "#6c94cc", "#5c74bc", "#4c94c4", "#3474ac", "#146c9c", "#14547c", 
        "#04446c", "#043c74", "#043464", "#1c345c", "#14244c", "#1c2c44", 
        "#645cac"
    ],
    'gray_neutral': [
        "#ececec", "#f4f4f4", "#ececf4", "#f4f4ec", "#ececec", "#ece4e4", 
        "#dce4dc", "#e4e4d4", "#dcd4c4", "#d4d4cc", "#d3cbc3", "#c4bcb4", 
        "#ececdc", "#f4fce4", "#ececdc", "#e4dcc4", "#c4c4b4", "#ccccc4", 
        "#bcc4bc", "#b4b4b4", "#bcc4b4", "#bcc4ac", "#a4aca4", "#aca49c", 
        "#7c7c7c", "#6c6c64", "#5c544c", "#44443c", "#34342c"
    ]
}

# All available color schemes
AVAILABLE_SCHEMES = {
    'primary': PRIMARY_CELL_TYPE_COLORS,
    'scientific': SCIENTIFIC_JOURNAL_SCHEME,
    'functional': FUNCTIONAL_GROUPING_SCHEME,
    'modern': MODERN_GRADIENT_SCHEME,
    'warm': WARM_TONE_SCHEME,
    'golden': GOLDEN_RATIO_SCHEME
}


# =============================================================================
# Color Assignment Functions
# =============================================================================

def get_cell_type_colors(cell_types: List[str], 
                        scheme: str = 'primary',
                        random_seed: Optional[int] = None) -> Dict[str, str]:
    """
    Assign colors to cell types using specified color scheme.
    
    For cell types not in the predefined scheme, colors are assigned from
    the extended palette to avoid conflicts with existing assignments.
    
    Args:
        cell_types (List[str]): List of unique cell type names
        scheme (str): Color scheme name ('primary', 'scientific', 'functional', etc.)
        random_seed (int, optional): Random seed for reproducible color assignment
        
    Returns:
        Dict[str, str]: Mapping of cell types to hex color codes
        
    Raises:
        ValueError: If insufficient unique colors available
        KeyError: If color scheme not found
    """
    if random_seed is not None:
        random.seed(random_seed)
    
    if scheme not in AVAILABLE_SCHEMES:
        available = list(AVAILABLE_SCHEMES.keys())
        raise KeyError(f"Color scheme '{scheme}' not found. Available schemes: {available}")
    
    fixed_colors = AVAILABLE_SCHEMES[scheme]
    cell_type_colors = {}
    
    # Get colors already used in the fixed scheme
    used_colors = set(fixed_colors.values())
    
    # Create pool of available colors from extended palette
    all_palette_colors = list(itertools.chain(*EXTENDED_COLOR_PALETTE.values()))
    available_colors = [color for color in all_palette_colors if color not in used_colors]
    
    # Assign colors to cell types
    for cell_type in cell_types:
        if cell_type in fixed_colors:
            # Use predefined color from scheme
            cell_type_colors[cell_type] = fixed_colors[cell_type]
        else:
            # Assign from available colors
            if available_colors:
                chosen_color = random.choice(available_colors)
                cell_type_colors[cell_type] = chosen_color
                available_colors.remove(chosen_color)
                
                logging.info(f"Assigned color {chosen_color} to cell type '{cell_type}' "
                           f"(not in {scheme} scheme)")
            else:
                # Fallback to matplotlib default colors if extended palette exhausted
                warnings.warn(f"Extended color palette exhausted for cell type '{cell_type}'. "
                            "Consider using fewer cell types or defining custom colors.")
                # Generate a random hex color as last resort
                fallback_color = f"#{random.randint(0, 0xFFFFFF):06x}"
                cell_type_colors[cell_type] = fallback_color
                
                logging.warning(f"Used fallback color {fallback_color} for '{cell_type}'")
    
    return cell_type_colors


def get_color_families() -> Dict[str, List[str]]:
    """
    Get organized color families for custom color selection.
    
    Returns:
        Dict[str, List[str]]: Color families organized by theme
    """
    return EXTENDED_COLOR_PALETTE.copy()


def list_available_schemes() -> List[str]:
    """
    List all available color schemes.
    
    Returns:
        List[str]: Names of available color schemes
    """
    return list(AVAILABLE_SCHEMES.keys())


def get_scheme_description(scheme: str) -> str:
    """
    Get description of a color scheme.
    
    Args:
        scheme (str): Color scheme name
        
    Returns:
        str: Description of the color scheme
        
    Raises:
        KeyError: If scheme not found
    """
    descriptions = {
        'primary': "Current default scheme optimized for bone/cartilage lineages",
        'scientific': "High contrast, colorblind-friendly scheme for publications",
        'functional': "Colors grouped by biological function and cell relationships",
        'modern': "Contemporary gradient-based aesthetic for presentations",
        'warm': "Warm tone scheme comfortable for extended viewing",
        'golden': "Mathematical beauty based on golden ratio (1.618) color spacing"
    }
    
    if scheme not in descriptions:
        raise KeyError(f"Scheme '{scheme}' not found")
    
    return descriptions[scheme]


def validate_color_scheme(colors: Dict[str, str]) -> bool:
    """
    Validate that a color scheme has sufficient contrast and uniqueness.
    
    Args:
        colors (Dict[str, str]): Color mapping to validate
        
    Returns:
        bool: True if color scheme passes validation
    """
    # Check for duplicate colors
    color_values = list(colors.values())
    if len(color_values) != len(set(color_values)):
        logging.warning("Color scheme contains duplicate colors")
        return False
    
    # Check for valid hex format
    for cell_type, color in colors.items():
        if not (color.startswith('#') and len(color) == 7):
            try:
                int(color[1:], 16)  # Try to parse as hex
            except ValueError:
                logging.error(f"Invalid color format for '{cell_type}': {color}")
                return False
    
    logging.info(f"Color scheme validation passed for {len(colors)} cell types")
    return True

def get_color_scheme(scheme_name="default"):
    """Get color scheme for visualizations (simplified interface)"""
    
    schemes = {
        "default": [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
        ],
        "primary": list(PRIMARY_CELL_TYPE_COLORS.values()),
        "scientific": list(SCIENTIFIC_JOURNAL_SCHEME.values()),
        "functional": list(FUNCTIONAL_GROUPING_SCHEME.values()),
        "modern": list(MODERN_GRADIENT_SCHEME.values()),
        "warm": list(WARM_TONE_SCHEME.values()),
        "golden": list(GOLDEN_RATIO_SCHEME.values())
    }
    
    return schemes.get(scheme_name, schemes["default"])

# =============================================================================
# Backward Compatibility Aliases
# =============================================================================

# Maintain compatibility with original function name
get_type_to_color = get_cell_type_colors
fixed_type_to_color = PRIMARY_CELL_TYPE_COLORS


# =============================================================================
# Module Information
# =============================================================================

__version__ = "1.0.0"
__author__ = "Xinyan"
__all__ = [
    'PRIMARY_CELL_TYPE_COLORS',
    'AVAILABLE_SCHEMES',
    'EXTENDED_COLOR_PALETTE',
    'get_cell_type_colors',
	'get_color_scheme', 
    'get_color_families', 
    'list_available_schemes',
    'get_scheme_description',
    'validate_color_scheme',
    # Compatibility aliases
    'get_type_to_color',
    'fixed_type_to_color'
]

if __name__ == "__main__":
    # Example usage and testing
    logging.basicConfig(level=logging.INFO)
    
    print("Spatialcell Color Schemes Module")
    print(f"Version: {__version__}")
    print(f"Available schemes: {list_available_schemes()}")
    
    # Test color assignment
    test_cell_types = ['Gli1⁺-SuSC', 'Pre-osteoblast', 'Unknown_Cell_Type']
    colors = get_cell_type_colors(test_cell_types, scheme='primary')
    
    print(f"\nExample color assignment:")
    for cell_type, color in colors.items():
        print(f"  {cell_type}: {color}")
    
    print(f"\nValidation result: {validate_color_scheme(colors)}")# Version: 1.0.7
