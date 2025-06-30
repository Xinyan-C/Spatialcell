"""
SpatialCell: Integrated Spatial Transcriptomics Analysis Pipeline

A comprehensive pipeline for spatial transcriptomics analysis that combines
cell segmentation and automated cell type annotation using QuPath, Bin2cell, and TopAct.

Author: Xinyan
Email: keepandon@gmail.com
GitHub: https://github.com/Xinyan-C/Spatialcell
"""

__version__ = "1.0.0"
__author__ = "Xinyan"
__email__ = "keepandon@gmail.com"

# Import main classes and functions for easy access
try:
    from .workflows.complete_pipeline import SpatialCellPipeline
    from .utils.config_manager import load_config
    from .preprocessing.svg_to_npz import SVGToNPZConverter
    from .spatial_segmentation.spatial_processor import SpatialSegmentationProcessor
    from .cell_annotation.annotation_processor import CellAnnotationProcessor
    
    __all__ = [
        'SpatialCellPipeline',
        'load_config',
        'SVGToNPZConverter',
        'SpatialSegmentationProcessor', 
        'CellAnnotationProcessor',
        '__version__',
        '__author__',
        '__email__'
    ]
    
except ImportError as e:
    # Handle missing dependencies gracefully
    print(f"Warning: Some SpatialCell components could not be imported: {e}")
    print("Please ensure all dependencies are installed: pip install -r requirements.txt")
    
    __all__ = ['__version__', '__author__', '__email__']


def get_version():
    """Return the current version of SpatialCell."""
    return __version__


def citation():
    """Return citation information for SpatialCell."""
    return f"""
SpatialCell: Integrated Spatial Transcriptomics Analysis Pipeline
Author: {__author__}
Version: {__version__}
GitHub: https://github.com/Xinyan-C/Spatialcell

If you use SpatialCell in your research, please cite:
@software{{spatialcell{__version__.replace('.', '')},
  author = {{{__author__}}},
  title = {{SpatialCell: Integrated Spatial Transcriptomics Analysis Pipeline}},
  url = {{https://github.com/Xinyan-C/Spatialcell}},
  version = {{{__version__}}},
  year = {{2024}}
}}
"""


def info():
    """Print package information."""
    print(f"SpatialCell v{__version__}")
    print(f"Author: {__author__}")
    print(f"Email: {__email__}")
    print("GitHub: https://github.com/Xinyan-C/Spatialcell")
    print("\nFor citation information, run: spatialcell.citation()")