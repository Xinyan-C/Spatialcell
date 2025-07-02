"""
SpatialCell: Integrated Spatial Transcriptomics Analysis Pipeline

Author: Xinyan
Email: keepandon@gmail.com
GitHub: https://github.com/Xinyan-C/Spatialcell
"""

__version__ = "1.0.5"
__author__ = "Xinyan"
__email__ = "keepandon@gmail.com"

# Core imports (always available)
__all__ = ['__version__', '__author__', '__email__', 'get_version', 'citation', 'info']

# Check for optional dependencies
_topact_available = False
_bin2cell_available = False

try:
    import topact
    _topact_available = True
except ImportError:
    pass

try:
    import bin2cell
    _bin2cell_available = True
except ImportError:
    pass

# Import available modules
try:
    from .preprocessing.svg_to_npz import convert_svg_to_npz
    from .utils.roi_extractor import ROIExtractor
    __all__.extend(['convert_svg_to_npz', 'ROIExtractor'])
except ImportError:
    pass

if _bin2cell_available:
    try:
        from .spatial_segmentation.spatial_processor import SpatialSegmentationProcessor
        __all__.append('SpatialSegmentationProcessor')
    except ImportError:
        pass

if _topact_available:
    try:
        from .cell_annotation.annotation_processor import CellAnnotationProcessor
        from .utils.config_manager import load_config
        __all__.extend(['CellAnnotationProcessor', 'load_config'])
    except ImportError:
        pass

if _topact_available and _bin2cell_available:
    try:
        from .workflows.complete_pipeline import SpatialCellPipeline
        __all__.append('SpatialCellPipeline')
    except ImportError:
        pass


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
    """Print package information and installation status."""
    print(f"SpatialCell v{__version__}")
    print(f"Author: {__author__}")
    print(f"Email: {__email__}")
    print("GitHub: https://github.com/Xinyan-C/Spatialcell")
    print()
    
    # Show feature status
    print("Feature Status:")
    print(f"  Basic preprocessing: Available")
    print(f"  Spatial segmentation (bin2cell): {'Available' if _bin2cell_available else 'Missing'}")
    print(f"  Cell annotation (topact): {'Available' if _topact_available else 'Missing'}")
    
    # Show installation guidance for missing features
    if not _topact_available:
        print()
        print("For complete functionality, install TopAct:")
        print("  pip install git+https://gitlab.com/kfbenjamin/topact.git")
        print("  Official: https://gitlab.com/kfbenjamin/topact")
        print("  Note: TopACT is licensed under GPL v3")


# Show guidance message on import if topact is missing
if not _topact_available:
    print("SpatialCell: TopAct not found. For complete cell annotation functionality, run:")
    print("pip install git+https://gitlab.com/kfbenjamin/topact.git")
    print("Official repository: https://gitlab.com/kfbenjamin/topact")# Version: 1.0.7
