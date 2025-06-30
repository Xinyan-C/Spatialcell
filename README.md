# SpatialCell: Integrated Spatial Transcriptomics Analysis Pipeline

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub issues](https://img.shields.io/github/issues/Xinyan-C/Spatialcell)](https://github.com/Xinyan-C/Spatialcell/issues)

**SpatialCell** is an integrated computational pipeline for spatial transcriptomics analysis that combines cell segmentation and automated cell type annotation. It seamlessly integrates **QuPath** for histological image analysis, **Bin2cell** for spatial cell segmentation, and **TopAct** for machine learning-based cell classification.

## 🚀 Key Features

- **Multi-scale Cell Segmentation**: QuPath + Bin2cell integration for accurate cell boundary detection
- **Automated Cell Annotation**: TopAct-based machine learning classification
- **ROI-aware Processing**: Region-of-interest focused analysis for large datasets
- **Scalable Pipeline**: Support for multiple developmental time points (E14.5, E18.5, P3)
- **Visualization Tools**: Comprehensive plotting and export capabilities
- **Modular Design**: Easy to customize and extend for specific research needs

## 🔧 Installation

### Prerequisites
- Python 3.8 or higher
- QuPath (for image analysis)
- Git

### Quick Install
```bash
# Clone the repository
git clone https://github.com/Xinyan-C/Spatialcell.git
cd Spatialcell

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
```

### Alternative: Install from PyPI (Coming Soon)
```bash
pip install spatialcell
```

## 📖 Quick Start

### 1. Basic Workflow
```python
from spatialcell.workflows import SpatialCellPipeline

# Initialize pipeline
pipeline = SpatialCellPipeline(
    sample_name="E18.5",
    input_dir="/path/to/visium/data",
    output_dir="/path/to/output"
)

# Run complete analysis
pipeline.run_full_analysis()
```

### 2. Step-by-Step Processing
```python
# Cell segmentation with Bin2cell
pipeline.run_segmentation(
    source_image_path="/path/to/histology.tif",
    roi_file="/path/to/regions.txt"
)

# Cell classification with TopAct  
pipeline.run_classification(
    classifier_path="/path/to/trained_model.joblib"
)

# Generate visualizations
pipeline.create_visualizations()
```

## 🗂️ Project Structure

```
Spatialcell/
├── spatialcell/                    # Main package
│   ├── qupath_scripts/             # QuPath cell detection scripts
│   ├── preprocessing/              # Data preprocessing modules
│   ├── spatial_segmentation/       # Bin2cell integration
│   ├── cell_annotation/            # TopAct classification
│   ├── utils/                      # Utility functions
│   └── workflows/                  # Complete pipelines
├── examples/                       # Usage examples and tutorials
├── requirements.txt                # Python dependencies
├── setup.py                       # Package installation
└── README.md                      # This file
```

## 📋 Supported Data Types

- **Spatial Transcriptomics**: 10x Visium, Slide-seq
- **Image Formats**: TIFF, SVG, PNG, JPEG
- **Development Stages**: E14.5, E18.5, P3 (extensible)
- **Cell Types**: Customizable classification schemes

## 🔬 Workflow Overview

1. **Histological Analysis**: QuPath-based nucleus detection
2. **Data Preprocessing**: SVG to NPZ conversion and filtering  
3. **Spatial Segmentation**: Bin2cell integration with cell boundaries
4. **Cell Classification**: TopAct machine learning annotation
5. **Visualization**: Multi-scale plotting and export

## 📚 Documentation

- **[Installation Guide](examples/installation.md)**: Detailed setup instructions
- **[Tutorial Notebooks](examples/)**: Step-by-step analysis examples
- **[API Reference](docs/)**: Complete function documentation
- **[FAQ](docs/faq.md)**: Common questions and troubleshooting

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Development Setup
```bash
# Clone and install in development mode
git clone https://github.com/Xinyan-C/Spatialcell.git
cd Spatialcell
pip install -e ".[dev]"

# Run tests
pytest tests/
```

## 📄 Citation

If you use SpatialCell in your research, please cite:

```bibtex
@software{spatialcell2024,
  author = {Xinyan},
  title = {SpatialCell: Integrated Spatial Transcriptomics Analysis Pipeline},
  url = {https://github.com/Xinyan-C/Spatialcell},
  year = {2024}
}
```

## 📧 Contact

- **Author**: Xinyan
- **Email**: keepandon@gmail.com
- **GitHub**: [@Xinyan-C](https://github.com/Xinyan-C)

## 📝 License

SpatialCell is licensed under the Apache License 2.0, which provides patent protection.

### Dependency Licenses:
- **bin2cell**: MIT License (automatically installed)
- **TopACT**: GPL v3 License (optional, user installs separately)

Apache 2.0 license includes patent protection clauses, providing additional legal protection for users.

## 📋 Patent Protection

The Apache 2.0 license provides:
- Patent grants from all contributors
- Patent retaliation protection
- Commercial use protection

For full license text, see [LICENSE](LICENSE) file.

## 🙏 Acknowledgments

- **QuPath**: For excellent histological image analysis tools
- **Bin2cell**: For spatial cell segmentation methods
- **TopAct**: For machine learning-based cell classification
- **Scanpy**: For single-cell analysis infrastructure