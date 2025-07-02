# SpatialCell: Integrated Spatial Transcriptomics Analysis Pipeline

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub issues](https://img.shields.io/github/issues/Xinyan-C/Spatialcell)](https://github.com/Xinyan-C/Spatialcell/issues)

**SpatialCell** is an integrated computational pipeline for spatial transcriptomics analysis that combines cell segmentation and automated cell type annotation. It seamlessly integrates **Stardist (applied as QuPath plugin for cell detection)** for histological image analysis, **Bin2cell** for spatial cell segmentation, and **TopAct** for machine learning-based cell classification.

## 🚀 Key Features

- **Multi-scale Cell Segmentation**: Stardist-enabled QuPath cell detection with Bin2cell spatial segmentation
- **Automated Cell Annotation**: TopAct-based machine learning classification
- **ROI-aware Processing**: Region-of-interest focused analysis for large datasets
- **Scalable Pipeline**: Support for multiple developmental time points (E14.5, E18.5, P3)
- **Visualization Tools**: Comprehensive plotting and export capabilities
- **Modular Design**: Easy to customize and extend for specific research needs

## 🔧 Installation

### Prerequisites
- Python 3.8 or higher (3.10 recommended)
- QuPath (for image analysis)
- Git


### Quick Install (Recommended)

```bash
pip install spatialcell
```
### Complete Functionality
```bash
pip install spatialcell
pip install git+https://gitlab.com/kfbenjamin/topact.git
```
### Alternative: Install from Source

```bash
# Clone the repository
git clone https://github.com/Xinyan-C/Spatialcell.git
cd Spatialcell

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
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
# Cell detection and spatial segmentation
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
│   ├── qupath_scripts/             # QuPath-Stardist integration scripts
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

1. **Histological Analysis**: Stardist-based cell detection via QuPath
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

## 🔗 Reference

- **QuPath**: Bankhead P, Loughrey MB, Fernández JA, et al. QuPath: Open source software for digital pathology image analysis. Sci Rep. 2017;7(1):16878. doi:10.1038/s41598-017-17204-5
- **Stardist**: Schmidt U, Weigert M, Broaddus C, Myers G. Cell detection with star-convex polygons. In: Frangi AF, Schnabel JA, Davatzikos C, Alberola-López C, Fichtinger G, eds. Medical Image Computing and Computer Assisted Intervention – MICCAI 2018. Springer International Publishing; 2018:265-273. doi:10.1007/978-3-030-00934-2_30
- **Bin2cell**: Polański K, Bartolomé-Casado R, Sarropoulos I, et al. Bin2cell reconstructs cells from high resolution visium HD data. Bioinformatics. 2024;40(9):btae546. doi:10.1093/bioinformatics/btae546
- **TopAct**: Benjamin K, Bhandari A, Kepple JD, et al. Multiscale topology classifies cells in subcellular spatial transcriptomics. Nature. 2024;630(8018):943-949. doi:10.1038/s41586-024-07563-1
- **Scanpy**: Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology. 2018;19(1):15. doi:10.1186/s13059-017-1382-0
