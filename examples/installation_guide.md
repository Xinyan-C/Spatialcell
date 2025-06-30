# SpatialCell Installation Guide

## System Requirements

- **Python**: 3.10 or higher
- **Memory**: At least 16GB RAM (32GB+ recommended for large datasets)
- **Storage**: 10GB+ free space for intermediate files
- **OS**: Linux, macOS, or Windows

## Step-by-Step Installation

### 1. Install Python Environment

#### Option A: Using Conda (Recommended)
```bash
# Create new environment
conda create -n spatialcell python=3.10
conda activate spatialcell
```

#### Option B: Using venv
```bash
python -m venv spatialcell_env
source spatialcell_env/bin/activate  # Linux/Mac
# or
spatialcell_env\Scripts\activate     # Windows
```

### 2. Install SpatialCell

#### From GitHub (Latest)
```bash
git clone https://github.com/Xinyan-C/Spatialcell.git
cd Spatialcell
pip install -r requirements.txt
pip install -e .
```

#### From PyPI (Coming Soon)
```bash
pip install spatialcell
```

### 3. Install QuPath (Required)

1. Download QuPath from: https://qupath.github.io/
2. Install following the official instructions
3. Note the installation path for later use

### 4. Verify Installation

```python
# Test import
import spatialcell
print(f"SpatialCell version: {spatialcell.__version__}")

# Test core modules
from spatialcell.workflows import SpatialCellPipeline
from spatialcell.utils import load_config
print("✅ Installation successful!")
```

## Troubleshooting

### Common Issues

#### 1. Missing Dependencies
```bash
# If you encounter missing packages
pip install --upgrade pip
pip install -r requirements.txt --force-reinstall
```

#### 2. CUDA/GPU Issues
```bash
# For GPU acceleration (optional)
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118
```

#### 3. Memory Issues
```bash
# Increase swap space or use smaller batch sizes
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
```

#### 4. QuPath Integration
- Ensure QuPath is properly installed
- Check that groovy scripts can access your data paths
- Verify image file formats are supported

### Platform-Specific Notes

#### Linux
```bash
# Install system dependencies if needed
sudo apt-get update
sudo apt-get install build-essential python3-dev
```

#### macOS
```bash
# Install Xcode command line tools if needed
xcode-select --install
```

#### Windows
- Use Anaconda/Miniconda for easier dependency management
- Consider using WSL2 for Linux-like environment

### Getting Help

1. **Check logs**: Look for error messages in console output
2. **GitHub Issues**: https://github.com/Xinyan-C/Spatialcell/issues
3. **Email**: keepandon@gmail.com

## Next Steps

After successful installation:

1. Copy `examples/config_example.yml` to your working directory
2. Modify paths in the configuration file
3. Run `python examples/basic_tutorial.py`
4. Follow the tutorial notebooks for detailed examples

## Development Installation

For contributing to SpatialCell:

```bash
git clone https://github.com/Xinyan-C/Spatialcell.git
cd Spatialcell
pip install -e ".[dev]"

# Run tests
pytest tests/

# Format code
black spatialcell/
flake8 spatialcell/
```