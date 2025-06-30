# Utilities Module

This module provides essential utility tools for the Spatialcell pipeline, focusing on ROI (Region of Interest) coordinate extraction and validation.

## Overview

The utilities module supports the spatial transcriptomics workflow by providing tools to:

1. **Extract ROI coordinates** from Loupe Browser exported CSV files
2. **Validate ROI coordinates** by visualizing them on source images
3. **Generate standardized coordinate files** for downstream analysis

## Tools

### 1. ROI Extractor (`roi_extractor.py`)

Extracts ROI coordinates from Loupe Browser exported CSV files and generates standardized coordinate range files.

#### Usage

**Command Line:**
```bash
# Basic usage
python roi_extractor.py --sample E14.5

# With custom directory
python roi_extractor.py --sample P3 --base_dir /data/coordinates/

# With verbose output
python roi_extractor.py --sample custom_sample --base_dir /data/ --verbose
```

**Python API:**
```python
from utils import extract_roi_coordinates

# Extract coordinates
roi_ranges = extract_roi_coordinates(
    sample_name="E14.5",
    sample_dir="/data/E14.5/",
    output_path="/data/E14.5_ranges.txt"
)
```

#### Input File Structure

The tool expects the following file structure from Loupe Browser exports:

```
sample_directory/
├── sample_name-all.csv      # Total coordinates (required)
├── sample_name-CS1.csv      # ROI 1 (optional)
├── sample_name-CS2.csv      # ROI 2 (optional)
├── sample_name-WT1.csv      # ROI 3 (optional)
└── ...                      # Additional ROI files
```

#### Output Format

Generates a standardized coordinate range file:
```
Sample - ROI1 Rectangle coordinates:
  X: 1000.0 - 2000.0
  Y: 1500.0 - 2500.0

Sample - ROI2 Rectangle coordinates:
  X: 2500.0 - 3500.0
  Y: 1000.0 - 2000.0
```

### 2. ROI Validator (`roi_validator.py`)

Validates ROI coordinates by overlaying them on source images and generating validation reports.

#### Usage

**Command Line:**
```bash
# Basic validation
python roi_validator.py --image_path tissue.tif --coords_file E14.5_ranges.txt --sample E14.5

# With custom output directory
python roi_validator.py \
  --image_path /data/images/tissue.tif \
  --coords_file /data/coordinates/P3_ranges.txt \
  --sample P3 \
  --output_dir /data/validation/
```

**Python API:**
```python
from utils import parse_coordinate_file, plot_image_with_rois
from utils.roi_validator import read_image

# Load and validate
image = read_image("tissue.tif")
rois = parse_coordinate_file("sample_ranges.txt")
plot_image_with_rois(image, rois, "validation_plot.png")
```

#### Output Files

1. **Validation Plot**: Visual overlay of ROIs on the source image
2. **Validation Report**: Detailed text report with coordinate validation

## Integration with Pipeline

### Workflow Integration

1. **Export from Loupe Browser**: Select ROIs and export as CSV files
2. **Extract Coordinates**: Use `roi_extractor.py` to generate coordinate ranges
3. **Validate Coordinates**: Use `roi_validator.py` to verify ROI placement
4. **Use in Pipeline**: Feed coordinate files to spatial processing modules

### Example Workflow

```bash
# Step 1: Extract ROI coordinates
python utils/roi_extractor.py --sample E14.5 --base_dir /data/loupe_exports/

# Step 2: Validate coordinates
python utils/roi_validator.py \
  --image_path /data/images/E14.5_tissue.tif \
  --coords_file /data/loupe_exports/E14.5/E14.5_ranges.txt \
  --sample E14.5

# Step 3: Use in spatial processing
python spatial_segmentation/spatial_processor.py \
  --path /data/visium/E14.5/ \
  --region_file /data/loupe_exports/E14.5/E14.5_ranges.txt \
  ...
```

## File Format Support

### Input Formats
- **CSV files**: Loupe Browser exports
- **Image files**: TIFF, PNG, JPEG (via tifffile or Pillow)

### Output Formats
- **Coordinate files**: Plain text with standardized format
- **Validation plots**: PNG images with ROI overlays
- **Reports**: Text files with validation details

## Dependencies

### Required
- pandas
- numpy
- matplotlib

### Optional (auto-detected)
- tifffile (preferred for TIFF files)
- Pillow (fallback for image reading)

## Error Handling

The tools include comprehensive error handling for:
- Missing input files
- Invalid coordinate formats
- Image reading errors
- Coordinate validation issues

## Configuration

### ROI Detection Patterns

The extractor automatically detects common ROI file patterns:
- `sample-CS1.csv`, `sample-CS2.csv`, `sample-CS3.csv`
- `sample-WT1.csv`, `sample-WT2.csv`, `sample-WT3.csv`
- Custom patterns: `sample-[label].csv`

### Validation Criteria

The validator checks for:
- Coordinates within image bounds
- Valid coordinate ranges (min < max)
- Reasonable ROI sizes (> 10×10 pixels)
- Coordinate data integrity

## Troubleshooting

### Common Issues

1. **No ROI files found**: Check file naming conventions
2. **Coordinate validation warnings**: Verify ROI selections in Loupe Browser
3. **Image loading errors**: Ensure tifffile or Pillow is installed
4. **Empty coordinate ranges**: Check CSV file format and barcode matching

### Debug Mode

Enable verbose logging for detailed troubleshooting:
```bash
python roi_extractor.py --sample E14.5 --verbose
python roi_validator.py --image_path tissue.tif --coords_file ranges.txt --sample E14.5 --verbose
```

## API Reference

### Key Functions

```python
# ROI Extraction
extract_roi_coordinates(sample_name, sample_dir, output_path=None)
calculate_roi_range(roi_csv_path, total_coords, sample_name, roi_label)
detect_roi_files(sample_dir, sample_name)

# ROI Validation  
parse_coordinate_file(coords_file)
plot_image_with_rois(image, rois, output_path=None)
validate_roi_coordinates(rois, image_shape)
```

See individual function docstrings for detailed parameter descriptions.