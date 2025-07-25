# Preprocessing Module

This module handles the conversion of QuPath-exported SVG files to NPZ format for downstream spatial transcriptomics analysis.

## Overview

After nucleus detection in QuPath, the detected objects need to be exported as SVG files and converted to a format suitable for spatial analysis. This module provides tools to:

1. Parse QuPath-exported SVG files
2. Convert SVG paths/polygons to binary label matrices
3. Save results as compressed NPZ files for efficient storage

## Files

- `svg_to_npz.py`: Main conversion script
- `__init__.py`: Python package initialization

## Usage

### Command Line Interface

```bash
python svg_to_npz.py --svg nuclei.svg --height 2048 --width 2048 --output labels.npz
```

### Parameters

- `--svg`: Path to QuPath-exported SVG file
- `--height`: Original image height in pixels
- `--width`: Original image width in pixels  
- `--output`: Output NPZ file path
- `--verbose`: Enable detailed logging (optional)

### Example Workflow

1. **In QuPath**: Run nucleus detection script and export results as SVG
2. **Get image dimensions**: Note the original image height and width
3. **Convert to NPZ**:
   ```bash
   python preprocessing/svg_to_npz.py \
     --svg /path/to/nuclei_detections.svg \
     --height 4096 \
     --width 4096 \
     --output /path/to/labels.npz
   ```

### Python API

```python
from preprocessing.svg_to_npz import convert_svg_to_npz

# Convert SVG to NPZ
num_objects = convert_svg_to_npz(
    svg_path="nuclei.svg",
    height=2048,
    width=2048,
    output_path="labels.npz",
    verbose=True
)

print(f"Converted {num_objects} nucleus objects")
```

## Input Format

The script expects SVG files exported from QuPath containing:
- `<path>` elements (preferred) with `d` attributes containing coordinate data
- OR `<polygon>` elements with `points` attributes
- Each element represents one detected nucleus

## Output Format

- **NPZ file**: Compressed sparse matrix format
- **Label matrix**: 2D array where each pixel has a label value
  - `0`: Background
  - `1, 2, 3, ...`: Individual nucleus labels

## Dependencies

- numpy
- opencv-python (cv2)
- scipy
- lxml (for XML parsing)

## Error Handling

The script includes robust error handling for:
- Invalid SVG files
- Missing namespace declarations
- Malformed coordinate data
- File I/O errors

## Performance Notes

- Uses sparse matrix format for memory efficiency
- Supports large images (tested up to 4096x4096 pixels)
- Processing time scales with number of detected objects