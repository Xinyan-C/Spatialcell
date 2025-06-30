# QuPath Scripts

This directory contains Groovy scripts for QuPath software to perform nucleus detection using StarDist.

## Prerequisites

1. **QuPath Software**: Download and install QuPath (https://qupath.github.io/)
2. **StarDist Extension**: Install the StarDist extension for QuPath
3. **StarDist Model**: Download a pre-trained StarDist model (e.g., `he_heavy_augment.pb`)

## Usage

### 1. Nucleus Detection

**File**: `nucleus_detection.groovy`

**Purpose**: Detect cell nuclei in H&E stained tissue images using StarDist

**Steps**:
1. Open your H&E image in QuPath
2. Create annotations for regions of interest (ROI)
3. Select the annotations you want to process
4. Update the `modelPath` variable in the script to point to your StarDist model
5. Run the script in QuPath's script editor
6. Export the detection results as SVG files for downstream processing

**Parameters to adjust**:
- `modelPath`: Path to your StarDist model file
- `threshold`: Probability threshold for detection (default: 0.3)
- `pixelSize`: Detection resolution (default: 0.3)

## Output

The script will detect nuclei within selected annotations and display the total count. The detected objects can then be exported as SVG files for further processing in the Spatialcell pipeline.

## Citations

If you use these scripts, please cite:
- StarDist paper: https://doi.org/10.48550/arXiv.1806.03535
- QuPath paper: https://doi.org/10.1038/s41598-017-17204-5
- Spatialcell pipeline: https://github.com/Xinyan-C/Spatialcell