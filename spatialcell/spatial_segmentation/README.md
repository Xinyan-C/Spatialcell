# Spatial Segmentation Module

This module integrates bin2cell functionality for comprehensive spatial transcriptomics analysis, providing tools for cell segmentation, label expansion, and spatial visualization.

## Overview

The spatial segmentation module is built upon the [bin2cell package](https://github.com/Teichlab/bin2cell.git) and provides:

1. **Spatial data processing** with destriping and normalization
2. **Label integration** from QuPath nucleus detection results
3. **Label expansion** to capture surrounding spatial context
4. **Multi-modal segmentation** combining morphology and gene expression
5. **ROI-based analysis** for region-specific processing
6. **Comprehensive visualization** of segmentation results

## Key Features

### **Multi-Modal Segmentation**
- **Morphology-based**: Integration of QuPath/StarDist nucleus detection
- **Expression-based**: Gene expression-driven segmentation using StarDist
- **Joint segmentation**: Combines both approaches for optimal results

###  **Advanced Processing**
- **Destriping**: Removes spatial artifacts from Visium data
- **Label expansion**: Extends nucleus labels to capture cell boundaries
- **Spatial filtering**: ROI-based analysis for focused studies

###  **Rich Visualization**
- **Multi-panel plots**: Before/after comparisons
- **ROI-specific views**: Detailed visualization per region
- **Quality control**: Comprehensive QC plots and metrics

## Dependencies

### Core Dependencies
- **bin2cell**: Main spatial analysis framework
- **scanpy**: Single-cell analysis tools
- **pandas**: Data manipulation
- **numpy**: Numerical computations
- **matplotlib**: Plotting and visualization
- **scipy**: Scientific computing

### Visualization Dependencies
- **colorcet**: Advanced colormaps
- **matplotlib.patches**: Shape drawing

## Usage

### Command Line Interface

```bash
# Basic usage
python spatial_processor.py \
  --path /data/visium_sample/ \
  --source_image_path /data/tissue_image.tif \
  --region_file /data/roi_coordinates.txt \
  --npz_path /data/nucleus_labels.npz \
  --output_dir /data/results/ \
  --sample E14.5

# With custom parameters
python spatial_processor.py \
  --path /data/visium_sample/ \
  --source_image_path /data/tissue_image.tif \
  --region_file /data/roi_coordinates.txt \
  --npz_path /data/nucleus_labels.npz \
  --output_dir /data/results/ \
  --sample P3 \
  --algorithm volume_ratio \
  --volume_ratio 3.0 \
  --prob_thresh 0.1 \
  --labels_key labels_qupath_expanded
```

### Python API

```python
from spatial_segmentation import process_spatial_data
import argparse

# Create arguments object
class Args:
    def __init__(self):
        self.path = "/data/visium_sample/"
        self.source_image_path = "/data/tissue_image.tif"
        self.region_file = "/data/roi_coordinates.txt"
        self.npz_path = "/data/nucleus_labels.npz"
        self.output_dir = "/data/results/"
        self.sample = "E14.5"
        # ... other parameters

args = Args()

# Process spatial data
adata, cdata = process_spatial_data(args)
```

## Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--path` | Path to Visium data directory | `/data/visium/` |
| `--source_image_path` | Path to source tissue image | `/data/tissue.tif` |
| `--region_file` | Path to ROI coordinates file | `/data/regions.txt` |
| `--npz_path` | Path to QuPath NPZ labels | `/data/labels.npz` |
| `--output_dir` | Output directory | `/data/results/` |
| `--sample` | Sample identifier | `E14.5` |

### StarDist Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--prob_thresh` | 0.05 | Probability threshold for detection |
| `--nms_thresh` | 0.5 | Non-maximum suppression threshold |

### Label Expansion Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--algorithm` | `max_bin_distance` | Expansion algorithm |
| `--max_bin_distance` | 2 | Maximum distance for expansion |
| `--volume_ratio` | 4.0 | Volume ratio for expansion |
| `--k` | 4 | Number of nearest neighbors |
| `--subset_pca` | True | Use PCA for mean bins only |

### Analysis Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--labels_key` | `labels_joint` | Final labels for cell calling |

## Workflow Steps

### 1. Data Loading and Preprocessing
```python
# Load Visium data using bin2cell
adata = b2c.read_visium(path, source_image_path=source_image_path)

# Filter genes and cells
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_counts=1)

# Calculate microns per pixel (MPP)
mpp = adata.uns['spatial'][library]['scalefactors']['microns_per_pixel']
```

### 2. Image Processing
```python
# Scale H&E image
b2c.scaled_he_image(adata, mpp=mpp, save_path="stardist/he.tiff")

# Remove spatial artifacts
b2c.destripe(adata, adjust_counts=True)
```

### 3. Nucleus Segmentation Integration
```python
# Insert QuPath nucleus labels
b2c.insert_labels(adata, labels_npz_path=npz_path, 
                  basis="spatial", spatial_key="spatial", 
                  mpp=mpp, labels_key="labels_qupath")

# Expand labels to capture cell boundaries
b2c.expand_labels(adata, labels_key='labels_qupath',
                  expanded_labels_key="labels_qupath_expanded",
                  algorithm=algorithm, max_bin_distance=max_bin_distance)
```

### 4. Expression-Based Segmentation
```python
# Generate gene expression grid
b2c.grid_image(adata, "n_counts_adjusted", mpp=mpp,
               sigma=5, save_path="stardist/gex.tiff")

# Run StarDist on expression data
b2c.stardist(image_path="stardist/gex.tiff",
             labels_npz_path="stardist/gex_labels.npz",
             stardist_model="2D_versatile_fluo")

# Insert expression-based labels
b2c.insert_labels(adata, labels_npz_path="stardist/gex_labels.npz",
                  basis="array", mpp=mpp, labels_key="labels_gex")
```

### 5. Joint Segmentation
```python
# Combine morphology and expression segmentation
b2c.salvage_secondary_labels(adata, 
                            primary_label="labels_qupath_expanded",
                            secondary_label="labels_gex", 
                            labels_key="labels_joint")

# Convert bins to cells
cdata = b2c.bin_to_cell(adata, labels_key=labels_key,
                       spatial_keys=["spatial", "spatial_cropped_150_buffer"])
```

## Output Files

### Data Files
- `{sample}_2um.h5ad`: Processed spot-level data
- `{sample}_b2c.h5ad`: Cell-level data from bin-to-cell conversion
- `ROI_Data/{roi_name}/{roi_name}_adata.h5ad`: ROI-specific spot data
- `ROI_Data/{roi_name}/{roi_name}_cdata.h5ad`: ROI-specific cell data

### Intermediate Files
- `stardist/he.tiff`: Scaled H&E image
- `stardist/gex.tiff`: Gene expression grid image
- `stardist/gex_{prob}_{nms}.npz`: Expression-based segmentation labels

### Visualization Files

#### Quality Control
- `destripe/{roi_name}_destripe.pdf`: Destriping results
- `npz_labels/{roi_name}_npz_labels.pdf`: Initial label visualization

#### Segmentation Results
- `segmentation/{roi_name}_segmentation.pdf`: QuPath label integration
- `expanded_labels/{roi_name}_expanded_labels.pdf`: Expanded labels
- `gex_labels/{roi_name}_gex_labels.pdf`: Expression-based labels
- `joint_labels/{roi_name}_joint_labels.pdf`: Combined segmentation

#### Rendered Views
- `render_labels/{roi_name}_render_labels.pdf`: Rendered nucleus labels
- `render_gex/{roi_name}_render_gex.pdf`: Rendered expression labels

## ROI Processing

The module supports ROI-based analysis using coordinate files:

### ROI Coordinate Format
```
Sample - ROI1 Rectangle coordinates:
  X: 1000.0 - 2000.0
  Y: 1500.0 - 2500.0

Sample - ROI2 Rectangle coordinates:
  X: 2500.0 - 3500.0
  Y: 1000.0 - 2000.0
```

### ROI-Specific Outputs
Each ROI generates separate visualization and data files, enabling:
- **Focused analysis** on specific tissue regions
- **Comparative studies** between different areas
- **Quality control** per region
- **Independent processing** of large datasets

## Algorithm Options

### Label Expansion Algorithms

#### Max Bin Distance
```bash
--algorithm max_bin_distance --max_bin_distance 2
```
- Expands labels to nearby spots within specified distance
- Good for dense tissues with clear cell boundaries

#### Volume Ratio
```bash
--algorithm volume_ratio --volume_ratio 4.0
```
- Expands based on cell volume estimates
- Better for sparse tissues or when cell sizes vary

### StarDist Parameters

#### Probability Threshold
- **Higher values** (0.1-0.3): More stringent detection, fewer false positives
- **Lower values** (0.01-0.05): More sensitive detection, may include noise

#### NMS Threshold
- **Higher values** (0.7-0.9): Allows closer detections, may merge nearby cells
- **Lower values** (0.3-0.5): Stricter separation, reduces over-segmentation

## Quality Control

### Automated QC Metrics
- Spot count before/after filtering
- Label distribution across ROIs
- Segmentation coverage statistics
- Cell size distributions

### Visual QC
- Multi-panel comparisons
- ROI-specific visualizations  
- Overlay validation plots
- Before/after processing views

## Integration with Pipeline

### Input Requirements
1. **Visium data**: Standard Space Ranger output
2. **Tissue image**: High-resolution H&E image
3. **Nucleus labels**: NPZ file from QuPath processing
4. **ROI coordinates**: Standardized coordinate file

### Output for Downstream Analysis
- **Cell-level data** ready for cell type annotation
- **Spatial coordinates** preserved for spatial analysis
- **Quality metrics** for method validation
- **Visualization files** for publication

## Troubleshooting

### Common Issues

1. **Memory errors**: Reduce image resolution or process ROIs separately
2. **Segmentation quality**: Adjust StarDist parameters
3. **Label expansion issues**: Try different algorithms or parameters
4. **ROI parsing errors**: Check coordinate file format

### Performance Optimization

- **Parallel processing**: Multiple ROIs processed independently
- **Memory management**: Efficient handling of large images
- **Intermediate caching**: Saves processing time for reruns
- **Modular design**: Process individual steps separately

## Citation

If you use this module, please cite:

1. **Spatialcell pipeline**: [Your publication]
2. **bin2cell**: Kueckelhaus et al. (2024)
3. **scanpy**: Wolf et al. (2018)
4. **StarDist**: Schmidt et al. (2018)

## Examples

See the `examples/` directory for complete workflow examples and sample data.