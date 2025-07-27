# Cell Annotation Module

This module provides comprehensive tools for spatial cell type annotation using the TopAct framework, enabling accurate cell type prediction in spatial transcriptomics data with cell boundary constraints from bin2cell segmentation.

## Overview

The cell annotation module integrates TopAct's powerful classification capabilities with Spatialcell's enhanced spatial analysis framework to provide:

1. **Reference-based classifier training** from single-cell RNA-seq data
2. **High-resolution spatial annotation** with floating-point coordinate support
3. **Cell boundary-constrained analysis** using bin2cell segmentation results
4. **Multi-scale classification** across different neighborhood sizes
5. **Comprehensive visualization** with publication-quality outputs

## Architecture

### Integration Framework

```
Single-cell Reference Data → TopAct Classifier Training
                                      ↓
Spatial Transcriptomics Data → Cell Boundary Segmentation (bin2cell)
                                      ↓
High-Definition Spatial Data → TopAct Classification + Cell Constraints
                                      ↓
Multi-scale Annotation Results → Comprehensive Visualization
```

## Components

### 1. Classifier Training (`classifier_trainer.py`)

Trains TopAct Support Vector Classifiers from single-cell reference data, supporting time-point-specific model generation.

#### Features
- **Multi-format data loading** from R-exported matrix market files
- **Time-point stratification** for developmental studies
- **Quality validation** with comprehensive error checking
- **Batch processing** for multiple time points
- **Model serialization** for downstream analysis

#### Usage

**Command Line:**
```bash
# Train classifiers for multiple time points
python classifier_trainer.py \
  --time_points E14.5 E18.5 P3 \
  --input_dir /data/reference_data/ \
  --output_dir /data/classifiers/ \
  --label_column celltype_merge

# Train with custom cell type labels
python classifier_trainer.py \
  --time_points P3 \
  --input_dir /data/reference/ \
  --output_dir /data/models/ \
  --label_column custom_annotation
```

**Python API:**
```python
from cell_annotation import train_multiple_time_points

# Train classifiers
classifier_paths = train_multiple_time_points(
    time_points=['E14.5', 'E18.5'],
    input_dir='/data/reference/',
    output_dir='/data/classifiers/',
    label_column='celltype_merge'
)
```

#### Input Requirements

The training module expects R-exported data with specific file structure:

```
reference_data/
├── all_samples_counts.mtx      # Sparse count matrix (genes x cells)
├── all_samples_genes.csv       # Gene identifiers
├── all_samples_barcodes.csv    # Cell barcodes  
├── all_samples_meta_data.csv   # Cell metadata with annotations
```

**File Formats:**
- **MTX**: Matrix Market format sparse matrix
- **CSV**: Standard comma-separated values with headers
- **Metadata**: Must contain `orig.ident` and specified label column

#### Output Files
- `clf_{time_point}.joblib`: Serialized TopAct classifier
- `{time_point}_TopAct_Classifier_Training.log`: Training logs
- Performance metrics and validation results

### 2. Spatial Annotation (`annotation_processor.py`)

Processes high-definition spatial transcriptomics data for cell type annotation using trained classifiers with cell boundary constraints.

#### Key Innovations
- **Cell boundary awareness** using bin2cell segmentation
- **Floating-point coordinate support** for high-resolution data
- **Multi-scale neighborhood analysis** with micrometer precision
- **Parallel processing** for computational efficiency
- **ROI-based analysis** for region-specific studies

#### Usage

**Command Line:**
```bash
# Basic spatial annotation
python annotation_processor.py \
  --sample E14.5 \
  --out_dir /data/results/ \
  --expr_path /data/spatial/matrix.h5 \
  --pos_path /data/spatial/positions.parquet \
  --source_image_path /data/tissue.tif \
  --roi_file /data/regions.txt \
  --clf_path /data/classifiers/clf_E14.5.joblib \
  --bin2cell_dir /data/segmentation.h5ad

# Advanced processing with custom parameters
python annotation_processor.py \
  --sample P3 \
  --out_dir /data/results/ \
  --expr_path /data/spatial/matrix.h5 \
  --pos_path /data/spatial/positions.parquet \
  --source_image_path /data/tissue.tif \
  --roi_file /data/regions.txt \
  --clf_path /data/classifiers/clf_P3.joblib \
  --bin2cell_dir /data/segmentation.h5ad \
  --min_scale 2.0 \
  --max_scale 12.0 \
  --num_proc 120 \
  --labels labels_joint
```

#### Processing Workflow

1. **Data Integration**
   - Load HD spatial expression data
   - Import bin2cell segmentation results
   - Parse ROI coordinate definitions
   - Extract classifier training genes

2. **Cell Label Assignment**
   - Map HD coordinates to cell boundaries
   - Assign cell labels from bin2cell results
   - Filter data within ROI boundaries
   - Prepare expression matrices

3. **Multi-scale Classification**
   - Build spatial count grids
   - Apply cell-constrained neighborhoods
   - Perform parallel TopAct classification
   - Generate confidence matrices

4. **Result Generation**
   - Save classification probabilities
   - Cache intermediate results
   - Generate processing logs

### 3. Visualization (`annotation_visualizer.py`)

Creates comprehensive visualizations of spatial classification results with publication-quality outputs.

#### Visualization Types

**Scatter Plot Visualization**
- Pure classification results overlay
- Cell type legend with custom colors
- High-resolution PDF output

**Background Image Integration**
- Tissue image with classification overlay
- Side-by-side comparison views
- ROI-specific cropping and scaling

**Data Export**
- CSV files with coordinates and predictions
- Confidence scores and metadata
- Ready for further analysis

#### Usage

**Command Line:**
```bash
# Basic visualization
python annotation_visualizer.py \
  --sample E14.5 \
  --sd_dir /data/results/ \
  --outfile_dir /data/results/ \
  --clf_dir /data/classifiers/ \
  --roi_file /data/E14.5_ranges.txt

# Advanced visualization with background
python annotation_visualizer.py \
  --sample P3 \
  --sd_dir /data/results/ \
  --outfile_dir /data/results/ \
  --clf_dir /data/classifiers/ \
  --roi_file /data/P3_ranges.txt \
  --background_image /data/tissue.tif \
  --color_scheme functional \
  --rename_cell_types True \
  --point_size 15 \
  --verbose
```

#### Color Schemes

Multiple scientifically-designed color schemes:

- **primary**: Default scheme optimized for bone/cartilage studies
- **scientific**: High contrast, publication-ready
- **functional**: Colors grouped by biological function
- **modern**: Contemporary gradient aesthetics
- **warm**: Comfortable for extended viewing
- **golden**: Mathematical golden ratio spacing

#### Output Files

**Per ROI:**
- `Spatial_Classification_{ROI}.pdf`: Scatter plot visualization
- `Spatial_Classification_{ROI}_overlay.pdf`: Background + overlay
- `Spatial_Classification_{ROI}_side_by_side.pdf`: Comparison view
- `Spatial_Classification_{ROI}_data.csv`: Classification data

## Integration with Pipeline

### Workflow Integration

```bash
# Step 1: Train classifiers from reference data
python cell_annotation/classifier_trainer.py \
  --time_points E14.5 E18.5 P3 \
  --input_dir /data/reference/ \
  --output_dir /data/classifiers/

# Step 2: Process spatial data with bin2cell constraints  
python cell_annotation/annotation_processor.py \
  --sample E14.5 \
  --out_dir /data/annotation/ \
  --expr_path /data/visium/matrix.h5 \
  --pos_path /data/visium/positions.parquet \
  --source_image_path /data/tissue.tif \
  --roi_file /data/regions.txt \
  --clf_path /data/classifiers/clf_E14.5.joblib \
  --bin2cell_dir /data/segmentation.h5ad

# Step 3: Generate comprehensive visualizations
python cell_annotation/annotation_visualizer.py \
  --sample E14.5 \
  --sd_dir /data/annotation/ \
  --outfile_dir /data/annotation/ \
  --clf_dir /data/classifiers/ \
  --roi_file /data/regions.txt \
  --background_image /data/tissue.tif
```

### Python API Integration

```python
import scanpy as sc
from cell_annotation import (
    train_topact_classifier, 
    process_sample_annotation,
    load_classification_data
)

# Train classifier
adata_ref = sc.read_h5ad('reference.h5ad')
classifier = train_topact_classifier(adata_ref, 'E14.5')

# Process spatial data
args = create_processing_args()  # Configure parameters
process_sample_annotation(args)

# Load and analyze results
sd, confidence_matrix, classes = load_classification_data(
    'sd_file.joblib', 'results.npy', 'classifier.joblib'
)
```

## Advanced Features

### Cell Boundary Constraints

The module implements custom neighborhood functions that respect cell boundaries from bin2cell segmentation:

```python
def create_cell_constrained_neighborhood(sd, center, scale):
    """Custom neighborhood constrained to single cells."""
    center_x, center_y = center
    center_label = get_cell_label(sd, center_x, center_y)
    
    if center_label == 0:
        return pd.DataFrame()  # Background
    
    # Filter to same cell within scale
    same_cell_data = sd.table[sd.table['cell_label'] == center_label]
    distances = calculate_distances(same_cell_data, center_x, center_y)
    return same_cell_data[distances <= scale]
```

### Multi-Scale Analysis

Classification performed across multiple spatial scales:

- **Scale Range**: 3-9 micrometers (default)
- **Scale Steps**: 2 micrometer increments
- **Physical Units**: Automatic conversion using MPP (microns per pixel)
- **Adaptive Classification**: Starts from first non-zero expression scale

### Parallel Processing

Efficient parallel implementation:

- **Process Pool**: Configurable number of worker processes
- **Memory Management**: Optimized for large datasets
- **Progress Tracking**: Verbose logging and status updates
- **Error Handling**: Robust recovery from individual failures




## Troubleshooting

### Common Issues

**Memory Errors**
- Reduce number of parallel processes
- Process ROIs individually
- Use smaller scale ranges
- Implement data chunking

**Classification Quality**
- Verify reference data quality
- Check gene overlap between reference and spatial
- Adjust confidence thresholds
- Validate cell boundary assignments

**Performance Issues**
- Enable parallel processing
- Use SSD storage for intermediate files
- Optimize scale parameters
- Cache intermediate results

### Debug Mode

Enable detailed logging:
```bash
python annotation_processor.py --verbose [other args]
python annotation_visualizer.py --verbose [other args]
```

## Dependencies

### Core Requirements
- **TopAct**: Classification framework
- **bin2cell**: Spatial segmentation
- **scanpy**: Single-cell analysis
- **pandas**: Data manipulation
- **numpy**: Numerical computing
- **scipy**: Scientific computing
- **scikit-learn**: Machine learning utilities

### Visualization Requirements
- **matplotlib**: Plotting framework
- **PIL/Pillow**: Image processing
- **natsort**: Natural sorting
- **colormaps**: Enhanced color palettes

### Optional Dependencies
- **joblib**: Model serialization
- **tifffile**: TIFF image support
- **dask**: Large dataset processing
- **numba**: Performance acceleration


## Citation

If you use this module, please cite:

1. **Spatialcell pipeline**: [Hold on to my publication]
2. **TopAct framework**: Andreatta et al. (2021)
3. **bin2cell package**: Kueckelhaus et al. (2024)
4. **scanpy**: Wolf et al. (2018)

## Examples

See the `examples/` directory for:
- Complete workflow tutorials
- Sample datasets and configurations
- Best practice guidelines
- Troubleshooting examples

## Contributing

Guidelines for extending the cell annotation module:
- Follow TopAct API conventions
- Maintain cell boundary constraint compatibility
- Include comprehensive documentation
- Add appropriate error handling
- Provide usage examples