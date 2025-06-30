# Workflows Module

This module provides high-level workflow orchestration for the complete Spatialcell analysis pipeline, enabling both single sample processing and efficient batch operations across multiple samples.

## Overview

The workflows module serves as the orchestration layer that integrates all Spatialcell components into cohesive, end-to-end analysis pipelines. It provides:

1. **Complete pipeline execution** for single samples
2. **Batch processing capabilities** for multiple samples  
3. **Configuration management** with validation
4. **Progress tracking and reporting** for all operations
5. **Error handling and recovery** mechanisms

## Components

### 1. Complete Pipeline (`complete_pipeline.py`)

Orchestrates the entire Spatialcell workflow from start to finish for individual samples.

#### Workflow Steps

1. **ROI Coordinate Extraction** (optional)
   - Extract regions of interest from Loupe Browser exports
   - Generate standardized coordinate files

2. **SVG to NPZ Conversion** (optional)
   - Convert QuPath nucleus detection results
   - Create binary label matrices

3. **Spatial Segmentation**
   - Process Visium data with bin2cell integration
   - Perform nucleus detection and label expansion
   - Generate cell-level aggregated data

4. **Cell Type Annotation** (optional)
   - Apply TopAct classifiers to spatial data
   - Multi-scale neighborhood analysis
   - Cell boundary-constrained classification

5. **Comprehensive Visualization** (optional)
   - Generate publication-quality figures
   - Multiple visualization formats
   - ROI-specific analysis outputs

#### Usage

**Command Line:**
```bash
# Run complete pipeline from configuration
python complete_pipeline.py --config sample_config.yaml

# Validate configuration only
python complete_pipeline.py --config sample_config.yaml --validate_only

# Override output directory
python complete_pipeline.py --config sample_config.yaml --output_dir /custom/output/
```

**Python API:**
```python
from workflows import run_complete_pipeline

# Execute pipeline
run_complete_pipeline('sample_config.yaml')
```

#### Configuration File Format

The pipeline uses YAML configuration files with comprehensive parameter control:

```yaml
# Basic sample information
sample_name: E14_5_sample
output_dir: ./spatialcell_output

# Pipeline step controls
enable_roi_extraction: true
enable_svg_conversion: true  
enable_spatial_segmentation: true
enable_cell_annotation: true
enable_visualization: true

# Input data paths
visium_data_path: ./data/visium/
source_image_path: ./data/tissue.tif
svg_file: ./data/nuclei.svg
classifier_path: ./models/clf_E14.5.joblib

# Processing parameters
stardist_prob_thresh: 0.05
min_scale: 3.0
max_scale: 9.0
num_processes: 80

# Visualization settings
color_scheme: primary
point_size: 10
rename_cell_types: true
```

### 2. Configuration Management (`pipeline_config.py`)

Provides robust configuration management with validation and default value handling.

#### PipelineConfig Class

Comprehensive configuration class with:
- **Dataclass-based structure** with type hints
- **Automatic path derivation** based on sample names
- **Input validation** with detailed error messages
- **File format support** for YAML and JSON
- **Default value management** for all parameters

#### Key Features

**Configuration Loading:**
```python
from workflows import PipelineConfig

# Load from file
config = PipelineConfig.from_file('config.yaml')

# Create from dictionary
config = PipelineConfig.from_dict(config_dict)

# Create with defaults
config = PipelineConfig(sample_name='test', output_dir='./output')
```

**Validation:**
```python
from workflows import validate_pipeline_inputs

# Comprehensive validation
validate_pipeline_inputs(config)
```

**Configuration Generation:**
```python
from workflows import create_example_config

# Generate example configuration
create_example_config('example_config.yaml', 'sample_name')
```

#### Configuration Parameters

**File Paths:**
- Input data locations (Visium, images, models)
- Output directory structure
- Intermediate file paths

**Processing Parameters:**
- StarDist detection thresholds
- Label expansion algorithms
- Multi-scale analysis ranges
- Computational resource limits

**Visualization Options:**
- Color schemes and palettes
- Point sizes and shapes
- Background image integration
- Output formats and quality

### 3. Batch Processing (`batch_processor.py`)

Enables efficient processing of multiple samples with parallel execution and comprehensive reporting.

#### Features

**Parallel Processing:**
- Configurable number of worker processes
- Automatic load balancing
- Memory-efficient execution

**Error Handling:**
- Continue-on-error options
- Failed configuration preservation
- Detailed error logging and reporting

**Progress Tracking:**
- Real-time status updates
- Processing time estimates
- Success/failure statistics

#### Usage

**Command Line:**
```bash
# Basic batch processing
python batch_processor.py \
  --config_dir ./sample_configs/ \
  --output_dir ./batch_output/

# Parallel processing with 8 workers
python batch_processor.py \
  --config_dir ./configs/ \
  --output_dir ./output/ \
  --max_parallel 8

# Validation only mode
python batch_processor.py \
  --config_dir ./configs/ \
  --validate_only

# Continue on errors
python batch_processor.py \
  --config_dir ./configs/ \
  --output_dir ./output/ \
  --continue_on_error
```

**Python API:**
```python
from workflows import BatchConfig, run_batch_processing

# Create batch configuration
batch_config = BatchConfig(
    base_output_dir='./batch_output',
    max_parallel=4,
    continue_on_error=True
)

# Run batch processing
summary = run_batch_processing('./configs/', batch_config)
```

#### Batch Configuration Structure

```
sample_configs/
├── sample_E14_5_config.yaml
├── sample_E18_5_config.yaml
├── sample_P3_config.yaml
└── template_config.yaml
```

#### Batch Output Structure

```
batch_output/
├── samples/
│   ├── E14_5/
│   │   ├── spatial_segmentation/
│   │   ├── cell_annotation/
│   │   └── visualizations/
│   ├── E18_5/
│   └── P3/
├── batch_logs/
├── batch_reports/
└── failed_configs/
```

## Advanced Features

### Configuration Templates

Create multiple sample configurations from templates:

```python
from workflows import create_batch_configs_from_template

# Sample-specific parameters
sample_info = [
    {
        'sample_name': 'E14_5',
        'visium_data_path': './data/E14_5/visium/',
        'classifier_path': './models/clf_E14.5.joblib'
    },
    {
        'sample_name': 'E18_5', 
        'visium_data_path': './data/E18_5/visium/',
        'classifier_path': './models/clf_E18.5.joblib'
    }
]

# Generate configurations
config_files = create_batch_configs_from_template(
    'template.yaml', sample_info, './generated_configs/'
)
```

### Custom Workflow Steps

Disable specific pipeline steps for focused analysis:

```yaml
# Only run spatial segmentation
enable_roi_extraction: false
enable_svg_conversion: false
enable_spatial_segmentation: true
enable_cell_annotation: false
enable_visualization: false
```

### Resource Management

Optimize computational resources:

```yaml
# High-performance configuration
num_processes: 120
num_threads: 120
memory_limit_gb: 800

# Memory-constrained configuration  
num_processes: 40
num_threads: 40
memory_limit_gb: 200
```

### Error Recovery

Handle processing failures gracefully:

```bash
# Continue batch processing despite failures
python batch_processor.py \
  --config_dir ./configs/ \
  --output_dir ./output/ \
  --continue_on_error

# Review failed configurations
ls ./output/failed_configs/
```

## Integration Examples

### Single Sample Analysis

```python
from workflows import PipelineConfig, run_complete_pipeline

# Create configuration
config = PipelineConfig(
    sample_name='developmental_sample',
    visium_data_path='./data/visium/',
    source_image_path='./data/tissue.tif',
    classifier_path='./models/classifier.joblib',
    output_dir='./results/'
)

# Save configuration
config.save('sample_config.yaml')

# Run pipeline
run_complete_pipeline('sample_config.yaml')
```

### Batch Analysis

```python
from workflows import BatchConfig, run_batch_processing

# Setup batch processing
batch_config = BatchConfig(
    base_output_dir='./developmental_study/',
    max_parallel=8,
    continue_on_error=True
)

# Process all samples
summary = run_batch_processing('./sample_configs/', batch_config)

# Review results
print(f"Success rate: {summary['success_rate']:.1%}")
for result in summary['processing_results']:
    if result['success']:
        print(f"✓ {result['sample_name']}: {result['processing_time']:.1f}s")
    else:
        print(f"✗ {result['sample_name']}: {result['error_message']}")
```

### Development Workflow

```python
from workflows import PipelineConfig, validate_pipeline_inputs

# Validate configuration during development
config = PipelineConfig.from_file('test_config.yaml')

try:
    validate_pipeline_inputs(config)
    print("Configuration is valid!")
except ValueError as e:
    print(f"Configuration error: {e}")
    
# Generate example for reference
create_example_config('reference_config.yaml', 'example_sample')
```

## Performance Optimization

### Computational Resources

**CPU Optimization:**
- Match `num_processes` to available CPU cores
- Use `num_threads` for memory-bound operations
- Consider NUMA topology for large systems

**Memory Management:**
- Set `memory_limit_gb` based on available RAM
- Monitor peak memory usage during processing
- Use batch processing for memory-constrained systems

**Storage Optimization:**
- Use SSD storage for intermediate files
- Configure appropriate temp directories
- Monitor disk space during batch operations

### Batch Processing Strategies

**Small Datasets (< 10 samples):**
- Use sequential processing (`max_parallel=1`)
- Enable all visualization options
- Generate comprehensive reports

**Medium Datasets (10-50 samples):**
- Use moderate parallelism (`max_parallel=4-8`)
- Balance resource usage across samples
- Enable progress monitoring

**Large Datasets (> 50 samples):**
- Use high parallelism (`max_parallel=16+`)
- Implement checkpointing strategies
- Use distributed processing if available

### Error Handling Strategies

**Development Phase:**
- Disable `continue_on_error` for immediate feedback
- Use `validate_only` mode for configuration testing
- Enable verbose logging for debugging

**Production Phase:**
- Enable `continue_on_error` for robustness
- Implement automated retry mechanisms
- Generate comprehensive error reports

## Troubleshooting

### Configuration Issues

**Invalid Paths:**
```bash
# Check file existence
python complete_pipeline.py --config config.yaml --validate_only
```

**Parameter Validation:**
```python
from workflows import validate_pipeline_inputs, PipelineConfig

config = PipelineConfig.from_file('config.yaml')
try:
    validate_pipeline_inputs(config)
except ValueError as e:
    print(f"Validation error: {e}")
```

### Processing Failures

**Memory Issues:**
- Reduce `num_processes` and `memory_limit_gb`
- Process samples sequentially
- Use smaller scale ranges

**File Access Errors:**
- Check file permissions
- Verify path accessibility
- Ensure sufficient disk space

**Resource Contention:**
- Reduce parallel processing load
- Stagger batch processing jobs
- Monitor system resource usage

### Batch Processing Issues

**Failed Samples:**
```bash
# Review failed configurations
ls batch_output/failed_configs/

# Check detailed logs
cat batch_output/batch_logs/batch_processing_*.log
```

**Incomplete Processing:**
```python
# Resume from batch summary
with open('batch_summary.json') as f:
    summary = json.load(f)
    
failed_samples = [r for r in summary['processing_results'] if not r['success']]
print(f"Failed samples: {[r['sample_name'] for r in failed_samples]}")
```

## Best Practices

### Configuration Management

1. **Use version control** for configuration files
2. **Template-based approach** for similar samples
3. **Validation before processing** to catch errors early
4. **Document parameter choices** in configuration comments

### Batch Processing

1. **Start with small batches** to validate settings
2. **Monitor resource usage** during processing
3. **Implement checkpointing** for long-running jobs
4. **Preserve failed configurations** for debugging

### Resource Planning

1. **Estimate processing time** based on data size
2. **Plan storage requirements** for outputs
3. **Consider network bandwidth** for remote data
4. **Implement monitoring** for long-running processes

## Examples

See the `examples/` directory for:
- Complete configuration templates
- Batch processing workflows  
- Integration examples
- Performance optimization guides