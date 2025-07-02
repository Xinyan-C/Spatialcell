#!/usr/bin/env python3
"""
TopAct Spatial Annotation Processor for Spatialcell Pipeline

This module processes high-definition spatial transcriptomics data (2μm bin) for cell type annotation
using TopAct framework with cell boundary constraints from bin2cell segmentation.

Built upon:
- TopAct framework: https://gitlab.com/kfbenjamin/topact.git
- bin2cell package: https://github.com/Teichlab/bin2cell.git

Author: Xinyan
License: MIT
"""

import os
import re
import resource
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import logging
from topact.countdata import CountMatrix  # Using TopAct framework
from topact.classifier import SVCClassifier  # Using TopAct classifiers
import bin2cell as b2c  # Using bin2cell for spatial processing
from joblib import load, dump
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import warnings

# Import enhanced spatial analysis from utils module
try:
    from ..utils.spatial_analysis import EnhancedCountGrid as spatial_CountGrid
    SPATIAL_AVAILABLE = True
except ImportError:
    SPATIAL_AVAILABLE = False
    warnings.warn("Enhanced spatial analysis module not found. Some functionality may be limited.")


def setup_logging(sample_name: str, log_filename: Optional[str] = None) -> None:
    """
    Configure logging for spatial annotation processing.
    
    Args:
        sample_name (str): Sample identifier for log naming
        log_filename (str, optional): Custom log filename
    """
    if log_filename is None:
        log_filename = f"{sample_name.replace('.', '_')}_spatial_annotation.log"
    
    # Clear existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_filename, mode="a", encoding="utf-8")
        ],
        force=True
    )
    
    logging.info(f"Spatial annotation processing started for sample: {sample_name}")


def setup_computational_resources(threads: str = "80", mem_gb: int = 600) -> None:
    """
    Configure computational resources for processing.
    
    Args:
        threads (str): Number of threads to use
        mem_gb (int): Memory limit in GB
    """
    # Set threading environment variables
    os.environ["OMP_NUM_THREADS"] = threads
    os.environ["MKL_NUM_THREADS"] = threads  
    os.environ["NUMEXPR_NUM_THREADS"] = threads
    logging.info(f"Thread limit set to: {threads}")
    
    # Set memory limit
    mem_limit = mem_gb * 1024**3
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        resource.setrlimit(resource.RLIMIT_AS, (mem_limit, hard))
        logging.info(f"Memory soft limit set to {mem_gb}GB")
    except Exception as e:
        logging.warning(f"Failed to set memory limit: {e}")


def parse_roi_coordinates(roi_file: str) -> Dict[str, Dict[str, float]]:
    """
    Parse ROI coordinates from coordinate file.
    
    Args:
        roi_file (str): Path to ROI coordinate file
        
    Returns:
        dict: ROI coordinates mapping
        
    Raises:
        FileNotFoundError: If ROI file doesn't exist
        ValueError: If file format is invalid
    """
    if not Path(roi_file).exists():
        raise FileNotFoundError(f"ROI file not found: {roi_file}")
    
    logging.info(f"Parsing ROI coordinates from: {roi_file}")
    rois = {}
    
    try:
        with open(roi_file, 'r', encoding='utf-8') as f:
            content = f.read()
            
        # Parse coordinate patterns (handle multiple formats)
        pattern = r'([A-Za-z]+\d*(?:\.\d+)?\s*-\s*\w+\d*)\s*(?:Rectangle coordinates|矩形框的坐标范围):\s*X:\s*([\d\.\-]+)\s*-\s*([\d\.\-]+)\s*Y:\s*([\d\.\-]+)\s*-\s*([\d\.\-]+)'
        matches = re.findall(pattern, content, re.MULTILINE | re.DOTALL)
        
        for match in matches:
            roi_name = match[0].strip()
            x_min, x_max, y_min, y_max = map(float, match[1:5])
            
            rois[roi_name] = {
                'x_min': x_min, 'x_max': x_max,
                'y_min': y_min, 'y_max': y_max
            }
            
    except Exception as e:
        raise ValueError(f"Failed to parse ROI file: {e}")
    
    logging.info(f"Parsed {len(rois)} ROI regions")
    return rois


def load_hd_spatial_data(expression_path: str, positions_path: str, 
                        train_genes: List[str]) -> Tuple[sc.AnnData, pd.DataFrame]:
    """
    Load high-definition spatial transcriptomics data.
    
    Args:
        expression_path (str): Path to expression matrix
        positions_path (str): Path to spatial positions
        train_genes (list): List of genes used in training
        
    Returns:
        tuple: (AnnData object, positions DataFrame)
    """
    logging.info(f"Loading HD expression data: {expression_path}")
    
    if not Path(expression_path).exists():
        raise FileNotFoundError(f"Expression data not found: {expression_path}")
    
    # Load expression data
    adata_hd = sc.read_10x_mtx(expression_path, var_names='gene_symbols', cache=True)
    
    logging.info(f"Loading spatial positions: {positions_path}")
    if not Path(positions_path).exists():
        raise FileNotFoundError(f"Positions file not found: {positions_path}")
        
    # Load spatial coordinates
    positions = pd.read_parquet(positions_path)
    positions.rename(columns={
        'pxl_col_in_fullres': 'x', 
        'pxl_row_in_fullres': 'y'
    }, inplace=True)
    
    logging.info(f"Loaded {adata_hd.n_obs} spots and {adata_hd.n_vars} genes")
    logging.info(f"Loaded {len(positions)} spatial positions")
    
    return adata_hd, positions


def process_hd_data_with_cell_labels(expression_path: str, positions_path: str, 
                                   train_genes: List[str], rois: Dict[str, Dict[str, float]],
                                   bin2cell_adata: sc.AnnData, 
                                   labels_column: str = 'labels_qupath_expanded') -> Dict[str, pd.DataFrame]:
    """
    Process HD data and assign cell labels from bin2cell segmentation.
    
    This function integrates bin2cell segmentation results with HD spatial data,
    enabling cell-boundary-constrained analysis.
    
    Args:
        expression_path (str): Path to HD expression matrix
        positions_path (str): Path to HD spatial positions  
        train_genes (list): Genes used for training classifiers
        rois (dict): ROI coordinate definitions
        bin2cell_adata (sc.AnnData): bin2cell processed data with cell labels
        labels_column (str): Column containing cell labels
        
    Returns:
        dict: ROI-specific spatial data with cell labels
    """
    # Load HD data
    adata_hd, positions = load_hd_spatial_data(expression_path, positions_path, train_genes)
    
    # Extract spatial coordinates and cell labels from bin2cell data
    bin2cell_spatial = bin2cell_adata.obsm['spatial']
    bin2cell_labels = bin2cell_adata.obs[labels_column]
    
    # Create lookup dictionary for cell labels (using float coordinates)
    label_dict = {
        (float(x), float(y)): label 
        for x, y, label in zip(bin2cell_spatial[:, 0], bin2cell_spatial[:, 1], bin2cell_labels)
    }
    
    logging.info(f"Created cell label lookup with {len(label_dict)} entries")
    
    roi_spatial_data = {}
    
    for roi_name, coords in rois.items():
        logging.info(f"Processing ROI: {roi_name}")
        
        # Filter positions within ROI
        roi_mask = (
            (positions['x'] >= coords['x_min']) & (positions['x'] <= coords['x_max']) &
            (positions['y'] >= coords['y_min']) & (positions['y'] <= coords['y_max'])
        )
        positions_roi = positions[roi_mask].copy()
        roi_barcodes = positions_roi['barcode']
        
        # Filter expression data for ROI
        adata_hd_roi = adata_hd[adata_hd.obs_names.isin(roi_barcodes)].copy()
        
        # Align gene list with training genes
        df_expr = pd.DataFrame.sparse.from_spmatrix(
            adata_hd_roi.X, 
            index=adata_hd_roi.obs_names, 
            columns=adata_hd_roi.var_names
        )
        df_expr_aligned = df_expr.reindex(columns=train_genes, fill_value=0).astype(pd.SparseDtype("int64", 0))
        df_expr_aligned.index.name = "barcode"
        
        # Convert to coordinate-oriented format
        coo = df_expr_aligned.sparse.to_coo()
        expr_long = pd.DataFrame({
            "barcode": df_expr_aligned.index[coo.row],
            "gene": df_expr_aligned.columns[coo.col], 
            "counts": coo.data
        })
        
        # Merge with spatial coordinates
        spatial_data = pd.merge(positions_roi, expr_long, on="barcode", how="inner")
        
        # Assign cell labels using bin2cell segmentation
        spatial_data['cell_label'] = spatial_data.apply(
            lambda row: label_dict.get((float(row['x']), float(row['y'])), 0), axis=1
        )
        
        roi_spatial_data[roi_name] = spatial_data
        logging.info(f"ROI {roi_name}: {len(spatial_data)} data points with cell labels")
    
    return roi_spatial_data


def build_spatial_countgrid(spatial_data: pd.DataFrame, genes: List[str]) -> Any:
    """
    Build CountGrid object for spatial analysis using TopAct framework.
    
    Args:
        spatial_data (pd.DataFrame): Spatial expression data
        genes (list): List of genes
        
    Returns:
        CountGrid: Spatial count grid object
    """
    if not SPATIAL_AVAILABLE:
        raise ImportError("Xinyantopact module required for CountGrid construction")
    
    logging.info("Building CountGrid object for spatial analysis...")
    
    try:
        # Use custom spatial module to build CountGrid
        sd = spatial_CountGrid.from_coordinate_table(
            spatial_data, 
            genes=genes, 
            count_col="counts", 
            gene_col="gene"
        )
        logging.info("CountGrid construction completed")
        return sd
        
    except Exception as e:
        raise RuntimeError(f"Failed to build CountGrid: {e}")


def create_cell_constrained_neighborhood(sd: Any, center: Tuple[float, float], scale: float) -> pd.DataFrame:
    """
    Custom neighborhood function constrained to cell boundaries.
    
    This function implements cell-boundary-aware neighborhood selection,
    ensuring that spatial analysis respects individual cell segmentation.
    
    Args:
        sd: CountGrid object
        center (tuple): Center coordinates (x, y)
        scale (float): Neighborhood scale
        
    Returns:
        pd.DataFrame: Filtered neighborhood data within same cell
    """
    center_x, center_y = center
    
    # Get cell label for center point
    center_matches = sd.table[(sd.table['x'] == center_x) & (sd.table['y'] == center_y)]
    
    if center_matches.empty:
        return pd.DataFrame()
        
    center_label = center_matches['cell_label'].iloc[0]
    
    # Return empty if no valid cell assignment
    if center_label == 0:
        return pd.DataFrame()
    
    # Filter to same cell and within distance
    same_cell_data = sd.table[sd.table['cell_label'] == center_label]
    distances = np.sqrt((same_cell_data['x'] - center_x)**2 + (same_cell_data['y'] - center_y)**2)
    within_scale = same_cell_data[distances <= scale]
    
    return within_scale


def perform_multiscale_classification(sd: Any, classifier: SVCClassifier, outfile: str,
                                    min_scale: float = 3.0, max_scale: float = 9.0,
                                    num_proc: int = 100, mpp: Optional[float] = None) -> np.ndarray:
    """
    Perform multi-scale spatial classification with cell boundary constraints.
    
    Uses TopAct framework for classification while respecting cell boundaries
    defined by bin2cell segmentation.
    
    Args:
        sd: CountGrid object
        classifier (SVCClassifier): Trained TopAct classifier
        outfile (str): Output file path for results
        min_scale (float): Minimum neighborhood scale in micrometers
        max_scale (float): Maximum neighborhood scale in micrometers
        num_proc (int): Number of parallel processes
        mpp (float, optional): Microns per pixel for scale conversion
        
    Returns:
        np.ndarray: Classification confidence matrix
    """
    logging.info(f"Starting multi-scale classification with cell boundary constraints")
    logging.info(f"Scale range: {min_scale} - {max_scale} micrometers")
    logging.info(f"Using {num_proc} processes")
    
    try:
        # Perform classification using TopAct with custom neighborhood function
        confidence_matrix = sd.classify_parallel(
            classifier, 
            min_scale=min_scale,
            max_scale=max_scale, 
            outfile=outfile,
            mpp=mpp,
            num_proc=num_proc,
            verbose=False,
            neighborhood_func=create_cell_constrained_neighborhood
        )
        
        logging.info(f"Classification completed, results saved to: {outfile}")
        return confidence_matrix
        
    except Exception as e:
        raise RuntimeError(f"Classification failed: {e}")


def get_sample_mpp(visium_path: str, source_image_path: str) -> float:
    """
    Extract microns per pixel (MPP) from Visium data.
    
    Args:
        visium_path (str): Path to Visium data directory
        source_image_path (str): Path to source image
        
    Returns:
        float: Microns per pixel value
    """
    try:
        # Load Visium data using bin2cell
        adata = b2c.read_visium(visium_path, source_image_path=source_image_path)
        library = list(adata.uns['spatial'].keys())[0]
        mpp = adata.uns['spatial'][library]['scalefactors']['microns_per_pixel']
        
        logging.info(f"Extracted MPP: {mpp} micrometers/pixel")
        return mpp
        
    except Exception as e:
        raise ValueError(f"Failed to extract MPP: {e}")


def process_sample_annotation(args: argparse.Namespace) -> None:
    """
    Main processing function for spatial cell type annotation.
    
    Args:
        args: Command line arguments containing processing parameters
    """
    setup_logging(args.sample)
    setup_computational_resources(args.threads, args.mem_gb)
    
    # Create output directory
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Load MPP from Visium data
    mpp = get_sample_mpp(os.path.dirname(args.expr_path), args.source_image_path)
    
    # Load bin2cell data for cell labels
    if not Path(args.bin2cell_dir).exists():
        raise FileNotFoundError(f"bin2cell data not found: {args.bin2cell_dir}")
    
    bin2cell_adata = sc.read_h5ad(args.bin2cell_dir)
    logging.info(f"Loaded bin2cell data with {bin2cell_adata.n_obs} cells")
    
    # Parse ROI coordinates
    rois = parse_roi_coordinates(args.roi_file)
    
    # Load classifier and get training genes
    if not Path(args.clf_path).exists():
        raise FileNotFoundError(f"Classifier not found: {args.clf_path}")
    
    clf_example = load(args.clf_path)
    train_genes = clf_example.train_genes
    logging.info(f"Loaded classifier with {len(train_genes)} training genes")
    
    # Process HD data and assign cell labels
    spatial_data_path = os.path.join(args.out_dir, f"spatial_data_{args.sample}_roi.joblib")
    
    if Path(spatial_data_path).exists():
        roi_spatial_data = load(spatial_data_path)
        logging.info(f"Loaded cached spatial data for {args.sample}")
    else:
        roi_spatial_data = process_hd_data_with_cell_labels(
            args.expr_path, args.pos_path, train_genes, rois, 
            bin2cell_adata, args.labels
        )
        dump(roi_spatial_data, spatial_data_path)
        logging.info(f"Cached spatial data for {args.sample}")
    
    # Process each ROI independently
    for roi_name in rois.keys():
        logging.info(f"Processing ROI: {roi_name}")
        
        spatial_data = roi_spatial_data[roi_name]
        safe_roi_name = roi_name.replace(' ', '_')
        
        # Define file paths
        sd_path = os.path.join(args.out_dir, f"sd_{args.sample}_{safe_roi_name}.joblib")
        outfile = os.path.join(args.out_dir, f"outfile_{args.sample}_{safe_roi_name}.npy")
        
        # Build CountGrid
        if Path(sd_path).exists():
            sd = load(sd_path)
            logging.info(f"Loaded cached CountGrid for {roi_name}")
        else:
            sd = build_spatial_countgrid(spatial_data, train_genes)
            dump(sd, sd_path)
            logging.info(f"Cached CountGrid for {roi_name}")
        
        # Load classifier for this ROI
        if not Path(args.clf_path).exists():
            raise FileNotFoundError(f"Classifier not found: {args.clf_path}")
        
        classifier = load(args.clf_path)
        logging.info(f"Using classifier: {args.clf_path}")
        
        # Perform classification if not already done
        if not Path(outfile).exists():
            confidence_matrix = perform_multiscale_classification(
                sd, classifier, outfile,
                min_scale=args.min_scale,
                max_scale=args.max_scale, 
                num_proc=args.num_proc,
                mpp=mpp
            )
            logging.info(f"Classification completed for {roi_name}")
        else:
            logging.info(f"Classification results already exist for {roi_name}")
    
    logging.info(f"Sample {args.sample} processing completed successfully!")


def main():
    """Main entry point for spatial annotation processing."""
    parser = argparse.ArgumentParser(
        description="Process spatial transcriptomics data for cell type annotation using TopAct",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python annotation_processor.py --sample E14.5 \\
    --out_dir /data/results/ \\
    --expr_path /data/visium/matrix.h5 \\
    --pos_path /data/positions.parquet \\
    --source_image_path /data/tissue.tif \\
    --roi_file /data/regions.txt \\
    --clf_path /data/classifiers/clf_E14.5.joblib \\
    --bin2cell_dir /data/bin2cell_results.h5ad

This tool integrates TopAct classification with bin2cell segmentation
for accurate spatial cell type annotation.
        """
    )
    
    # Required parameters
    parser.add_argument("--sample", type=str, required=True,
                        help="Sample identifier")
    parser.add_argument("--out_dir", type=str, required=True,
                        help="Output directory for results")
    parser.add_argument("--expr_path", type=str, required=True,
                        help="Path to HD expression matrix")
    parser.add_argument("--pos_path", type=str, required=True,
                        help="Path to HD spatial positions")
    parser.add_argument("--source_image_path", type=str, required=True,
                        help="Path to source tissue image")
    parser.add_argument("--roi_file", type=str, required=True,
                        help="Path to ROI coordinates file")
    parser.add_argument("--clf_path", type=str, required=True,
                        help="Path to trained classifier")
    parser.add_argument("--bin2cell_dir", type=str, required=True,
                        help="Path to bin2cell processed data")
    
    # Optional parameters
    parser.add_argument("--min_scale", type=float, default=3.0,
                        help="Minimum neighborhood radius in micrometers (default: 3.0)")
    parser.add_argument("--max_scale", type=float, default=9.0,
                        help="Maximum neighborhood radius in micrometers (default: 9.0)")
    parser.add_argument("--num_proc", type=int, default=80,
                        help="Number of parallel processes (default: 80)")
    parser.add_argument("--threads", type=str, default="80",
                        help="Number of threads (default: 80)")
    parser.add_argument("--mem_gb", type=int, default=600,
                        help="Memory limit in GB (default: 600)")
    parser.add_argument("--labels", type=str, default="labels_qupath_expanded",
                        help="Cell label column from bin2cell (default: labels_qupath_expanded)")
    
    args = parser.parse_args()
    
    # Construct classifier path if not absolute
    if not os.path.isabs(args.clf_path):
        clf_filename = f"clf_{args.sample}.joblib"
        args.clf_path = os.path.join(os.path.dirname(args.clf_path), clf_filename)
    
    try:
        process_sample_annotation(args)
        print(f"Spatial annotation processing completed for {args.sample}!")
        
    except Exception as e:
        logging.error(f"Processing failed: {e}")
        print(f"Processing failed: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())# Version: 1.0.7
