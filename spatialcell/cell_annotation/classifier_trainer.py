#!/usr/bin/env python3
"""
TopAct Classifier Training Module for Spatialcell Pipeline

This module trains TopAct classifiers for spatial cell type annotation from single-cell 
reference data, supporting time-point-specific classifier generation.

Built upon the TopAct framework: https://gitlab.com/kfbenjamin/topact.git
Training uses Support Vector Classifier (SVC) from TopAct's implementation.

Author: Xinyan
License: MIT
"""

import os
import argparse
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csr_matrix
from topact.countdata import CountMatrix  # Using TopAct framework
from topact.classifier import SVCClassifier, train_from_countmatrix  # Using TopAct SVC implementation
from joblib import dump
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np


def setup_logging(sample_name: str, output_dir: str) -> None:
    """
    Setup logging configuration for classifier training.
    
    Args:
        sample_name (str): Sample/time-point name for log file naming
        output_dir (str): Directory for log file output
    """
    log_file = os.path.join(output_dir, f"{sample_name}_TopAct_Classifier_Training.log")
    
    # Clear any existing handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='a'),
            logging.StreamHandler(sys.stdout)
        ],
        force=True
    )
    
    logging.info(f"Logging initialized for sample: {sample_name}")


def validate_file_path(file_path: str) -> str:
    """
    Validate that a file exists and return the path.
    
    Args:
        file_path (str): Path to validate
        
    Returns:
        str: Validated file path
        
    Raises:
        FileNotFoundError: If file doesn't exist
    """
    if not Path(file_path).exists():
        raise FileNotFoundError(f"Required file not found: {file_path}")
    return file_path


def load_reference_data(input_dir: str, label_column: str = "celltype_merge") -> sc.AnnData:
    """
    Load reference single-cell data from R-exported matrix market format.
    
    This function loads data exported from R scripts containing:
    - Sparse count matrix (MTX format)
    - Gene names (CSV)
    - Cell barcodes (CSV) 
    - Cell metadata (CSV)
    
    Args:
        input_dir (str): Directory containing exported R data files
        label_column (str): Column name for cell type labels
        
    Returns:
        sc.AnnData: Loaded reference data as AnnData object
        
    Raises:
        FileNotFoundError: If required files are missing
        ValueError: If data dimensions don't match
    """
    # Construct expected file paths
    mtx_file = validate_file_path(os.path.join(input_dir, "all_samples_counts.mtx"))
    genes_file = validate_file_path(os.path.join(input_dir, "all_samples_genes.csv"))
    barcodes_file = validate_file_path(os.path.join(input_dir, "all_samples_barcodes.csv"))
    meta_file = validate_file_path(os.path.join(input_dir, "all_samples_meta_data.csv"))

    # Load sparse count matrix
    logging.info(f"Loading sparse count matrix: {mtx_file}")
    counts = mmread(mtx_file)
    # Transpose matrix from (genes, cells) to (cells, genes) for AnnData
    X = csr_matrix(counts.T)
    
    # Load gene names
    logging.info(f"Loading gene names: {genes_file}")
    genes_df = pd.read_csv(genes_file)
    genes = genes_df["Gene"].tolist()
    
    # Load cell barcodes
    logging.info(f"Loading cell barcodes: {barcodes_file}")
    barcodes_df = pd.read_csv(barcodes_file)
    cells = barcodes_df["Barcode"].tolist()
    
    # Validate matrix dimensions
    if X.shape != (len(cells), len(genes)):
        raise ValueError(f"Matrix dimensions {X.shape} don't match cells ({len(cells)}) and genes ({len(genes)})")
    
    # Load cell metadata
    logging.info(f"Loading cell metadata: {meta_file}")
    meta_df = pd.read_csv(meta_file, index_col=0)
    
    # Ensure barcode consistency
    if not all(cell in meta_df.index for cell in cells):
        raise ValueError("Cell barcodes in matrix don't match metadata")
    meta_df = meta_df.loc[cells]  # Ensure consistent ordering
    
    # Validate label column
    if label_column not in meta_df.columns:
        raise ValueError(f"Label column '{label_column}' not found in metadata")
    
    if meta_df[label_column].isna().any():
        na_count = meta_df[label_column].isna().sum()
        logging.warning(f"Found {na_count} missing values in label column '{label_column}'")
    
    # Create AnnData object
    adata = sc.AnnData(X=X, obs=meta_df, var=pd.DataFrame(index=genes))
    
    logging.info(f"Reference data loaded successfully:")
    logging.info(f"  Cells: {adata.n_obs}")
    logging.info(f"  Genes: {adata.n_vars}")
    logging.info(f"  Cell types: {len(adata.obs[label_column].unique())}")
    logging.info(f"Metadata preview:\n{adata.obs.head()}")
    
    return adata


def get_time_point_mapping() -> Dict[str, List[str]]:
    """
    Define mapping from time points to sample identifiers.
    
    Returns:
        dict: Mapping of time points to sample IDs
    """
    return {
        "E14.5": ["E14-CS1", "E14-CS2", "E14-WT1", "E14-WT2"],
        "E18.5": ["E18-CS1", "E18-CS2", "E18-WT1", "E18-WT2"], 
        "P3": ["P3-CS1", "P3-CS2", "P3-CS3", "P3-WT1", "P3-WT2", "P3-WT3"]
    }


def subset_by_time_point(adata: sc.AnnData, time_point: str) -> sc.AnnData:
    """
    Extract subset of reference data for specific time point.
    
    Args:
        adata (sc.AnnData): Complete reference dataset
        time_point (str): Time point identifier
        
    Returns:
        sc.AnnData: Subset data for specified time point
        
    Raises:
        ValueError: If time point is not supported
    """
    time_point_mapping = get_time_point_mapping()
    
    if time_point not in time_point_mapping:
        available_points = list(time_point_mapping.keys())
        raise ValueError(f"Time point '{time_point}' not supported. Available: {available_points}")
    
    # Extract cells for specified time point
    sample_ids = time_point_mapping[time_point]
    subset_mask = adata.obs['orig.ident'].isin(sample_ids)
    subset_adata = adata[subset_mask].copy()
    
    logging.info(f"Time point {time_point} subset:")
    logging.info(f"  Samples included: {sample_ids}")
    logging.info(f"  Cells extracted: {subset_adata.n_obs}")
    
    return subset_adata


def train_topact_classifier(adata: sc.AnnData, time_point: str, 
                           label_column: str = "celltype_merge") -> SVCClassifier:
    """
    Train TopAct SVC classifier for specified time point.
    
    This function uses TopAct's Support Vector Classifier implementation
    for spatial cell type annotation.
    
    Args:
        adata (sc.AnnData): Reference data for training
        time_point (str): Time point identifier
        label_column (str): Column containing cell type labels
        
    Returns:
        SVCClassifier: Trained TopAct classifier
        
    Raises:
        ValueError: If insufficient data for training
    """
    logging.info(f"Training TopAct classifier for {time_point}...")
    
    if adata.n_obs == 0:
        raise ValueError(f"No cells available for time point {time_point}")
    
    # Extract cell type labels
    labels = adata.obs[label_column]
    
    if labels.empty or labels.isna().all():
        raise ValueError(f"No valid labels found in column '{label_column}'")
    
    # Check for sufficient diversity
    unique_labels = labels.dropna().unique()
    if len(unique_labels) < 2:
        raise ValueError(f"Need at least 2 cell types for training, found {len(unique_labels)}")
    
    logging.info(f"Training data summary:")
    logging.info(f"  Total cells: {adata.n_obs}")
    logging.info(f"  Cell types: {len(unique_labels)}")
    logging.info(f"  Cell type distribution:\n{labels.value_counts()}")
    
    # Convert to TopAct CountMatrix format
    scdata = CountMatrix(
        adata.X, 
        genes=adata.var_names.tolist(), 
        samples=adata.obs.index.tolist()
    )
    scdata.add_metadata("cluster", labels.values)
    
    # Initialize and train TopAct SVC classifier
    clf = SVCClassifier()
    train_from_countmatrix(clf, scdata, "cluster")
    
    # Store training genes for consistency
    clf.train_genes = adata.var_names.tolist()
    
    logging.info(f"TopAct classifier training completed:")
    logging.info(f"  Model classes: {clf.classes}")
    logging.info(f"  Training genes: {len(clf.train_genes)}")
    
    return clf


def save_classifier(classifier: SVCClassifier, output_path: str) -> None:
    """
    Save trained classifier to disk.
    
    Args:
        classifier (SVCClassifier): Trained TopAct classifier
        output_path (str): Output file path
    """
    # Ensure output directory exists
    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)
    
    # Save classifier
    dump(classifier, output_path)
    logging.info(f"Classifier saved to: {output_path}")


def train_multiple_time_points(time_points: List[str], input_dir: str, 
                              output_dir: str, label_column: str = "celltype_merge") -> Dict[str, str]:
    """
    Train classifiers for multiple time points.
    
    Args:
        time_points (list): List of time points to process
        input_dir (str): Input directory with reference data
        output_dir (str): Output directory for classifiers
        label_column (str): Label column name
        
    Returns:
        dict: Mapping of time points to classifier file paths
    """
    # Load complete reference dataset once
    logging.info("Loading complete reference dataset...")
    adata = load_reference_data(input_dir, label_column)
    
    classifier_paths = {}
    
    for time_point in time_points:
        try:
            # Setup logging for this time point
            setup_logging(time_point, output_dir)
            
            # Extract subset for time point
            subset_adata = subset_by_time_point(adata, time_point)
            
            # Train classifier
            classifier = train_topact_classifier(subset_adata, time_point, label_column)
            
            # Save classifier
            output_path = os.path.join(output_dir, f"clf_{time_point}.joblib")
            save_classifier(classifier, output_path)
            
            classifier_paths[time_point] = output_path
            
            logging.info(f"Successfully trained classifier for {time_point}")
            
        except Exception as e:
            logging.error(f"Failed to train classifier for {time_point}: {e}")
            continue
    
    return classifier_paths


def main():
    """Main entry point for TopAct classifier training."""
    parser = argparse.ArgumentParser(
        description="Train TopAct classifiers for spatial cell type annotation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Train classifiers for specific time points
  python classifier_trainer.py --time_points E14.5 E18.5 P3 \\
    --input_dir /data/reference/ --output_dir /data/classifiers/

  # Train with custom label column
  python classifier_trainer.py --time_points P3 \\
    --input_dir /data/reference/ --output_dir /data/classifiers/ \\
    --label_column custom_celltype

Input directory should contain:
  - all_samples_counts.mtx (sparse count matrix)
  - all_samples_genes.csv (gene names)
  - all_samples_barcodes.csv (cell barcodes)
  - all_samples_meta_data.csv (cell metadata with labels)
        """
    )
    
    parser.add_argument("--time_points", type=str, nargs="+", required=True,
                        help="Time points to train classifiers for (e.g., E14.5 E18.5 P3)")
    parser.add_argument("--input_dir", type=str, required=True,
                        help="Directory containing R-exported reference data files")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory to save trained classifiers and logs")
    parser.add_argument("--label_column", type=str, default="celltype_merge",
                        help="Metadata column containing cell type labels (default: celltype_merge)")
    
    args = parser.parse_args()
    
    try:
        # Validate input directory
        if not Path(args.input_dir).exists():
            raise FileNotFoundError(f"Input directory not found: {args.input_dir}")
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Train classifiers
        classifier_paths = train_multiple_time_points(
            args.time_points, args.input_dir, args.output_dir, args.label_column
        )
        
        # Summary
        if classifier_paths:
            print(f"\nClassifier training completed successfully!")
            print(f"Trained classifiers for {len(classifier_paths)} time points:")
            for time_point, path in classifier_paths.items():
                print(f"  {time_point}: {path}")
        else:
            print("No classifiers were successfully trained.")
            return 1
            
    except Exception as e:
        logging.error(f"Training failed: {e}")
        print(f"Classifier training failed: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())# Version: 1.0.7
