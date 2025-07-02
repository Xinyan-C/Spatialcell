#!/usr/bin/env python3
"""
Spatial Processing Module for Spatialcell Pipeline

This module integrates bin2cell functionality for comprehensive spatial transcriptomics analysis,
including cell segmentation, label expansion, and spatial visualization.

Built upon the bin2cell package: https://github.com/BayraktarLab/bin2cell

Author: Xinyan
License: MIT
"""

import os
import argparse
import logging
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scanpy as sc
import colorcet as cc
import matplotlib.colors as mcolors
import bin2cell as b2c  
from scipy.sparse import load_npz
import numpy as np
from pathlib import Path


def setup_logging(output_dir):
    """
    Configure logging to output to both console and a file.
    
    Args:
        output_dir (str): Directory to save log file
    """
    log_file = os.path.join(output_dir, "spatial_processing.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )


def read_regions(region_file):
    """
    Read and parse region coordinates from a file.
    
    Args:
        region_file (str): Path to region coordinates file
        
    Returns:
        dict: Dictionary containing region coordinates
    """
    logging.info("Reading region file...")
    
    if not Path(region_file).exists():
        raise FileNotFoundError(f"Region file not found: {region_file}")
    
    with open(region_file, 'r', encoding='utf-8') as file:
        doc_content = file.read()
    
    regions = {}
    lines = doc_content.strip().split('\n')
    
    for i in range(0, len(lines), 4):
        if i + 2 < len(lines):
            # Parse region name - handle both Chinese and English formats
            if "矩形框的坐标范围:" in lines[i]:
                region_name = lines[i].split('矩形框的坐标范围:')[0].strip()
            else:
                region_name = lines[i].split('Rectangle coordinates:')[0].strip()
            
            if "整张图像" in region_name or "whole image" in region_name.lower():
                continue
                
            x_line = lines[i + 1].strip()
            y_line = lines[i + 2].strip()
            
            # Parse coordinates
            x_min, x_max = map(float, x_line.split(':')[1].strip().split(' - '))
            y_min, y_max = map(float, y_line.split(':')[1].strip().split(' - '))
            
            regions[region_name] = {
                "col_min": x_min,
                "col_max": x_max,
                "row_min": y_min,
                "row_max": y_max
            }
    
    logging.info(f"Found {len(regions)} regions.")
    return regions


def destripe_and_save_regions(adata, region_file, save_dir, mpp):
    """
    Visualize and save destripe results for each region.
    
    Args:
        adata: AnnData object
        region_file (str): Path to region coordinates file
        save_dir (str): Directory to save visualization results
        mpp (float): Microns per pixel
    """
    logging.info(f"Visualizing destripe results in {save_dir}...")
    regions = read_regions(region_file)
    os.makedirs(save_dir, exist_ok=True)
    basename = os.path.basename(save_dir.rstrip(os.sep))

    for region_name, region_coords in regions.items():
        mask = (
            (adata.obsm['spatial'][:, 0] >= region_coords["col_min"]) &
            (adata.obsm['spatial'][:, 0] <= region_coords["col_max"]) &
            (adata.obsm['spatial'][:, 1] >= region_coords["row_min"]) &
            (adata.obsm['spatial'][:, 1] <= region_coords["row_max"])
        )
        bdata = adata[mask]
        logging.info(f"{region_name} - Number of spots selected: {bdata.n_obs}")

        titles = [
            f"{region_name} - Image",
            f"{region_name} - Raw Counts",
            f"{region_name} - Adjusted Counts"
        ]
        
        sc.pl.spatial(
            bdata,
            color=[None, "n_counts", "n_counts_adjusted"],
            color_map="OrRd",
            img_key=f"{mpp}_mpp_150_buffer",
            basis="spatial_cropped_150_buffer",
            title=titles,
            show=False
        )
        
        save_path = os.path.join(save_dir, f"{region_name}_{basename}.pdf")
        plt.tight_layout()
        plt.savefig(save_path, dpi=1200, bbox_inches="tight", pad_inches=0.5)
        plt.close()
        logging.info(f"Saved: {save_path}")


def visualize_labels_in_regions(npz_path, region_file, save_dir):
    """
    Visualize labels from NPZ file for each region.
    
    Args:
        npz_path (str): Path to NPZ labels file
        region_file (str): Path to region coordinates file
        save_dir (str): Directory to save visualization results
    """
    logging.info(f"Visualizing NPZ labels in {save_dir}...")
    
    if not Path(npz_path).exists():
        raise FileNotFoundError(f"NPZ file not found: {npz_path}")
    
    labels_sparse = load_npz(npz_path)
    labels_matrix = labels_sparse.toarray()
    regions = read_regions(region_file)
    os.makedirs(save_dir, exist_ok=True)
    basename = os.path.basename(save_dir.rstrip(os.sep))

    for region_name, region_coords in regions.items():
        row_min, row_max = int(region_coords["row_min"]), int(region_coords["row_max"])
        col_min, col_max = int(region_coords["col_min"]), int(region_coords["col_max"])
        cropped_labels = labels_matrix[row_min:row_max, col_min:col_max]
        
        logging.info(f"{region_name} - Cropped labels matrix shape: {cropped_labels.shape}")

        masked_labels = np.ma.masked_where(cropped_labels == 1, cropped_labels)
        unique_labels = np.unique(cropped_labels[cropped_labels != 1])
        colors = [cc.glasbey[i % len(cc.glasbey)] for i in range(len(unique_labels))]
        cmap = mcolors.ListedColormap(colors)
        cmap.set_bad('white')

        plt.imshow(masked_labels, cmap=cmap, interpolation='none')
        plt.colorbar()
        plt.title(f"Labels in {region_name}")
        save_path = os.path.join(save_dir, f"{region_name}_{basename}.pdf")
        plt.tight_layout()
        plt.savefig(save_path, dpi=1200, bbox_inches="tight", pad_inches=0.5)
        plt.close()
        logging.info(f"Saved: {save_path}")


def visualize_after_insert_labels(adata, region_file, save_dir, mpp):
    """
    Visualize labels after insertion for each region.
    
    Args:
        adata: AnnData object
        region_file (str): Path to region coordinates file
        save_dir (str): Directory to save visualization results
        mpp (float): Microns per pixel
    """
    logging.info(f"Visualizing segmentation in {save_dir}...")
    regions = read_regions(region_file)
    os.makedirs(save_dir, exist_ok=True)
    basename = os.path.basename(save_dir.rstrip(os.sep))

    for region_name, region_coords in regions.items():
        mask = (
            (adata.obsm['spatial'][:, 0] >= region_coords["col_min"]) &
            (adata.obsm['spatial'][:, 0] <= region_coords["col_max"]) &
            (adata.obsm['spatial'][:, 1] >= region_coords["row_min"]) &
            (adata.obsm['spatial'][:, 1] <= region_coords["row_max"])
        )
        bdata = adata[mask]
        logging.info(f"{region_name} - Number of spots selected: {bdata.n_obs}")

        bdata = bdata[bdata.obs['labels_qupath'] > 0]
        bdata.obs['labels_qupath'] = bdata.obs['labels_qupath'].astype(int)
        unique_labels = sorted(bdata.obs['labels_qupath'].unique())
        bdata.obs['labels_qupath'] = pd.Categorical(
            bdata.obs['labels_qupath'], categories=unique_labels, ordered=True)
        colors = [cc.glasbey[i % len(cc.glasbey)] for i in range(len(unique_labels))]
        bdata.uns['labels_qupath_colors'] = colors

        titles = [f"{region_name} - Image", f"{region_name} - Labels QuPath"]
        fig = sc.pl.spatial(
            bdata,
            color=[None, "labels_qupath"],
            img_key=f"{mpp}_mpp_150_buffer",
            basis="spatial_cropped_150_buffer",
            title=titles,
            show=False,
            legend_loc=None,
            return_fig=True
        )
        
        axes = fig.axes
        ax = axes[1]
        labels_source = bdata.obs['labels_qupath'].cat.categories
        colors_source = bdata.uns['labels_qupath_colors']
        handles = [Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=c, markersize=10) for c in colors_source]
        ax.legend(handles, labels_source, loc='center left',
                  bbox_to_anchor=(2, 0.5), ncol=15, fontsize=8)

        save_path = os.path.join(save_dir, f"{region_name}_{basename}.pdf")
        plt.tight_layout()
        plt.savefig(save_path, dpi=1200, bbox_inches="tight", pad_inches=0.5)
        plt.close()
        logging.info(f"Saved: {save_path}")


def render_labels_after_insert(adata, region_file, save_dir, npz_path, mpp, source_image_path):
    """
    Render labels after insertion for each region.
    
    Args:
        adata: AnnData object
        region_file (str): Path to region coordinates file
        save_dir (str): Directory to save visualization results
        npz_path (str): Path to NPZ labels file
        mpp (float): Microns per pixel
        source_image_path (str): Path to source image
    """
    logging.info(f"Rendering labels in {save_dir}...")
    regions = read_regions(region_file)
    os.makedirs(save_dir, exist_ok=True)
    basename = os.path.basename(save_dir.rstrip(os.sep))

    for region_name, region_coords in regions.items():
        mask = (
            (adata.obsm['spatial'][:, 0] >= region_coords["col_min"]) &
            (adata.obsm['spatial'][:, 0] <= region_coords["col_max"]) &
            (adata.obsm['spatial'][:, 1] >= region_coords["row_min"]) &
            (adata.obsm['spatial'][:, 1] <= region_coords["row_max"])
        )
        bdata = adata[mask]
        logging.info(f"{region_name} - Number of spots selected: {bdata.n_obs}")

        bdata = bdata[bdata.obs['labels_qupath'] > 0]
        bdata.obs['labels_qupath'] = bdata.obs['labels_qupath'].astype(str)
        crop = b2c.get_crop(bdata, basis="spatial", spatial_key="spatial", mpp=mpp)
        rendered = b2c.view_labels(image_path=source_image_path, labels_npz_path=npz_path, crop=crop)

        plt.figure()
        plt.imshow(rendered)
        plt.axis("off")
        save_path = os.path.join(save_dir, f"{region_name}_{basename}.pdf")
        plt.tight_layout()
        plt.savefig(save_path, dpi=1200, bbox_inches="tight", pad_inches=0.5)
        plt.close()
        logging.info(f"Saved: {save_path}")


def visualize_after_expansion(adata, region_file, save_dir, mpp):
    """
    Visualize expanded labels for each region.
    
    Args:
        adata: AnnData object
        region_file (str): Path to region coordinates file
        save_dir (str): Directory to save visualization results
        mpp (float): Microns per pixel
    """
    logging.info(f"Visualizing expanded labels in {save_dir}...")
    regions = read_regions(region_file)
    os.makedirs(save_dir, exist_ok=True)
    basename = os.path.basename(save_dir.rstrip(os.sep))

    for region_name, region_coords in regions.items():
        mask = (
            (adata.obsm['spatial'][:, 0] >= region_coords["col_min"]) &
            (adata.obsm['spatial'][:, 0] <= region_coords["col_max"]) &
            (adata.obsm['spatial'][:, 1] >= region_coords["row_min"]) &
            (adata.obsm['spatial'][:, 1] <= region_coords["row_max"])
        )
        bdata = adata[mask]
        logging.info(f"{region_name} - Number of spots selected: {bdata.n_obs}")

        bdata = bdata[bdata.obs['labels_qupath_expanded'] > 0]
        bdata.obs['labels_qupath_expanded'] = bdata.obs['labels_qupath_expanded'].astype(int)
        unique_labels = sorted(bdata.obs['labels_qupath_expanded'].unique())
        bdata.obs['labels_qupath_expanded'] = pd.Categorical(
            bdata.obs['labels_qupath_expanded'], categories=unique_labels, ordered=True)
        colors = [cc.glasbey[i % len(cc.glasbey)] for i in range(len(unique_labels))]
        bdata.uns['labels_qupath_expanded_colors'] = colors

        titles = [f"{region_name} - Image", f"{region_name} - Expanded Labels"]
        fig = sc.pl.spatial(
            bdata,
            color=[None, "labels_qupath_expanded"],
            img_key=f"{mpp}_mpp_150_buffer",
            basis="spatial_cropped_150_buffer",
            title=titles,
            show=False,
            legend_loc=None,
            return_fig=True
        )
        
        axes = fig.axes
        ax = axes[1]
        labels_source = bdata.obs['labels_qupath_expanded'].cat.categories
        colors_source = bdata.uns['labels_qupath_expanded_colors']
        handles = [Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=c, markersize=10) for c in colors_source]
        ax.legend(handles, labels_source, loc='center left',
                  bbox_to_anchor=(2, 0.5), ncol=15, fontsize=8)

        save_path = os.path.join(save_dir, f"{region_name}_{basename}.pdf")
        plt.tight_layout()
        plt.savefig(save_path, dpi=1200, bbox_inches="tight", pad_inches=0.5)
        plt.close()
        logging.info(f"Saved: {save_path}")


def visualize_gex_labels(adata, region_file, save_dir, mpp):
    """
    Visualize GEX (gene expression) labels for each region.
    
    Args:
        adata: AnnData object
        region_file (str): Path to region coordinates file
        save_dir (str): Directory to save visualization results
        mpp (float): Microns per pixel
    """
    logging.info(f"Visualizing GEX labels in {save_dir}...")
    regions = read_regions(region_file)
    os.makedirs(save_dir, exist_ok=True)
    basename = os.path.basename(save_dir.rstrip(os.sep))

    for region_name, region_coords in regions.items():
        mask = (
            (adata.obsm['spatial'][:, 0] >= region_coords["col_min"]) &
            (adata.obsm['spatial'][:, 0] <= region_coords["col_max"]) &
            (adata.obsm['spatial'][:, 1] >= region_coords["row_min"]) &
            (adata.obsm['spatial'][:, 1] <= region_coords["row_max"])
        )
        bdata = adata[mask]
        logging.info(f"{region_name} - Number of spots selected: {bdata.n_obs}")

        bdata = bdata[bdata.obs['labels_gex'] > 0]
        bdata.obs['labels_gex'] = bdata.obs['labels_gex'].astype(int)
        unique_labels = sorted(bdata.obs['labels_gex'].unique())
        bdata.obs['labels_gex'] = pd.Categorical(
            bdata.obs['labels_gex'], categories=unique_labels, ordered=True)
        colors = [cc.glasbey[i % len(cc.glasbey)] for i in range(len(unique_labels))]
        bdata.uns['labels_gex_colors'] = colors

        titles = [f"{region_name} - Image", f"{region_name} - GEX Labels"]
        fig = sc.pl.spatial(
            bdata,
            color=[None, "labels_gex"],
            img_key=f"{mpp}_mpp_150_buffer",
            basis="spatial_cropped_150_buffer",
            title=titles,
            show=False,
            legend_loc=None,
            return_fig=True
        )
        
        axes = fig.axes
        ax = axes[1]
        labels_source = bdata.obs['labels_gex'].cat.categories
        colors_source = bdata.uns['labels_gex_colors']
        handles = [Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=c, markersize=10) for c in colors_source]
        ax.legend(handles, labels_source, loc='center left',
                  bbox_to_anchor=(2, 0.5), ncol=15, fontsize=8)

        save_path = os.path.join(save_dir, f"{region_name}_{basename}.pdf")
        plt.tight_layout()
        plt.savefig(save_path, dpi=1200, bbox_inches="tight", pad_inches=0.5)
        plt.close()
        logging.info(f"Saved: {save_path}")


def render_gex_labels(adata, region_file, save_dir, prob_thresh, nms_thresh, mpp):
    """
    Render GEX labels for each region.
    
    Args:
        adata: AnnData object
        region_file (str): Path to region coordinates file
        save_dir (str): Directory to save visualization results
        prob_thresh (float): Probability threshold for StarDist
        nms_thresh (float): NMS threshold for StarDist
        mpp (float): Microns per pixel
    """
    logging.info(f"Rendering GEX labels in {save_dir}...")
    regions = read_regions(region_file)
    os.makedirs(save_dir, exist_ok=True)
    basename = os.path.basename(save_dir.rstrip(os.sep))

    for region_name, region_coords in regions.items():
        mask = (
            (adata.obsm['spatial'][:, 0] >= region_coords["col_min"]) &
            (adata.obsm['spatial'][:, 0] <= region_coords["col_max"]) &
            (adata.obsm['spatial'][:, 1] >= region_coords["row_min"]) &
            (adata.obsm['spatial'][:, 1] <= region_coords["row_max"])
        )
        bdata = adata[mask]
        logging.info(f"{region_name} - Number of spots selected: {bdata.n_obs}")

        bdata = bdata[bdata.obs['labels_gex'] > 0]
        bdata.obs['labels_gex'] = bdata.obs['labels_gex'].astype(str)
        crop = b2c.get_crop(bdata, basis="array", mpp=mpp)
        rendered = b2c.view_labels(
            image_path="stardist/gex.tiff",
            labels_npz_path=f"stardist/gex_{prob_thresh}_{nms_thresh}.npz",
            crop=crop
        )

        plt.figure()
        plt.imshow(rendered)
        plt.axis("off")
        save_path = os.path.join(save_dir, f"{region_name}_{basename}.pdf")
        plt.tight_layout()
        plt.savefig(save_path, dpi=1200, bbox_inches="tight", pad_inches=0.5)
        plt.close()
        logging.info(f"Saved: {save_path}")


def visualize_joint_labels(adata, region_file, save_dir, mpp):
    """
    Visualize joint labels for each region.
    
    Args:
        adata: AnnData object
        region_file (str): Path to region coordinates file
        save_dir (str): Directory to save visualization results
        mpp (float): Microns per pixel
    """
    logging.info(f"Visualizing joint labels in {save_dir}...")
    regions = read_regions(region_file)
    os.makedirs(save_dir, exist_ok=True)
    basename = os.path.basename(save_dir.rstrip(os.sep))

    for region_name, region_coords in regions.items():
        mask = (
            (adata.obsm['spatial'][:, 0] >= region_coords["col_min"]) &
            (adata.obsm['spatial'][:, 0] <= region_coords["col_max"]) &
            (adata.obsm['spatial'][:, 1] >= region_coords["row_min"]) &
            (adata.obsm['spatial'][:, 1] <= region_coords["row_max"])
        )
        bdata = adata[mask]
        logging.info(f"{region_name} - Number of spots selected: {bdata.n_obs}")

        bdata = bdata[bdata.obs['labels_joint'] > 0]
        bdata.obs['labels_joint'] = bdata.obs['labels_joint'].astype(int)
        unique_labels_joint = sorted(bdata.obs['labels_joint'].unique())
        bdata.obs['labels_joint'] = pd.Categorical(
            bdata.obs['labels_joint'], categories=unique_labels_joint, ordered=True)
        colors_joint = [cc.glasbey[i % len(cc.glasbey)] for i in range(len(unique_labels_joint))]
        bdata.uns['labels_joint_colors'] = colors_joint

        if 'labels_joint_source' in bdata.obs.columns:
            unique_labels_source = sorted(bdata.obs['labels_joint_source'].unique())
            bdata.obs['labels_joint_source'] = pd.Categorical(
                bdata.obs['labels_joint_source'], categories=unique_labels_source, ordered=True)
            colors_source = ['#2c5ca4', '#ec5414']
            bdata.uns['labels_joint_source_colors'] = colors_source
        else:
            logging.error("Column 'labels_joint_source' not found in bdata.obs")
            raise ValueError("Column 'labels_joint_source' not found")

        titles = [f"{region_name} - Image", f"{region_name} - Joint Source", f"{region_name} - Joint Labels"]
        fig = sc.pl.spatial(
            bdata,
            color=[None, "labels_joint_source", "labels_joint"],
            img_key=f"{mpp}_mpp_150_buffer",
            basis="spatial_cropped_150_buffer",
            title=titles,
            show=False,
            legend_loc=None,
            return_fig=True
        )
        
        axes = fig.axes
        ax_source = axes[1]
        labels_source = bdata.obs['labels_joint_source'].cat.categories
        colors_source = bdata.uns['labels_joint_source_colors']
        handles_source = [Line2D([0], [0], marker='o', color='w',
                                 markerfacecolor=c, markersize=10) for c in colors_source]
        ax_source.legend(handles_source, labels_source, loc='center left',
                         bbox_to_anchor=(2, 0.5), ncol=1, fontsize=8)

        ax_joint = axes[2]
        labels_joint = bdata.obs['labels_joint'].cat.categories
        colors_joint = bdata.uns['labels_joint_colors']
        handles_joint = [Line2D([0], [0], marker='o', color='w',
                                markerfacecolor=c, markersize=10) for c in colors_joint]
        ax_joint.legend(handles_joint, labels_joint, loc='center left',
                        bbox_to_anchor=(3, 0.5), ncol=15, fontsize=8)

        save_path = os.path.join(save_dir, f"{region_name}_{basename}.pdf")
        plt.tight_layout()
        plt.savefig(save_path, dpi=1200, bbox_inches="tight", pad_inches=0.5)
        plt.close()
        logging.info(f"Saved: {save_path}")


def visualize_joint_all(cdata, region_file, save_dir, mpp, labels_key):
    """
    Visualize all joint labels for each region.
    
    Args:
        cdata: Cell-level AnnData object
        region_file (str): Path to region coordinates file
        save_dir (str): Directory to save visualization results
        mpp (float): Microns per pixel
        labels_key (str): Label key for analysis
    """
    logging.info(f"Visualizing all joint labels in {save_dir}...")
    regions = read_regions(region_file)
    os.makedirs(save_dir, exist_ok=True)
    basename = os.path.basename(save_dir.rstrip(os.sep))

    for region_name, region_coords in regions.items():
        spatial_coords = cdata.obsm['spatial']
        mask = (
            (spatial_coords[:, 0] >= region_coords["col_min"]) &
            (spatial_coords[:, 0] <= region_coords["col_max"]) &
            (spatial_coords[:, 1] >= region_coords["row_min"]) &
            (spatial_coords[:, 1] <= region_coords["row_max"])
        )
        ddata = cdata[mask]
        logging.info(f"{region_name} - Number of spots selected: {ddata.n_obs}")

        titles = [f"{region_name} - Bin Count", f"{region_name} - Joint Source"]
        
        # Conditional color assignment based on labels_key
        if labels_key != "labels_joint":
            color = [None, "bin_count"]
        else:
            color = ["bin_count", "labels_joint_source"]
            
        sc.pl.spatial(
            ddata,
            color=color,
            img_key=f"{mpp}_mpp_150_buffer",
            basis="spatial_cropped_150_buffer",
            color_map="OrRd",
            title=titles,
            size=0.5,
            show=False
        )

        save_path = os.path.join(save_dir, f"{region_name}_{basename}.pdf")
        plt.tight_layout()
        plt.savefig(save_path, dpi=1200, bbox_inches="tight", pad_inches=0.5)
        plt.close()
        logging.info(f"Saved: {save_path}")


def generate_roi_data(adata, cdata, region_file, output_dir):
    """
    Generate separate adata and cdata files for each ROI based on spatial coordinates.

    Args:
        adata: AnnData object containing spot-level data with spatial coordinates
        cdata: AnnData object containing cell-level data with spatial coordinates  
        region_file (str): Path to the file containing ROI coordinates
        output_dir (str): Directory to save the new adata and cdata files
    """
    regions = read_regions(region_file)
    logging.info(f"Generating ROI-specific adata and cdata for {len(regions)} regions...")

    for region_name, region_coords in regions.items():
        # Filter adata spots
        mask_adata = (
            (adata.obsm['spatial'][:, 0] >= region_coords["col_min"]) &
            (adata.obsm['spatial'][:, 0] <= region_coords["col_max"]) &
            (adata.obsm['spatial'][:, 1] >= region_coords["row_min"]) &
            (adata.obsm['spatial'][:, 1] <= region_coords["row_max"])
        )
        adata_roi = adata[mask_adata].copy()
        logging.info(f"{region_name} - Number of spots in adata: {adata_roi.n_obs}")

        # Filter cdata cells
        mask_cdata = (
            (cdata.obsm['spatial'][:, 0] >= region_coords["col_min"]) &
            (cdata.obsm['spatial'][:, 0] <= region_coords["col_max"]) &
            (cdata.obsm['spatial'][:, 1] >= region_coords["row_min"]) &
            (cdata.obsm['spatial'][:, 1] <= region_coords["row_max"])
        )
        cdata_roi = cdata[mask_cdata].copy()
        logging.info(f"{region_name} - Number of cells in cdata: {cdata_roi.n_obs}")

        # Create save directory
        save_dir = os.path.join(output_dir, "ROI_Data", region_name)
        os.makedirs(save_dir, exist_ok=True)

        # Save new adata and cdata files
        adata_roi.write_h5ad(os.path.join(save_dir, f"{region_name}_adata.h5ad"))
        cdata_roi.write_h5ad(os.path.join(save_dir, f"{region_name}_cdata.h5ad"))
        logging.info(f"Saved ROI {region_name} data to {save_dir}")


def process_spatial_data(args):
    """
    Main processing function for spatial transcriptomics data.
    
    Args:
        args: Command line arguments
        
    Returns:
        tuple: (adata, cdata) - processed AnnData objects
    """
    setup_logging(args.output_dir)
    logging.info("Starting spatial transcriptomics processing with bin2cell integration...")

    # Load and filter data using bin2cell
    adata = b2c.read_visium(args.path, source_image_path=args.source_image_path)
    adata.var_names_make_unique()
    logging.info(f"Loaded AnnData: {adata.n_obs} observations, {adata.n_vars} variables")
    
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_counts=1)
    logging.info(f"After filtering: {adata.n_obs} observations, {adata.n_vars} variables")

    # Get microns per pixel (MPP)
    library = list(adata.uns['spatial'].keys())[0]
    mpp = adata.uns['spatial'][library]['scalefactors']['microns_per_pixel']
    logging.info(f"MPP (microns per pixel): {mpp}")

    # Process HE image using bin2cell
    os.makedirs("stardist", exist_ok=True)
    b2c.scaled_he_image(adata, mpp=mpp, save_path="stardist/he.tiff")
    logging.info("Scaled HE image saved to stardist/he.tiff")

    # Execute destriping using bin2cell
    b2c.destripe(adata, adjust_counts=True)
    logging.info("Destriping completed")

    # Visualize destripe results
    destripe_and_save_regions(adata, args.region_file, "destripe", mpp)

    # Visualize NPZ labels
    visualize_labels_in_regions(args.npz_path, args.region_file, "npz_labels")

    # Insert QuPath labels using bin2cell
    b2c.insert_labels(adata, labels_npz_path=args.npz_path, basis="spatial",
                      spatial_key="spatial", mpp=mpp, labels_key="labels_qupath")
    adata.obs["labels_qupath"] = adata.obs["labels_qupath"].apply(
        lambda x: 0 if x == 1 else x - 1)
    logging.info(f"Inserted and adjusted QuPath labels: {adata.obs['labels_qupath'].unique()}")

    # Visualize and render QuPath labels
    visualize_after_insert_labels(adata, args.region_file, "segmentation", mpp)
    render_labels_after_insert(adata, args.region_file, "render_labels", 
                              args.npz_path, mpp, args.source_image_path)

    # Expand QuPath labels using bin2cell
    b2c.expand_labels(
        adata,
        labels_key='labels_qupath',
        expanded_labels_key="labels_qupath_expanded",
        algorithm=args.algorithm,
        max_bin_distance=args.max_bin_distance,
        volume_ratio=args.volume_ratio,
        k=args.k,
        subset_pca=args.subset_pca
    )
    logging.info("QuPath labels expanded")
    visualize_after_expansion(adata, args.region_file, "expanded_labels", mpp)

    # Generate GEX grid image using bin2cell
    b2c.grid_image(adata, "n_counts_adjusted", mpp=mpp,
                   sigma=5, save_path="stardist/gex.tiff")
    logging.info("GEX grid image generated")

    # Run StarDist using bin2cell
    b2c.stardist(
        image_path="stardist/gex.tiff",
        labels_npz_path=f"stardist/gex_{args.prob_thresh}_{args.nms_thresh}.npz",
        stardist_model="2D_versatile_fluo",
        prob_thresh=args.prob_thresh,
        nms_thresh=args.nms_thresh
    )
    logging.info("StarDist segmentation completed")

    # Insert GEX labels using bin2cell
    b2c.insert_labels(
        adata, 
        labels_npz_path=f"stardist/gex_{args.prob_thresh}_{args.nms_thresh}.npz", 
        basis="array", 
        mpp=mpp, 
        labels_key="labels_gex"
    )
    logging.info("GEX labels inserted")

    # Visualize and render GEX labels
    visualize_gex_labels(adata, args.region_file, "gex_labels", mpp)
    render_gex_labels(adata, args.region_file, "render_gex",
                      args.prob_thresh, args.nms_thresh, mpp)

    # Merge primary and secondary labels using bin2cell
    b2c.salvage_secondary_labels(adata, primary_label="labels_qupath_expanded",
                                 secondary_label="labels_gex", labels_key="labels_joint")
    logging.info("Salvaged secondary labels into joint labels")
    visualize_joint_labels(adata, args.region_file, "joint_labels", mpp)

    # Bin to cell conversion using bin2cell
    cdata = b2c.bin_to_cell(adata, labels_key=args.labels_key, 
                           spatial_keys=["spatial", "spatial_cropped_150_buffer"])
    
    if args.labels_key != "labels_joint":
        logging.warning("⚠️  Note: HE segmentation optimized, applying 'labels_qupath_expanded'")
    
    visualize_joint_all(cdata, args.region_file, "joint_labels_all", mpp, args.labels_key)

    # Save overall results
    save_dir = os.path.join(args.output_dir, "Data")
    os.makedirs(save_dir, exist_ok=True)
    cdata.write_h5ad(os.path.join(save_dir, f"{args.sample}_b2c.h5ad"))
    adata.write_h5ad(os.path.join(save_dir, f"{args.sample}_2um.h5ad"))
    logging.info(f"Results saved to {save_dir}")
    
    # Generate and save ROI-specific data
    generate_roi_data(adata, cdata, args.region_file, args.output_dir)

    logging.info("Spatial processing completed successfully! ✨")
    return adata, cdata


def main():
    """Main entry point for spatial processing."""
    parser = argparse.ArgumentParser(
        description="Process spatial transcriptomics data using bin2cell integration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python spatial_processor.py --path /data/visium --source_image_path /data/image.tif \\
    --region_file regions.txt --npz_path labels.npz --output_dir results \\
    --sample E14.5

  python spatial_processor.py --path /data/visium --source_image_path /data/image.tif \\
    --region_file regions.txt --npz_path labels.npz --output_dir results \\
    --sample P3 --algorithm volume_ratio --volume_ratio 3.0
        """
    )
    
    # Required arguments
    parser.add_argument("--path", required=True,
                        help="Path to Visium data directory")
    parser.add_argument("--source_image_path", required=True,
                        help="Path to source image file")
    parser.add_argument("--region_file", required=True,
                        help="Path to region coordinates file")
    parser.add_argument("--npz_path", required=True,
                        help="Path to QuPath NPZ labels file")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to save output files")
    parser.add_argument("--sample", required=True,
                        help="Sample type")

    # StarDist parameters
    parser.add_argument("--prob_thresh", type=float, default=0.05,
                        help="Probability threshold for StarDist (default: 0.05)")
    parser.add_argument("--nms_thresh", type=float, default=0.5,
                        help="NMS threshold for StarDist (default: 0.5)")
    
    # Labels parameters
    parser.add_argument("--labels_key", default="labels_joint",
                        help="Labels key for bin_to_cell (default: 'labels_joint')")

    # Label expansion parameters
    parser.add_argument("--algorithm", default="max_bin_distance", 
                        choices=["max_bin_distance", "volume_ratio"],
                        help="Expansion algorithm: 'max_bin_distance' or 'volume_ratio' (default: 'max_bin_distance')")
    parser.add_argument("--max_bin_distance", type=int, default=2,
                        help="Maximum bin distance (default: 2)")
    parser.add_argument("--volume_ratio", type=float, default=4.0,
                        help="Volume ratio (default: 4.0)")
    parser.add_argument("--k", type=int, default=4,
                        help="Number of nearest neighbors (default: 4)")
    parser.add_argument("--subset_pca", type=str, default="True",
                        choices=["True", "False"],
                        help="Whether to compute PCA only for mean bins (default: 'True')")

    args = parser.parse_args()

    # Convert subset_pca from string to boolean
    args.subset_pca = args.subset_pca == "True"

    try:
        adata, cdata = process_spatial_data(args)
        print("✓ Spatial processing completed successfully!")
        
    except Exception as e:
        logging.error(f"Processing failed: {e}")
        raise


if __name__ == "__main__":
    main()# Version: 1.0.7
