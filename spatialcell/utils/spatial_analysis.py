#!/usr/bin/env python3
"""
Enhanced Spatial Analysis Module for Spatialcell Pipeline

This module extends TopAct framework functionality to support array-based spatial 
transcriptomics analysis with floating-point coordinates, enabling high-resolution
spatial cell type annotation.

Built upon TopAct framework: https://gitlab.com/kfbenjamin/topact.git
Enhanced for Spatialcell pipeline with float coordinate support.

Author: Xinyan
License: MIT
"""

import itertools
from typing import Iterable, Iterator, Sequence, cast, Callable, Tuple, Dict, Optional, Any
from multiprocessing import Process, Queue
from warnings import simplefilter
import logging

import numpy as np
import numpy.typing as npt
from scipy import sparse
import pandas as pd
from topact.countdata import CountTable  # Using TopAct framework
from topact.classifier import Classifier  # Using TopAct classifiers
from topact import densetools  # Using TopAct utilities


# =============================================================================
# Utility Functions for Float Coordinate Handling
# =============================================================================

def combine_coordinates(coords: Iterable[float]) -> str:
    """
    Combine float coordinate tuple into unique string identifier.
    
    Args:
        coords (Iterable[float]): Coordinate values to combine
        
    Returns:
        str: Unique string identifier for coordinates
        
    Example:
        >>> combine_coordinates([123.45, 678.90])
        '123.45,678.9'
    """
    return ','.join(map(str, coords))


def split_coordinates(identifier: str) -> Tuple[float, ...]:
    """
    Split unique string identifier back to corresponding float coordinates.
    
    Args:
        identifier (str): String identifier to split
        
    Returns:
        tuple: Float coordinates
        
    Example:
        >>> split_coordinates('123.45,678.9')
        (123.45, 678.9)
    """
    return tuple(map(float, identifier.split(','))) if identifier else ()


def get_x_coordinate(identifier: str) -> float:
    """
    Extract X coordinate (first coordinate) from unique identifier.
    
    Args:
        identifier (str): Coordinate identifier string
        
    Returns:
        float: X coordinate value
    """
    return split_coordinates(identifier)[0]


def get_y_coordinate(identifier: str) -> float:
    """
    Extract Y coordinate (second coordinate) from unique identifier.
    
    Args:
        identifier (str): Coordinate identifier string
        
    Returns:
        float: Y coordinate value
    """
    return split_coordinates(identifier)[1]


def extract_classifications_from_confidence(confidence_matrix: npt.NDArray,
                                          threshold: float) -> Dict[Tuple[float, float], int]:
    """
    Extract classification dictionary for all points from confidence matrix based on threshold.
    
    Args:
        confidence_matrix (np.ndarray): Classification confidence scores
        threshold (float): Confidence threshold for classification
        
    Returns:
        dict: Mapping of coordinates to cell type classifications
    """
    confident_indices = zip(*np.where(confidence_matrix.max(axis=-1) >= threshold))
    confident_indices = cast(Iterator[Tuple[int, int]], confident_indices)
    
    classifications: Dict[Tuple[float, float], int] = {}
    for i, j in confident_indices:
        cell_type = np.argmax(confidence_matrix[i, j])
        classifications[(i, j)] = cell_type
        
    return classifications


def extract_classification_image(confidence_matrix: npt.NDArray,
                               threshold: float) -> npt.NDArray:
    """
    Extract classification image from confidence matrix based on threshold.
    
    Args:
        confidence_matrix (np.ndarray): Classification confidence scores
        threshold (float): Confidence threshold
        
    Returns:
        np.ndarray: 2D classification image with NaN for unclassified regions
    """
    classifications = extract_classifications_from_confidence(confidence_matrix, threshold)
    
    image = np.empty(confidence_matrix.shape[:2])
    image[:] = np.nan
    
    for (i, j), classification in classifications.items():
        image[i, j] = classification
        
    return image


# =============================================================================
# Enhanced Expression Grid Class
# =============================================================================

class FloatExpressionGrid:
    """
    Spatial grid class with gene expression support for floating-point coordinates.
    
    This class handles spatial transcriptomics data with precise floating-point
    coordinates, enabling high-resolution spatial analysis.
    """
    
    def __init__(self,
                 table: pd.DataFrame,
                 genes: Sequence[str],
                 gene_col: str = "gene",
                 count_col: str = "counts") -> None:
        """
        Initialize grid from DataFrame containing expression data with float coordinates.
        
        Args:
            table (pd.DataFrame): Expression data with 'x', 'y' coordinates
            genes (Sequence[str]): List of gene names
            gene_col (str): Column name for gene identifiers
            count_col (str): Column name for expression counts
        """
        self.table = table.copy()
        self.genes = list(genes)
        self.gene_col = gene_col
        self.count_col = count_col
        
        # Calculate spatial bounds
        self.x_min, self.x_max = table['x'].min(), table['x'].max()
        self.y_min, self.y_max = table['y'].min(), table['y'].max()
        self.num_genes = len(genes)
        
        # Create gene index mapping for efficient lookup
        self.gene_to_index = {gene: i for i, gene in enumerate(genes)}
        
        logging.debug(f"Initialized FloatExpressionGrid: {self.num_genes} genes, "
                     f"spatial bounds: X[{self.x_min:.2f}, {self.x_max:.2f}], "
                     f"Y[{self.y_min:.2f}, {self.y_max:.2f}]")
    
    def get_expression_vector(self, sub_table: pd.DataFrame) -> sparse.spmatrix:
        """
        Calculate total gene expression vector for all data points in sub-table.
        
        Args:
            sub_table (pd.DataFrame): Subset of expression data
            
        Returns:
            sparse.spmatrix: Sparse expression vector
        """
        expression = np.zeros(self.num_genes, dtype=np.float32)
        
        for _, row in sub_table.iterrows():
            gene = row[self.gene_col]
            count = row[self.count_col]
            
            if gene in self.gene_to_index:
                gene_idx = self.gene_to_index[gene]
                expression[gene_idx] += count
        
        return sparse.csr_matrix(expression)
    
    def get_square_neighborhood(self, x: float, y: float, scale: float) -> pd.DataFrame:
        """
        Get all actual data points within square neighborhood of specified point and scale.
        
        Args:
            x (float): Center X coordinate
            y (float): Center Y coordinate  
            scale (float): Neighborhood radius (half-width of square)
            
        Returns:
            pd.DataFrame: Data points within neighborhood
        """
        x_min = x - scale
        x_max = x + scale
        y_min = y - scale
        y_max = y + scale
        
        neighborhood_mask = (
            (self.table['x'] >= x_min) & (self.table['x'] <= x_max) &
            (self.table['y'] >= y_min) & (self.table['y'] <= y_max)
        )
        
        return self.table[neighborhood_mask]
    
    def get_circular_neighborhood(self, x: float, y: float, radius: float) -> pd.DataFrame:
        """
        Get all data points within circular neighborhood of specified point.
        
        Args:
            x (float): Center X coordinate
            y (float): Center Y coordinate
            radius (float): Neighborhood radius
            
        Returns:
            pd.DataFrame: Data points within circular neighborhood
        """
        distances = np.sqrt((self.table['x'] - x)**2 + (self.table['y'] - y)**2)
        neighborhood_mask = distances <= radius
        
        return self.table[neighborhood_mask]


# =============================================================================
# Parallel Worker Process
# =============================================================================

class SpatialClassificationWorker(Process):
    """
    Worker process for parallel spatial classification with TopAct framework.
    
    Handles multi-scale classification of spatial points using custom neighborhood
    functions and floating-point coordinates.
    """
    
    def __init__(self,
                 expression_grid: FloatExpressionGrid,
                 count_grid: 'EnhancedCountGrid',
                 min_scale: float,
                 max_scale: float,
                 scale_step: float,
                 classifier: Classifier,
                 job_queue: Queue,
                 result_queue: Queue,
                 process_id: int,
                 verbose: bool = False,
                 neighborhood_func: Optional[Callable] = None) -> None:
        """
        Initialize worker process.
        
        Args:
            expression_grid (FloatExpressionGrid): Expression grid object
            count_grid (EnhancedCountGrid): Count grid for neighborhood functions
            min_scale (float): Minimum scale for analysis
            max_scale (float): Maximum scale for analysis  
            scale_step (float): Step size between scales
            classifier (Classifier): TopAct classifier
            job_queue (Queue): Job queue for coordinates to process
            result_queue (Queue): Result queue for classification outputs
            process_id (int): Worker process identifier
            verbose (bool): Enable verbose logging
            neighborhood_func (Callable, optional): Custom neighborhood function
        """
        super().__init__()
        self.expression_grid = expression_grid
        self.count_grid = count_grid
        self.min_scale = min_scale
        self.max_scale = max_scale
        self.scale_step = scale_step
        self.classifier = classifier
        self.job_queue = job_queue
        self.result_queue = result_queue
        self.process_id = process_id
        self.verbose = verbose
        self.neighborhood_func = neighborhood_func
        
    def run(self) -> None:
        """
        Main worker process execution.
        
        Processes coordinates from job queue and performs multi-scale classification
        using TopAct framework with custom neighborhood constraints.
        """
        # Suppress warnings in worker processes
        simplefilter(action='ignore', category=FutureWarning)
        
        if self.verbose:
            print(f'Spatial classification worker {self.process_id} started')
        
        num_classes = len(self.classifier.classes)
        scales = np.arange(self.min_scale, self.max_scale + self.scale_step, self.scale_step)
        num_scales = len(scales)
        
        # Pre-allocate expression array for efficiency
        expressions = np.zeros((num_scales, self.expression_grid.num_genes), dtype=np.float32)
        
        # Process jobs until sentinel None is received
        for center_x, center_y in iter(self.job_queue.get, None):
            if self.verbose:
                print(f"Worker {self.process_id} processing point ({center_x:.2f}, {center_y:.2f})")
            
            # Calculate expression at each scale
            for scale_idx, scale in enumerate(scales):
                if self.neighborhood_func is not None:
                    # Use custom neighborhood function (e.g., cell-constrained)
                    sub_table = self.neighborhood_func(self.count_grid, (center_x, center_y), scale)
                else:
                    # Use default square neighborhood
                    sub_table = self.expression_grid.get_square_neighborhood(center_x, center_y, scale)
                
                expression_vector = self.expression_grid.get_expression_vector(sub_table)
                expressions[scale_idx] = expression_vector.toarray()[0]
            
            # Find first non-zero scale and classify
            first_nonzero = densetools.first_nonzero_1d(expressions.sum(axis=1))
            probabilities = np.full((num_scales, num_classes), -1.0)
            
            if 0 <= first_nonzero < num_scales:
                # Classify scales with non-zero expression
                to_classify = expressions[first_nonzero:]
                classification_confidences = self.classifier.classify(to_classify)
                probabilities[first_nonzero:] = classification_confidences
            
            # Send result back to main process
            self.result_queue.put((center_x, center_y, probabilities.tolist()))
        
        # Signal completion
        self.result_queue.put(None)
        
        if self.verbose:
            print(f'Spatial classification worker {self.process_id} finished')


# =============================================================================
# Enhanced Count Grid Class
# =============================================================================

class EnhancedCountGrid(CountTable):
    """
    Enhanced spatial transcriptomics object with float coordinate support and TopAct integration.
    
    Extends TopAct's CountTable with high-resolution spatial analysis capabilities,
    supporting floating-point coordinates and custom neighborhood functions.
    """
    
    def __init__(self, *args, **kwargs) -> None:
        """
        Initialize spatial data from DataFrame with TopAct integration.
        
        Args:
            *args: Arguments passed to CountTable
            **kwargs: Keyword arguments passed to CountTable
        """
        super().__init__(*args, **kwargs)
        self.generate_expression_grid()
        
    @classmethod
    def from_coordinate_table(cls, 
                            table: pd.DataFrame, 
                            **kwargs) -> 'EnhancedCountGrid':
        """
        Create EnhancedCountGrid from table containing float coordinates.
        
        Args:
            table (pd.DataFrame): Data with 'x', 'y' columns and expression data
            **kwargs: Additional arguments for CountTable initialization
            
        Returns:
            EnhancedCountGrid: Initialized count grid object
        """
        # Create sample identifiers from coordinates
        coordinate_samples = table[['x', 'y']].drop_duplicates()
        sample_ids = coordinate_samples.apply(
            lambda row: combine_coordinates((row['x'], row['y'])), axis=1
        )
        
        # Add sample column to table
        enhanced_table = table.copy()
        enhanced_table['sample'] = enhanced_table.apply(
            lambda row: combine_coordinates((row['x'], row['y'])), axis=1
        )
        
        # Initialize count grid
        count_grid = cls(enhanced_table, samples=list(sample_ids), **kwargs)
        
        # Add coordinate metadata
        samples = count_grid.samples
        x_coordinates = {sample: get_x_coordinate(sample) for sample in samples}
        y_coordinates = {sample: get_y_coordinate(sample) for sample in samples}
        
        count_grid.add_metadata('x', x_coordinates)
        count_grid.add_metadata('y', y_coordinates)
        
        logging.info(f"Created EnhancedCountGrid with {len(samples)} unique spatial locations")
        
        return count_grid
    
    def calculate_pseudobulk(self) -> npt.NDArray:
        """
        Calculate pseudobulk expression across all samples.
        
        Returns:
            np.ndarray: Pseudobulk expression vector
        """
        return self.table.groupby('gene')[self.count_col].sum().reindex(
            self.genes, fill_value=0
        ).values
    
    def get_count_matrix(self) -> npt.NDArray:
        """
        Generate count matrix for float coordinates.
        
        Note: This method needs adjustment for float coordinates.
        Currently raises NotImplementedError as placeholder.
        
        Raises:
            NotImplementedError: Method needs implementation for float coordinates
        """
        raise NotImplementedError(
            "Count matrix generation needs adjustment for float coordinates. "
            "Consider using pseudobulk or spatial aggregation methods."
        )
    
    def create_density_mask(self, radius: float, threshold: int) -> npt.NDArray:
        """
        Create density mask for float coordinates.
        
        Note: This method needs adjustment for float coordinates.
        Currently raises NotImplementedError as placeholder.
        
        Args:
            radius (float): Neighborhood radius
            threshold (int): Density threshold
            
        Raises:
            NotImplementedError: Method needs implementation for float coordinates
        """
        raise NotImplementedError(
            "Density mask creation needs adjustment for float coordinates. "
            "Consider implementing spatial density estimation methods."
        )
    
    def generate_expression_grid(self) -> None:
        """
        Generate FloatExpressionGrid for high-resolution spatial analysis.
        """
        self.grid = FloatExpressionGrid(
            self.table,
            genes=self.genes,
            gene_col=self.gene_col,
            count_col=self.count_col
        )
        logging.debug("Generated FloatExpressionGrid for spatial analysis")
    
    def classify_parallel(self,
                         classifier: Classifier,
                         min_scale: float,
                         max_scale: float,
                         outfile: str,
                         mpp: float,
                         mask: Optional[npt.NDArray] = None,
                         num_processes: int = 1,
                         verbose: bool = False,
                         neighborhood_func: Optional[Callable] = None) -> npt.NDArray:
        """
        Perform parallel spatial classification with TopAct framework.
        
        Supports floating-point coordinates and custom neighborhood functions.
        Scale parameters are in micrometers and converted to pixels using MPP.
        
        Args:
            classifier (Classifier): Trained TopAct classifier
            min_scale (float): Minimum neighborhood scale in micrometers
            max_scale (float): Maximum neighborhood scale in micrometers
            outfile (str): Output file path for results
            mpp (float): Microns per pixel conversion factor
            mask (np.ndarray, optional): Spatial mask for analysis
            num_processes (int): Number of parallel processes
            verbose (bool): Enable verbose logging
            neighborhood_func (Callable, optional): Custom neighborhood function
            
        Returns:
            np.ndarray: Classification confidence matrix
        """
        # Convert scales from micrometers to pixels
        min_scale_pixels = min_scale / mpp
        max_scale_pixels = max_scale / mpp
        scale_step_pixels = 2.0 / mpp  # Default 2μm step converted to pixels
        
        logging.info(f"Starting parallel classification:")
        logging.info(f"  Scale range: {min_scale:.1f}-{max_scale:.1f} μm "
                    f"({min_scale_pixels:.1f}-{max_scale_pixels:.1f} pixels)")
        logging.info(f"  Processes: {num_processes}")
        
        # Get unique coordinates for classification centers
        unique_coordinates = self.table[['x', 'y']].drop_duplicates().values
        num_points = len(unique_coordinates)
        num_classes = len(classifier.classes)
        
        scales = np.arange(min_scale_pixels, max_scale_pixels + scale_step_pixels, scale_step_pixels)
        num_scales = len(scales)
        
        logging.info(f"  Points: {num_points}, Classes: {num_classes}, Scales: {num_scales}")
        
        # Initialize result matrix
        confidence_matrix = np.zeros((num_points, num_scales, num_classes), dtype=np.float32)
        
        # Create job and result queues
        job_queue = Queue()
        result_queue = Queue()
        
        # Add jobs to queue
        for x, y in unique_coordinates:
            job_queue.put((x, y))
        
        # Add sentinel values for workers
        for _ in range(num_processes):
            job_queue.put(None)
        
        # Start worker processes
        workers = []
        for i in range(num_processes):
            worker = SpatialClassificationWorker(
                self.grid, self, min_scale_pixels, max_scale_pixels, scale_step_pixels,
                classifier, job_queue, result_queue, i, verbose, neighborhood_func
            )
            workers.append(worker)
            worker.start()
        
        # Collect results
        finished_workers = 0
        coordinate_to_index = {(x, y): i for i, (x, y) in enumerate(unique_coordinates)}
        
        while finished_workers < num_processes:
            result = result_queue.get()
            if result is None:
                finished_workers += 1
            else:
                center_x, center_y, probabilities = result
                point_index = coordinate_to_index[(center_x, center_y)]
                confidence_matrix[point_index] = probabilities
        
        # Save results
        np.save(outfile, confidence_matrix)
        logging.info(f"Classification results saved to: {outfile}")
        
        # Wait for all processes to complete
        for worker in workers:
            worker.join()
        
        return confidence_matrix
    
    def annotate_with_classifications(self,
                                    confidence_matrix: npt.NDArray,
                                    threshold: float,
                                    cell_type_labels: Tuple[str, ...],
                                    column_name: str = "predicted_cell_type") -> pd.DataFrame:
        """
        Annotate data points with cell type classifications based on confidence matrix.
        
        Args:
            confidence_matrix (np.ndarray): Classification confidence scores
            threshold (float): Confidence threshold for classification
            cell_type_labels (tuple): Cell type label names
            column_name (str): Column name for annotations
            
        Returns:
            pd.DataFrame: Annotated table with cell type predictions
        """
        classifications = extract_classifications_from_confidence(confidence_matrix, threshold)
        
        annotated_table = self.table.copy()
        annotated_table[column_name] = np.nan
        
        # Map classifications to table
        for (x, y), cell_type_index in classifications.items():
            sample_identifier = combine_coordinates((x, y))
            cell_type_label = cell_type_labels[cell_type_index]
            
            sample_mask = annotated_table['sample'] == sample_identifier
            annotated_table.loc[sample_mask, column_name] = cell_type_label
        
        # Update internal table
        self.table = annotated_table
        
        logging.info(f"Annotated {len(classifications)} points with cell type predictions")
        logging.info(f"Classification distribution:\n{annotated_table[column_name].value_counts()}")
        
        return annotated_table


# =============================================================================
# Module Aliases for Backward Compatibility
# =============================================================================

# Create aliases for the enhanced classes to maintain API compatibility
CountGrid = EnhancedCountGrid
ExpressionGrid = FloatExpressionGrid
Worker = SpatialClassificationWorker

# Export main functions
extract_classifications = extract_classifications_from_confidence
extract_image = extract_classification_image
combine_coords = combine_coordinates
split_coords = split_coordinates
first_coord = get_x_coordinate
second_coord = get_y_coordinate


# =============================================================================
# Module Information
# =============================================================================

__version__ = "1.0.0"
__author__ = "Xinyan"
__all__ = [
    'EnhancedCountGrid', 'FloatExpressionGrid', 'SpatialClassificationWorker',
    'CountGrid', 'ExpressionGrid', 'Worker',  # Aliases for compatibility
    'combine_coordinates', 'split_coordinates', 
    'get_x_coordinate', 'get_y_coordinate',
    'extract_classifications_from_confidence', 'extract_classification_image'
]

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logging.info("Enhanced Spatial Analysis Module for Spatialcell Pipeline")
    logging.info(f"Version: {__version__}")
    logging.info("Built upon TopAct framework with float coordinate support")# Version: 1.0.7
