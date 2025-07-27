"""Xinyan-recode-Classes and methods for array-based spatial transcriptomics analysis with float coordinates support."""

import itertools
from typing import Iterable, Iterator, Sequence, cast, Callable, Tuple
from multiprocessing import Process, Queue
from warnings import simplefilter

import numpy as np
import numpy.typing as npt
from scipy import sparse
import pandas as pd
from topact.countdata import CountTable
from topact.classifier import Classifier
from topact import densetools

# 辅助函数
def combine_coords(coords: Iterable[float]) -> str:
    """将浮点坐标元组组合成唯一字符串标识符。"""
    return ','.join(map(str, coords))

def split_coords(ident: str) -> Tuple[float, ...]:
    """将唯一字符串标识符拆分为对应的浮点坐标。"""
    return tuple(map(float, ident.split(','))) if ident else ()

def first_coord(ident: str) -> float:
    """从唯一标识符中获取第一个浮点坐标。"""
    return split_coords(ident)[0]

def second_coord(ident: str) -> float:
    """从唯一标识符中获取第二个浮点坐标。"""
    return split_coords(ident)[1]

def extract_classifications(confidence_matrix: npt.NDArray,
                            threshold: float
                            ) -> dict[Tuple[float, float], int]:
    """根据阈值从置信度矩阵中提取所有点的分类字典。"""
    confident = zip(*np.where(confidence_matrix.max(axis=-1) >= threshold))
    confident = cast(Iterator[Tuple[int, int]], confident)
    classifications: dict[Tuple[float, float], int] = {}
    for i, j in confident:
        cell_type = np.argmax(confidence_matrix[i, j])
        classifications[(i, j)] = cell_type
    return classifications

def extract_image(confidence_matrix: npt.NDArray,
                  threshold: float
                  ) -> npt.NDArray:
    """根据阈值从置信度矩阵中提取图像。"""
    classifications = extract_classifications(confidence_matrix, threshold)
    image = np.empty(confidence_matrix.shape[:2])
    image[:] = np.nan
    for (i, j), c in classifications.items():
        image[i, j] = c
    return image

# ExpressionGrid 类
class ExpressionGrid:
    """带有基因表达的空间网格类，支持浮点坐标。"""
    def __init__(self,
                 table: pd.DataFrame,
                 genes: Sequence[str],
                 gene_col: str = "gene",
                 count_col: str = "count"
                 ):
        """从数据框初始化网格，包含浮点坐标的表达数据。"""
        self.table = table
        self.genes = genes
        self.gene_col = gene_col
        self.count_col = count_col
        self.x_min, self.x_max = table.x.min(), table.x.max()
        self.y_min, self.y_max = table.y.min(), table.y.max()
        self.num_genes = len(genes)

    def expression_vec(self, sub_table: pd.DataFrame) -> sparse.spmatrix:
        """返回子表中所有数据点的总基因表达向量。"""
        expr = np.zeros(self.num_genes)
        for _, row in sub_table.iterrows():
            gene = row[self.gene_col]
            count = row[self.count_col]
            expr[self.genes.index(gene)] += count
        return sparse.csr_matrix(expr)

    def square_nbhd(self, x: float, y: float, scale: float) -> pd.DataFrame:
        """返回指定点 (x, y) 在给定尺度内的方形邻域内的所有实际数据点。"""
        x_min = x - scale
        x_max = x + scale
        y_min = y - scale
        y_max = y + scale
        sub_table = self.table[(self.table.x >= x_min) & (self.table.x <= x_max) &
                               (self.table.y >= y_min) & (self.table.y <= y_max)]
        return sub_table

# Worker 类
class Worker(Process):
    def __init__(self,
                 grid: ExpressionGrid,
                 count_grid: 'CountGrid',
                 min_scale: float,
                 max_scale: float,
                 scale_step: float,  # 新增参数
                 classifier: Classifier,
                 job_queue: Queue,
                 res_queue: Queue,
                 procid: int,
                 verbose: bool,
                 neighborhood_func: Callable[['CountGrid', Tuple[float, float], float], pd.DataFrame] = None
                 ):
        super().__init__()
        self.grid = grid
        self.count_grid = count_grid
        self.min_scale = min_scale
        self.max_scale = max_scale
        self.scale_step = scale_step  # 使用传入的步长
        self.classifier = classifier
        self.job_queue = job_queue
        self.res_queue = res_queue
        self.procid = procid
        self.verbose = verbose
        self.neighborhood_func = neighborhood_func

    def run(self):
        simplefilter(action='ignore', category=FutureWarning)
        if self.verbose:
            print(f'Worker {self.procid} started')

        num_classes = len(self.classifier.classes)
        scales = np.arange(self.min_scale, self.max_scale + self.scale_step, self.scale_step)  # 使用动态步长
        num_scales = len(scales)
        exprs = np.zeros((num_scales, self.grid.num_genes))

        for center_x, center_y in iter(self.job_queue.get, None):
            if self.verbose:
                print(f"Worker {self.procid} got job ({center_x}, {center_y})")
            for scale_idx, scale in enumerate(scales):
                if self.neighborhood_func is not None:
                    sub_table = self.neighborhood_func(self.count_grid, (center_x, center_y), scale)
                else:
                    sub_table = self.grid.square_nbhd(center_x, center_y, scale)
                expr = self.grid.expression_vec(sub_table)
                exprs[scale_idx] = expr.toarray()[0]

            first_nonzero = densetools.first_nonzero_1d(exprs.sum(axis=1))
            probs = np.empty((num_scales, num_classes))
            probs[:] = -1

            if 0 <= first_nonzero < num_scales:
                to_classify = np.vstack(exprs[first_nonzero:])
                all_confidences = self.classifier.classify(to_classify)
                probs[first_nonzero:] = all_confidences

            self.res_queue.put((center_x, center_y, probs.tolist()))
        self.res_queue.put(None)
        if self.verbose:
            print(f'Worker {self.procid} finished')

# CountGrid 类
class CountGrid(CountTable):
    """支持浮点坐标的空间转录组学对象及其相关方法。"""
    def __init__(self, *args, **kwargs):
        """从数据框初始化空间数据。"""
        super().__init__(*args, **kwargs)
        self.generate_expression_grid()

    @classmethod
    def from_coord_table(cls, table, **kwargs):
        """从包含浮点坐标的表格创建 CountGrid。"""
        samples = table[['x', 'y']].drop_duplicates().apply(lambda row: combine_coords((row['x'], row['y'])), axis=1)
        new_table = table.copy()
        new_table['sample'] = new_table.apply(lambda row: combine_coords((row['x'], row['y'])), axis=1)
        count_grid = cls(new_table, samples=list(samples), **kwargs)
        samples = count_grid.samples
        x_coords = {sample: first_coord(sample) for sample in samples}
        y_coords = {sample: second_coord(sample) for sample in samples}
        count_grid.add_metadata('x', x_coords)
        count_grid.add_metadata('y', y_coords)
        return count_grid

    def pseudobulk(self) -> npt.NDArray:
        """计算所有样本的伪批量表达。"""
        return self.table.groupby('gene')[self.count_col].sum().reindex(self.genes, fill_value=0).values

    def count_matrix(self) -> npt.NDArray:
        """此方法暂未针对浮点坐标调整，保留占位符。"""
        raise NotImplementedError("count_matrix needs adjustment for float coordinates.")

    def density_mask(self, radius: float, threshold: int) -> npt.NDArray:
        """此方法暂未针对浮点坐标调整，保留占位符。"""
        raise NotImplementedError("density_mask needs adjustment for float coordinates.")

    def generate_expression_grid(self):
        """生成支持浮点坐标的 ExpressionGrid。"""
        self.grid = ExpressionGrid(self.table,
                                   genes=self.genes,
                                   gene_col=self.gene_col,
                                   count_col=self.count_col)

    def classify_parallel(self,
                          classifier: Classifier,
                          min_scale: float,  # in μm
                          max_scale: float,  # in μm
                          outfile: str,
                          mpp: float,  # microns per pixel
                          mask: npt.NDArray | None = None,
                          num_proc: int = 1,
                          verbose: bool = False,
                          neighborhood_func: Callable[['CountGrid', Tuple[float, float], float], pd.DataFrame] = None
                          ) -> npt.NDArray:
        """并行分类方法，支持浮点坐标，scale 参数以μm为单位，通过 MPP 转换为像素单位。"""
        # 将尺度从 μm 转换为像素单位
        min_scale_pixels = min_scale / mpp
        max_scale_pixels = max_scale / mpp
        scale_step_pixels = 2.0 / mpp  # 默认步长 2μm，转换为像素单位

        # 获取所有独特坐标点作为分类中心
        unique_coords = self.table[['x', 'y']].drop_duplicates().values
        num_points = len(unique_coords)
        num_classes = len(classifier.classes)
        scales = np.arange(min_scale_pixels, max_scale_pixels + scale_step_pixels, scale_step_pixels)
        num_scales = len(scales)

        # 初始化结果数组
        confidence_matrix = np.zeros((num_points, num_scales, num_classes))

        # 创建任务和结果队列
        job_queue = Queue()
        res_queue = Queue()

        # 将任务放入队列
        for x, y in unique_coords:
            job_queue.put((x, y))
        for _ in range(num_proc):
            job_queue.put(None)

        # 启动工作进程，传递 scale_step_pixels
        workers = []
        for i in range(num_proc):
            worker = Worker(self.grid, self, min_scale_pixels, max_scale_pixels, scale_step_pixels, 
                            classifier, job_queue, res_queue, i, verbose, neighborhood_func)
            workers.append(worker)
            worker.start()

        # 收集结果
        finished = 0
        coord_to_idx = {(x, y): i for i, (x, y) in enumerate(unique_coords)}
        while finished < num_proc:
            result = res_queue.get()
            if result is None:
                finished += 1
            else:
                center_x, center_y, probs = result
                idx = coord_to_idx[(center_x, center_y)]
                confidence_matrix[idx] = probs

        # 保存结果
        np.save(outfile, confidence_matrix)
    
        # 等待所有进程结束
        for worker in workers:
            worker.join()

        return confidence_matrix

    def annotate(self,
                 confidence_matrix: npt.NDArray,
                 threshold: float,
                 labels: Tuple[str, ...],
                 column_label: str = "cell type"):
        """根据置信度矩阵为数据点添加分类标注。"""
        classifications = extract_classifications(confidence_matrix, threshold)
        annotated_table = self.table.copy()
        annotated_table[column_label] = np.nan

        # 将分类结果映射到表格
        for (x, y), cell_type_idx in classifications.items():
            sample_id = combine_coords((x, y))
            annotated_table.loc[annotated_table['sample'] == sample_id, column_label] = labels[cell_type_idx]

        self.table = annotated_table
        return annotated_table