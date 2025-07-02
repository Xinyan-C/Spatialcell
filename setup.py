#!/usr/bin/env python3
"""
Setup script for SpatialCell: Integrated Spatial Transcriptomics Analysis Pipeline
"""

from setuptools import setup, find_packages
import os

def read_requirements():
    """Read requirements from requirements.txt"""
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    with open(requirements_path, 'r', encoding='utf-8') as f:
        requirements = []
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                requirements.append(line)
    return requirements

def read_long_description():
    """Read long description from README.md"""
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ""

setup(
    name="spatialcell",
    version="1.0.0",
    author="Xinyan",
    author_email="keepandon@gmail.com",
    description="Integrated pipeline for spatial transcriptomics cell segmentation and annotation using QuPath, Bin2cell, and TopAct",
    long_description=read_long_description(),
    long_description_content_type="text/markdown",
    url="https://github.com/Xinyan-C/Spatialcell",
    project_urls={
        "Bug Reports": "https://github.com/Xinyan-C/Spatialcell/issues",
        "Source": "https://github.com/Xinyan-C/Spatialcell",
        "Documentation": "https://github.com/Xinyan-C/Spatialcell#readme",
    },
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "black>=22.0",
            "flake8>=4.0",
            "mypy>=0.950",
        ],
        "docs": [
            "sphinx>=4.0",
            "sphinx-rtd-theme>=1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "spatialcell=spatialcell.workflows.main:main",
        ],
    },
    keywords=[
        "spatial-transcriptomics",
        "cell-segmentation", 
        "cell-annotation",
        "QuPath",
        "Bin2cell",
        "TopAct",
        "bioinformatics",
        "single-cell",
    ],
    zip_safe=False,
    include_package_data=True,
    package_data={
        "spatialcell": [
            "qupath_scripts/*.groovy",
            "examples/*.py",
            "examples/*.yml",
            "utils/*.py",
        ],
    },
)
