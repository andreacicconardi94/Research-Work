# SEM Patch Extractor

A lightweight Python tool to extract **N µm × N µm** square patches from SEM `.tif` images.

Supports two modes:
- **Scalebar mode**: detect 1 µm length from a reference scalebar image.
- **Metadata mode**: read the pixel size from embedded TIFF metadata.

## Features

- Automatic detection of scalebar length (in pixels)
- Patch extraction with configurable physical size (in µm)
- Metadata-based patching (via TIFF `Image Pixel Size`)
- Automatic margin cropping (removes dark footers/logos)
- Optional overlap and filtering
- Easy integration in Jupyter notebooks

## Installation

Clone this repository and install in editable mode:

```bash
git clone https://github.com/andreacicconardi94/Research-Work.git
cd Research-Work
pip install -e .
