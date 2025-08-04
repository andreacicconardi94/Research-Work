# SEM Patch Extractor

A lightweight Python tool to extract **N µm × N µm** square patches from SEM `.tif` images.

Supports two modes:
- 🖱️ **User-guided scalebar mode**: detect the pixel length of 1 µm by clicking near a visible scalebar.
- 🧠 **Metadata mode**: read pixel size from embedded SEM TIFF metadata.

---

## ✨ Features

- 🔍 **Interactive scalebar selection**: click near the scalebar and confirm detection.
- 📏 Patch extraction based on **physical size** (in µm).
- 🧪 **Dual visual feedback**: see detected scalebar in ROI and full image.
- 🧼 Automatic cropping (removes dark margins/logos).
- 🔁 Overlap and patch filtering options.
- ✅ Designed for **Jupyter notebooks**.
- 🔧 Modular code for integration and extension.

---

## 📦 Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/andreacicconardi94/Research-Work.git
cd Research-Work
pip install -e .
```

---

## 🧑‍💻 Usage

### 1. **Extract patches (user-guided scalebar detection)**

```python
from sem_patch_extractor import extract_patches_by_scalebar

patches = extract_patches_by_scalebar(
    scalebar_image="images/AIX_03_scalebar.tif",
    full_image="images/PlainImages/AIX_03.tif",
    save_dir="patches_output/AIX_03",
    patch_um=3,      # patch size in microns
    overlap=0        # optional pixel overlap
)
```

You will:
- Click near the scalebar.
- See a preview of the detected bar in the ROI and in the full image.
- Get patches saved automatically (any existing ones will be overwritten).

---

### 2. **Alternative: extract from TIFF metadata**

```python
from sem_patch_extractor import extract_patches_from_tif

patches = extract_patches_from_tif(
    image_path="images/SEM_AIX_10.tif",
    save_dir="patches_output/AIX_10",
    patch_um=2,
    overlap=0,
    crop_margins=True,
    viz_crop=True
)
```

Requires that TIFF metadata include `"Image Pixel Size"` in µm.

---

## 📝 Notes

- All patches are saved as `.png` in the specified output folder.
- Patch size is computed in pixels based on the detected 1 µm scale.
- Any `.png` files already in the output folder will be **overwritten**.

---

## 📂 Project structure

```
sem_patch_extractor/
├── sem_patch_extractor/
│   ├── __init__.py
│   ├── patches.py
├── patch_extractor_example.ipynb
├── README.md
├── pyproject.toml
```

