# sem_patch_extractor.py
"""SEM Patch Extractor
======================
Extract square **1 µm × 1 µm** patches from single‑frame SEM ``.tif`` images
acquired on Zeiss GEMINI or Thermo‑Fisher Helios microscopes.

The script can work with **raw micrographs** that still contain the black footer
/ logo *or* with images that have already been cropped – set
``auto_crop=False`` in that case.

Features
--------
* Reads the *Image Pixel Size* tag → converts to **px / µm** automatically.
* Optional **auto‑crop** removes dark footer / logo margins.
* **Patch filter** discards tiles that are mostly background (< 5 % "bright"
  pixels).
* **Preview**: saves a quick PNG with a red rectangle showing the crop.
* Helper ``read_sem_metadata`` lets you inspect pixel size & unit easily.

Quick example
--------------
>>> from sem_patch_extractor import extract_patches_from_tif, read_sem_metadata
>>> meta = read_sem_metadata(r"C:/SEM/AIX_19.tif")
>>> print(meta["pixel_size"], meta["pixel_unit"])
8.0 nm
>>> extract_patches_from_tif(
...     file_path=r"C:/SEM/AIX_19_cropped.tif",   # already cropped
...     save_dir=r"C:/SEM/patches",
...     auto_crop=False,
...     overlap=0,
...     viz_crop=True,
... )
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np
from PIL import Image
from scipy import ndimage  # connected‑component labelling

# ---------------------------------------------------------------------------
# Low‑level helpers
# ---------------------------------------------------------------------------

def _read_sem_image(file_path: str) -> Tuple[np.ndarray, dict]:
    """Return *(image_array, metadata)* for a GEMINI / Helios SEM .tif file."""
    from PIL import TiffTags, Image as PILImage

    tags: list[tuple[int, str, object]] = []
    with PILImage.open(file_path) as img:
        if hasattr(img, "tag_v2"):
            for tag_id, value in img.tag_v2.items():
                tag_name = TiffTags.TAGS.get(tag_id, f"Unknown Tag ({tag_id})")
                tags.append((tag_id, tag_name, value))
                
                
        arr = np.array(img)

    meta: dict[str, Optional[float | str]] = {
        "pixel_size": None,
        "pixel_unit": None,
    }


    pixel_size_candidates = []
    try:
        long_block: str = tags[6][2]  # type: ignore[index]
        meta["pixel_size_raw"] = None  # 'Pixel Size'
        for line in long_block.split("\r\n"):
            if "Pixel Size =" in line:
                val = float(re.findall(r"\d+\.\d+|\d+", line)[0])
                unit = "nm" if "nm" in line.lower() else "µm"
                if val > 0:
                    pixel_size_candidates.append((val, unit))

    # pick the largest valid Pixel Size value (e.g. Detector pixel size)
        if pixel_size_candidates:
            best_val, best_unit = max(pixel_size_candidates, key=lambda x: x[0])
            meta["pixel_size_raw"] = f"{best_val} {best_unit}"

    except Exception:
        pass  # leave as None if parsing fails

    return arr, meta


def read_sem_metadata(file_path: str) -> dict[str, Optional[float | str]]:
    """Public helper: return the metadata dict (pixel size + unit, etc.)."""
    _, meta = _read_sem_image(file_path)
    return meta.copy()

# ---------------------------------------------------------------------------
# Calibration + auto‑crop
# ---------------------------------------------------------------------------

def _px_per_um(pixel_size: float | None, unit: str | None) -> Optional[float]:
    """Convert *pixel_size* (given in nm or µm) into pixels / µm."""
    if pixel_size is None or unit is None:
        return None
    unit = unit.lower()
    if unit == "nm":
        return 1000.0 / pixel_size  # 1 µm = 1000 nm
    if unit in {"µm", "um"}:
        return 1.0 / pixel_size
    return None


def _auto_crop(img: np.ndarray, thresh_ratio: float = 0.15) -> np.ndarray:
    """Remove dark footer/logo via adaptive threshold + largest component."""
    img_f = img.astype(np.float32)
    thr = img_f.min() + thresh_ratio * (img_f.max() - img_f.min())
    mask = img_f > thr

    labels, nlab = ndimage.label(mask)
    if nlab == 0:
        return img  # no foreground detected → return full image

    counts = ndimage.sum(mask, labels, index=range(1, nlab + 1))
    main_lab = int(np.argmax(counts) + 1)
    fg = labels == main_lab

    rows = np.flatnonzero(fg.any(axis=1))
    cols = np.flatnonzero(fg.any(axis=0))
    r0, r1 = rows[[0, -1]]
    c0, c1 = cols[[0, -1]]
    return img[r0 : r1 + 1, c0 : c1 + 1]

# ---------------------------------------------------------------------------
# Patch extraction
# ---------------------------------------------------------------------------

def _extract_patches(
    arr: np.ndarray,
    patch_px: int,
    overlap: int = 0,
    min_fg_ratio: float = 0.05,
) -> List[np.ndarray]:
    """Slide a window and keep patches with ≥ *min_fg_ratio* foreground pixels."""
    patches: list[np.ndarray] = []
    step = patch_px - overlap
    h, w = arr.shape
    for y in range(0, h - patch_px + 1, step):
        for x in range(0, w - patch_px + 1, step):
            tile = arr[y : y + patch_px, x : x + patch_px]
            if (tile > tile.mean()).mean() < min_fg_ratio:
                continue  # skip almost‑black tiles
            patches.append(tile)
    return patches

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def extract_patches_from_tif(
    file_path: str,
    save_dir: str,
    *,
    auto_crop: bool = True,
    overlap: int = 0,
    viz_crop: bool = False,
) -> List[str]:
    """Create 1 µm × 1 µm patches from a SEM image and save them as PNG.

    Parameters
    ----------
    file_path : str
        Path to the input ``.tif`` file.
    save_dir : str
        Directory where patch PNGs will be written.
    auto_crop : bool, optional (default **True**)
        If *False* the image is used as‑is (assumed already trimmed).
    overlap : int, optional (default 0)
        Overlap between adjacent patches, in **pixels**.
    viz_crop : bool, optional (default **False**)
        Save an extra PNG showing the crop rectangle (quick QA).

    Returns
    -------
    list[str]
        Absolute paths of the saved patch images.
    """
    img, meta = _read_sem_image(file_path)

    img_use = _auto_crop(img) if auto_crop else img

    px_per_um = _px_per_um(meta.get("pixel_size"), meta.get("pixel_unit"))
    if px_per_um is None:
        raise ValueError("Missing Image Pixel Size metadata – cannot determine patch size.")
    patch_px = max(1, round(px_per_um))

    patches = _extract_patches(img_use, patch_px, overlap)

    save_dir_path = Path(save_dir)
    save_dir_path.mkdir(parents=True, exist_ok=True)
    stem = Path(file_path).stem

    saved_paths: list[str] = []
    for idx, tile in enumerate(patches):
        out_path = save_dir_path / f"{stem}_patch_{idx:04d}.png"
        Image.fromarray(tile).save(out_path)
        saved_paths.append(str(out_path))

    if viz_crop:
        _save_crop_preview(img, img_use, save_dir_path / f"{stem}_crop_preview.png")

    return saved_paths

# ---------------------------------------------------------------------------
# Preview helper
# ---------------------------------------------------------------------------

def _save_crop_preview(full_img: np.ndarray,
                       crop: np.ndarray,
                       out_png: Path) -> None:
    """Save a PNG with a red rectangle indicating the retained region."""
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    rgb = Image.fromarray(full_img).convert("RGB")

    # crude centre alignment (good enough for quick visual QA)
    h0, w0 = full_img.shape
    h1, w1 = crop.shape
    offset_y = (h0 - h1) // 2
    offset_x = (w0 - w1) // 2

    fig, ax = plt.subplots()
    ax.imshow(rgb, cmap="gray")
    rect = patches.Rectangle(
        (offset_x, offset_y),
        w1,
        h1,
        linewidth=2,
        edgecolor="red",
        facecolor="none",
    )
    ax.add_patch(rect)
    ax.axis("off")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)