# sem_patch_extractor.py
"""SEM Patch Extractor
=====================
Extract square **N µm × N µm** patches from SEM `.tif` images.

Supports two modes:
1. **Metadata mode**: use embedded Image Pixel Size to compute patch size.
2. **Scalebar mode**: input a scalebar image to measure 1 µm in pixels.

Features
--------
* `detect_scalebar_length` → returns pixel length of 1 µm from scalebar image.
* `extract_patches_by_scalebar` → splits full SEM image into patches of given µm side.
* `extract_patches_from_tif` → original metadata‐based extractor.
* Auto‐crop raw images to remove footer/logo.

Example (scalebar mode)
-----------------------
>>> from sem_patch_extractor import extract_patches_by_scalebar
>>> patches = extract_patches_by_scalebar(
...     scalebar_image="scalebar.tif",
...     full_image="sample.tif",
...     save_dir="patches/",
...     patch_um=2,
...     overlap=0
... )
>>> len(patches)

Example (metadata mode)
-----------------------
>>> from sem_patch_extractor import extract_patches_from_tif
>>> patches = extract_patches_from_tif(
...     file_path="sample.tif",
...     save_dir="patches/",
...     auto_crop=True,
...     patch_um=1.5
... )
>>> len(patches)
"""
from __future__ import annotations
import re
from typing import List, Optional

# ---------------------------------------------------------------------------
# Scalebar‐based functions
# ---------------------------------------------------------------------------

# Improved version of detect_scalebar_length_user_guided
# Supports scalebars with vertical end caps (e.g., [——])

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from PIL import Image
from scipy import ndimage
from pathlib import Path


def detect_scalebar_length_user_guided(image_path: str, roi_size=(300, 60)) -> int:
    """
    Allows user to click near the scalebar and detects it from that region.

    Parameters:
        image_path (str): Path to the SEM image containing a scalebar.
        roi_size (tuple): Size (width, height) of the region around the click to search the scalebar.

    Returns:
        int: Detected scalebar length in pixels (1 micron), or raises an error if detection fails.
    """
    img = np.array(Image.open(image_path).convert("L"))
    h, w = img.shape

    # Request user click
    fig, ax = plt.subplots()
    ax.imshow(img, cmap="gray")
    ax.set_title("Click near the scalebar")
    plt.axis("off")
    coords = plt.ginput(1, timeout=0)
    plt.close(fig)

    if not coords:
        raise ValueError("No click received. Aborting scalebar detection.")

    cx, cy = map(int, coords[0])
    rw, rh = roi_size

    # Extract ROI
    x0 = max(cx - rw // 2, 0)
    x1 = min(cx + rw // 2, w)
    y0 = max(cy - rh // 2, 0)
    y1 = min(cy + rh // 2, h)
    roi = img[y0:y1, x0:x1]

    # Threshold ROI
    thr = roi.mean() + 0.5 * (roi.max() - roi.mean())
    mask = roi > thr
    lbl, n = ndimage.label(mask)

    best_rect = None
    best_width = 0
    min_distance = float('inf')

    for lab in range(1, n + 1):
        coords = np.column_stack(np.where(lbl == lab))
        if coords.size == 0:
            continue
        y_min, x_min = coords.min(axis=0)
        y_max, x_max = coords.max(axis=0)
        width = x_max - x_min + 1
        height = y_max - y_min + 1
        cx_r = (x_min + x_max) / 2
        cy_r = (y_min + y_max) / 2
        dist = (cx_r - rw/2)**2 + (cy_r - rh/2)**2
        if dist < min_distance and width > roi.shape[1] * 0.3:
            best_width = width
            best_rect = (x_min, y_min, width, height)
            min_distance = dist

    if not best_rect:
        raise ValueError("Scalebar could not be detected in the selected region.")

    # Plot ROI with detected region
    fig1, ax1 = plt.subplots()
    ax1.imshow(roi, cmap="gray")
    rect1 = mpatches.Rectangle((best_rect[0], best_rect[1]), best_rect[2], best_rect[3],
                              linewidth=2, edgecolor="red", facecolor="none")
    ax1.add_patch(rect1)
    ax1.set_title(f"Scalebar in ROI: {best_width} px")
    ax1.axis("off")

    # Plot full image with selection location
    fig2, ax2 = plt.subplots()
    ax2.imshow(img, cmap="gray")
    rect2 = mpatches.Rectangle((x0 + best_rect[0], y0 + best_rect[1]), best_rect[2], best_rect[3],
                               linewidth=2, edgecolor="cyan", facecolor="none")
    ax2.add_patch(rect2)
    ax2.set_title("Scalebar position in full image")
    ax2.axis("off")

    plt.show()

    return best_width


def extract_patches_by_scalebar(
    scalebar_image: str,
    full_image: str,
    save_dir: str,
    patch_um: int = 1,
    overlap: int = 0,
) -> list[str]:
    """
    Split SEM `full_image` into non-overlapping patches of side `patch_um` microns,
    using a user-guided scalebar detection from `scalebar_image`.

    Parameters:
        scalebar_image (str): Image file with visible scalebar.
        full_image (str): Image to divide into patches.
        save_dir (str): Output directory to save patches.
        patch_um (int): Side of the patch in microns.
        overlap (int): Overlap between patches in pixels.

    Returns:
        list[str]: List of saved patch file paths.
    """
    L = detect_scalebar_length_user_guided(scalebar_image)
    if L <= 0:
        raise ValueError("Failed to detect scalebar length.")

    patch_px = patch_um * L
    img = np.array(Image.open(full_image).convert("L"))
    h, w = img.shape
    out_dir = Path(save_dir)

    # Ensure output directory is clean
    if out_dir.exists():
        for file in out_dir.glob("*.png"):
            file.unlink()
    out_dir.mkdir(parents=True, exist_ok=True)

    paths: list[str] = []
    step = patch_px - overlap
    idx = 0
    for y in range(0, h - patch_px + 1, step):
        for x in range(0, w - patch_px + 1, step):
            tile = img[y : y + patch_px, x : x + patch_px]
            out_path = out_dir / f"patch_{idx:04d}.png"
            Image.fromarray(tile).save(out_path)
            paths.append(str(out_path))
            idx += 1
    return paths



# ---------------------------------------------------------------------------
# Original metadata‐based extractor
# ---------------------------------------------------------------------------

def _read_sem_image(
    file_path: str,
) -> tuple[np.ndarray, dict[str, Optional[float | str]]]:
    """Read SEM .tif and return (image_array, metadata dict)."""
    from PIL import TiffTags, Image as PILImage
    tags: list[tuple[int, str, object]] = []
    with PILImage.open(file_path) as img:
        if hasattr(img, "tag_v2"):
            for tag_id, value in img.tag_v2.items():
                name = TiffTags.TAGS.get(tag_id, f"Unknown Tag ({tag_id})")
                tags.append((tag_id, name, value))
        arr = np.array(img)
    meta: dict[str, Optional[float | str]] = {"pixel_size": None, "pixel_unit": None}
    try:
        long_block = tags[6][2]  # type: ignore[index]
        for line in long_block.split("\r\n"):
            if "Image Pixel Size" in line:
                nums = re.findall(r"\d+\.\d+|\d+", line)
                if nums:
                    val = float(nums[0])
                    unit = "nm" if "nm" in line.lower() else "µm"
                    meta["pixel_size"] = val
                    meta["pixel_unit"] = unit
                break
    except Exception:
        pass
    return arr, meta


def read_sem_metadata(
    file_path: str,
) -> dict[str, Optional[float | str]]:
    """Return metadata dict containing 'pixel_size' and 'pixel_unit'."""
    _, meta = _read_sem_image(file_path)
    return meta.copy()


def _px_per_um(
    pixel_size: Optional[float],
    unit: Optional[str],
) -> Optional[float]:
    """Convert pixel size in nm/µm to px per µm."""
    if pixel_size is None or unit is None:
        return None
    unit = unit.lower()
    if unit == "nm":
        return 1000.0 / pixel_size
    if unit in {"µm", "um"}:
        return 1.0 / pixel_size
    return None


def _auto_crop(
    img: np.ndarray,
    thresh_ratio: float = 0.15,
) -> np.ndarray:
    """Remove dark margins via adaptive threshold + largest component."""
    imf = img.astype(np.float32)
    thr = imf.min() + thresh_ratio * (imf.max() - imf.min())
    mask = imf > thr
    lbl, n = ndimage.label(mask)
    if n == 0:
        return img
    counts = ndimage.sum(mask, lbl, index=range(1, n + 1))
    main = int(np.argmax(counts) + 1)
    fg = lbl == main
    rows = np.flatnonzero(fg.any(axis=1))
    cols = np.flatnonzero(fg.any(axis=0))
    r0, r1 = rows[[0, -1]]
    c0, c1 = cols[[0, -1]]
    return img[r0 : r1 + 1, c0 : c1 + 1]


def _extract_patches(
    arr: np.ndarray,
    patch_px: int,
    overlap: int = 0,
    min_fg_ratio: float = 0.05,
) -> list[np.ndarray]:
    """Slide window and keep patches with ≥ min_fg_ratio fg pixels."""
    patches: list[np.ndarray] = []
    step = patch_px - overlap
    h, w = arr.shape
    for y in range(0, h - patch_px + 1, step):
        for x in range(0, w - patch_px + 1, step):
            tile = arr[y : y + patch_px, x : x + patch_px]
            if (tile > tile.mean()).mean() < min_fg_ratio:
                continue
            patches.append(tile)
    return patches


def extract_patches_from_tif(
    file_path: str,
    save_dir: str,
    *,
    auto_crop: bool = True,
    overlap: int = 0,
    viz_crop: bool = False,
    patch_um: float = 1.0,
) -> list[str]:
    """Create N µm × N µm patches from SEM image via metadata."""
    img, meta = _read_sem_image(file_path)
    img_use = _auto_crop(img) if auto_crop else img
    px_per_um = _px_per_um(meta.get("pixel_size"), meta.get("pixel_unit"))
    if px_per_um is None:
        raise ValueError("Missing Image Pixel Size metadata — cannot compute patch size.")
    patch_px = max(1, round(patch_um * px_per_um))
    patches = _extract_patches(img_use, patch_px, overlap)
    out_dir = Path(save_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = Path(file_path).stem
    saved: list[str] = []
    for i, tile in enumerate(patches):
        p = out_dir / f"{stem}_patch_{i:04d}.png"
        Image.fromarray(tile).save(p)
        saved.append(str(p))
    if viz_crop:
        _save_crop_preview(img, img_use, out_dir / f"{stem}_crop_preview.png")
    return saved


def _save_crop_preview(
    full_img: np.ndarray,
    crop: np.ndarray,
    out_png: Path,
) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    rgb = Image.fromarray(full_img).convert("RGB")
    h0, w0 = full_img.shape
    h1, w1 = crop.shape
    off_y = (h0 - h1) // 2
    off_x = (w0 - w1) // 2
    fig, ax = plt.subplots()
    ax.imshow(rgb, cmap="gray")
    rect = patches.Rectangle((off_x, off_y), w1, h1,
                             edgecolor="red", facecolor="none",
                             linewidth=2)
    ax.add_patch(rect)
    ax.axis("off")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
