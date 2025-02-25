import cv2
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from PIL import Image
import json


def load_image(image_path):
    """Loads the image and converts it into a NumPy array."""
    image = Image.open(image_path)
    image = image.convert("RGB")  # Ensure it is in RGB format
    return np.array(image)

def resize_image(image, target_shape):
    """Resizes the image to match the specified target shape (height, width)."""
    return cv2.resize(image, (target_shape[1], target_shape[0]), interpolation=cv2.INTER_NEAREST)

def cluster_image(image, n_clusters):
    """Segments the image into n_clusters using K-Means clustering."""
    pixels = image.reshape(-1, 3)
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = kmeans.fit_predict(pixels)
    segmented_image = labels.reshape(image.shape[:2])
    return segmented_image, kmeans.cluster_centers_

def create_cluster_masks(clustered_image, n_clusters):
    """Creates binary masks for each cluster.
    
    Each mask is a 2D array with value 1 for pixels belonging to the cluster and 0 elsewhere.
    """
    masks = {}
    for cluster in range(n_clusters):
        # (clustered_image == cluster) returns a boolean array; convert it to uint8 (0 and 1)
        mask = (clustered_image == cluster).astype(np.uint8)
        masks[cluster] = mask
    return masks

def save_cluster_masks(masks, base_filename="cluster_mask"):
    """Saves each cluster mask as a separate .npy file."""
    for cluster, mask in masks.items():
        filename = f"{base_filename}_{cluster}.npy"
        np.save(filename, mask)
        print(f"Saved mask for cluster {cluster} to {filename}")

def visualize_mask(mask):
    """Visualizes a binary mask."""
    plt.imshow(mask, cmap='gray')
    plt.title("Cluster Mask")
    plt.show()


if __name__=='__main__':
    image_path = "mappa.png"  # Change to your file path
    image = load_image(image_path)
    
    # Display the original image
    plt.imshow(image)
    plt.title("Original Image")
    plt.show()
    
    # Set target shape (height, width) as extracted from the hypermap (e.g., 387x1324)
    height = int(input("Please insert the height to rescale the picture"))
    width = int(input("Please insert the width to rescale the picture"))
    target_shape = (height, width)
    resized_image = resize_image(image, target_shape)
    
    # Define number of clusters manually
    n_clusters = int(input("Please insert the number of clusters in the map"))
    
    # Perform clustering
    clustered_image, cluster_centers = cluster_image(resized_image, n_clusters)
    
    # Display the segmented image (clusters)
    plt.imshow(clustered_image, cmap='tab10')
    plt.title("Segmented Map")
    plt.colorbar()
    plt.show()
    
    # Create binary masks for each cluster
    masks = create_cluster_masks(clustered_image, n_clusters)
    
    # Save the masks as .npy files (each file contains a 2D array with 0's and 1's)
    save_cluster_masks(masks)
    
    # Optionally, visualize one mask (e.g. for cluster 2)
    visualize_mask(masks[2])
