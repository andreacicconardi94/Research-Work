# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 15:54:50 2023

@author: acicconardi
"""

import os
import statistics
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from tensorflow.keras.models import load_model
from skimage import measure
from sklearn.metrics import jaccard_score


def Test(desired_height, desired_width, model, folder):
    # Load the input in PNG format
    input_list = [file for file in os.listdir(folder) if 'mask_' not in file]
    mask_list = [file for file in os.listdir(folder) if 'mask_' in file]
    
    Hyperparameters = [i for i in np.arange(0.0,1.0,0.1)]
    # Variance = [[], [], [], []]
    Threshold = [[], [], [], []]
    Jacc_sim = []
    
    for idx, file in enumerate(input_list):
        test_image = Image.open(folder+"\\"+file)
        mask_image = Image.open(folder+'\\'+'mask_'+file)
        # Convert the image in a greyscale, if needed
        test_image = test_image.convert("L")
        mask_image = mask_image.convert("L")
        # Preprocess image as in training
        test_array = np.array(test_image) / 255.0
        mask_array = np.array(mask_image) / 255.0
        
        pad_height = max(desired_height - test_array.shape[0], 0)
        pad_width = max(desired_width - test_array.shape[1], 0)
    
        # Apply zero-padding
        padded_test = np.pad(test_array, ((0, pad_height), (0, pad_width)), mode='constant', constant_values=0.)
        padded_mask = np.pad(mask_array, ((0, pad_height), (0, pad_width)), mode='constant', constant_values=0.)
    
        # Add a batch dimension
        padded_test = np.expand_dims(padded_test, axis=0)
        padded_mask = np.expand_dims(padded_mask, axis=0)
        
        # Execute the prediction
        prediction = model.predict(padded_test)
        

        threshold = 0.75
        binary_mask = (prediction >= threshold).astype(np.uint8)

        # Visualize the input and the binary mask
        plt.imshow(binary_mask[0, :, :, 0], cmap='gray')
        plt.title("Binary Mask "+ file)
        plt.show()


        # Calculate Jaccard similarity (Intersection over Union)
        jaccard_similarity = jaccard_score(padded_mask.flatten(), binary_mask[0, :, :, 0].flatten())

        Threshold.append(threshold)
        Jacc_sim.append(jaccard_similarity)
        
        print("Accuracy (Jaccard Similarity):", jaccard_similarity)
    
    
    return Jacc_sim, Hyperparameters, Threshold
    
    
def count_nanocrystal(mask):
    # Reveal the border of nanocrystals
    border = measure.find_contours(mask, 0.5)
    
    # Apply label algorithm to connected components
    labels = measure.label(mask)
    
    # Count the nanocrystals
    num_nanocristalli = np.max(labels)
    
    return num_nanocristalli, border
    
    
    
if __name__=='__main__':
    folder = r"C:\Users\acicconardi\Documents\Corsi di dottorato\Machine Learning\Progetto CNN\Dataset\Test_5"

    desired_height = 512  # Desired height for padding
    desired_width = 320   # Desired width for padding
    # Loading the model
    model = load_model(r'C:\Users\acicconardi\Documents\Corsi di dottorato\Machine Learning\Progetto CNN\Dataset\U_net_5.h5')
    
    JS, Iper, Thr = Test(desired_height, desired_width, model, folder)


    
    plt.figure(figsize=(10, 6))
    x_labels = ['Test 1', 'Test 2', 'Test 3', 'Test 4']
    x = np.arange(len(x_labels))
    plt.xticks(x, x_labels)
    plt.xlabel('Test')
    plt.ylabel('Jaccard Similarity')
    plt.yticks(np.arange(0,1,0.05))
    plt.ylim(0,1)
    bar_width = 0.4
    plt.bar(x, JS, color='blue', width=bar_width, alpha=0.7)
    plt.legend()
    plt.title('Results of Jaccard Similarity for each Test')
    plt.grid(True)
    plt.show()

    

    
    
    
    
    
    
    
    