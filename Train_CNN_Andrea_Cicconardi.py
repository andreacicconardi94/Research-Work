# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 14:08:09 2023

@author: acicconardi
"""

import os
from PIL import Image
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers
import matplotlib.pyplot as plt
from tensorflow.keras.layers import Dropout
from tensorflow.keras import regularizers
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score



def unet(input_shape):
    inputs = tf.keras.Input(shape=input_shape)

    # Encoding
    conv1 = layers.Conv2D(64, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(inputs)
    conv1 = layers.Conv2D(64, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv1)
    pool1 = layers.MaxPooling2D(pool_size=(2, 2))(conv1)

    conv2 = layers.Conv2D(128, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(pool1)
    conv2 = layers.Conv2D(128, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv2)
    pool2 = layers.MaxPooling2D(pool_size=(2, 2))(conv2)
    
    conv3 = layers.Conv2D(256, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(pool2)
    conv3 = layers.Conv2D(256, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv3)
    # drop4 = Dropout(0.5)(conv3)
    pool3 = layers.MaxPooling2D(pool_size=(2, 2))(conv3)
    
    conv4 = layers.Conv2D(512, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(pool3)
    conv4 = layers.Conv2D(512, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv4)
    drop4 = Dropout(0.5)(conv4)
    pool4 = layers.MaxPooling2D(pool_size=(2, 2))(drop4)


    # Decoding
    conv5 = layers.Conv2D(256, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(pool4)
    conv5 = layers.Conv2D(256, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv5)
    up1 = layers.UpSampling2D(size=(2, 2))(conv5)
    concat1 = layers.Concatenate()([up1, conv4])
    conv6 = layers.Conv2D(512, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(concat1)
    conv6 = layers.Conv2D(512, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv6)

    conv7 = layers.Conv2D(256, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(pool3)
    conv7 = layers.Conv2D(256, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv7)
    up2 = layers.UpSampling2D(size=(2, 2))(conv7)
    concat2 = layers.Concatenate()([up2, conv3])
    conv8 = layers.Conv2D(128, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(concat2)
    conv8 = layers.Conv2D(128, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv8)
    
    conv9 = layers.Conv2D(128, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv8)
    conv9 = layers.Conv2D(128, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv9)
    up3 = layers.UpSampling2D(size=(2, 2))(conv9)
    concat3 = layers.Concatenate()([up3, conv2])
    conv10 = layers.Conv2D(64, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(concat3)
    conv10 = layers.Conv2D(64, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv10)
    
    conv11 = layers.Conv2D(64, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv10)
    conv11 = layers.Conv2D(64, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv11)
    up4 = layers.UpSampling2D(size=(2, 2))(conv11)
    concat4 = layers.Concatenate()([up4, conv1])
    conv12 = layers.Conv2D(32, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(concat4)
    conv12 = layers.Conv2D(32, 3, activation="relu", padding="same", kernel_regularizer=regularizers.l2(0.01))(conv12)
    

    decoded = layers.Conv2D(1, 1, activation="sigmoid")(conv12)



    # Create the model
    model = tf.keras.Model(inputs=inputs, outputs=decoded)
    return model

# Generate the dataset for training
def Dataset(desired_height, desired_width):
    folder = r'C:\Users\acicconardi\Documents\Corsi di dottorato\Machine Learning\Progetto CNN\Dataset\Train_5'
    input_list = [file for file in os.listdir(folder) if 'mask_' not in file]
    dataset_list = []
    for idx, file in enumerate(input_list):
        input_image = Image.open(folder+"\\"+file)
        mask_image = Image.open(folder+"\mask_"+file)
        # Convert the image in a greyscale
        input_image = input_image.convert("L")
        mask_image = mask_image.convert("L")
        # Normalize the image
        input_array = np.array(input_image) / 255.0
        mask_array = np.array(mask_image) / 255
        # Padding dimension
        pad_height = max(desired_height - input_array.shape[0], 0)
        pad_width = max(desired_width - input_array.shape[1], 0)
        padded_input = np.pad(input_array, ((0, pad_height), (0, pad_width)), mode='constant', constant_values=0.)
        padded_mask = np.pad(mask_array, ((0, pad_height), (0, pad_width)), mode='constant', constant_values=0)
        padded_input = np.expand_dims(padded_input, axis=0)
        padded_mask = np.expand_dims(padded_mask, axis=0)
        dataset_1 = tf.data.Dataset.from_tensor_slices((padded_input, padded_mask))
        dataset_list.append(dataset_1)
        
    dataset = dataset_list[0]
    for i in dataset_list[1:]:
        dataset = dataset.concatenate(i)
        

    # Shuffle e batch of dataset
    batch_size = 1
    dataset = dataset.shuffle(buffer_size=1000).batch(batch_size)

    return dataset, padded_input


if __name__=='__main__':
    desired_height = 512  # Desired height of image
    desired_width = 320   # Desired width of image
    
    dataset, padded_input = Dataset(desired_height, desired_width)
    
    # Define input dimensions
    input_shape = (padded_input.shape[1], padded_input.shape[2], 1)
    
    # Create U-Net instance
    model = unet(input_shape)
    
    # Compile model
    learning_rate = 0.001  # Set the learning rate
    model.compile(optimizer=Adam(learning_rate), loss="binary_crossentropy")
    
    history = model.fit(dataset, epochs=10, batch_size=1)
        
    # Get the loss during epochs
    loss = history.history['loss']
    
    # Plot the loss function
    plt.plot(loss)
    plt.title('Loss durante le epoche')
    plt.xlabel('Epoca')
    plt.ylabel('Perdita')
    plt.show()
    
    model.save(r'C:\Users\acicconardi\Documents\Corsi di dottorato\Machine Learning\Progetto CNN\Dataset\U_net_5.h5')
   

    







