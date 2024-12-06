# -*- coding: utf-8 -*-
"""
@author: acicconardi
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt
from scipy.optimize import curve_fit
import pandas as pd


def get_folder_path():
    """Prompt the user to input the folder path containing the data files."""
    while True:
        folder = input("Please enter the folder path containing the data files: ").strip()
        if os.path.isdir(folder):
            print(f"Folder '{folder}' found.")
            return folder
        else:
            print(f"Invalid folder path: '{folder}'. Please try again.")


def get_bounds_from_user(num_curves):
    """Prompt the user to input bounds for Lorentzian parameters."""
    bounds = {"lower": [], "initial": [], "upper": []}
    for i in range(num_curves):
        print(f"\nEnter bounds for Lorentzian peak {i + 1}:")

        for param in ["Amplitude", "Center", "Width"]:
            while True:
                try:
                    lower = float(input(f"  Lower bound for {param}: "))
                    initial = float(input(f"  Initial value for {param}: "))
                    upper = float(input(f"  Upper bound for {param}: "))

                    if lower <= initial <= upper:
                        bounds["lower"].append(lower)
                        bounds["initial"].append(initial)
                        bounds["upper"].append(upper)
                        break
                    else:
                        print("  Invalid input. Ensure that: lower <= initial <= upper.")
                except ValueError:
                    print("  Invalid input. Please enter numeric values.")
    return bounds


def single_lorentzian(x, amplitude, center, width):
    """Calculate a single Lorentzian peak."""
    return amplitude * width**2 / ((x - center) ** 2 + width**2)


def multiple_lorentzian(x, *params):
    """Sum of multiple Lorentzian peaks."""
    num_params = len(params) // 3
    result = np.zeros_like(x)
    for i in range(num_params):
        amplitude, center, width = params[i * 3:(i + 1) * 3]
        result += single_lorentzian(x, amplitude, center, width)
    return result


def fit_spectrum(wavelengths, intensity, num_curves, use_bounds, bounds, label, title):
    """
    Fit the spectrum with a specified number of Lorentzian peaks.
    
    Args:
        wavelengths (np.array): Array of wavelengths.
        intensity (np.array): Array of intensity values.
        num_curves (int): Number of Lorentzian peaks to fit.
        use_bounds (bool): Whether to use bounds in the fitting.
        bounds (dict): Dictionary containing bounds for the fitting.
        label (str): Label for the plot.
        title (str): Title for the plot.
        
    Returns:
        dict: Fitting results including peak centers and ratios.
    """
    # Normalize intensity
    intensity = intensity - np.min(intensity)

    # Prepare fitting parameters
    if use_bounds:
        p0 = bounds["initial"]
        lower_bounds = bounds["lower"]
        upper_bounds = bounds["upper"]
        bounds_tuple = (lower_bounds, upper_bounds)
    else:
        p0 = bounds["initial"]
        bounds_tuple = (-np.inf, np.inf)  # No bounds

    # Fit the data
    try:
        popt, _ = curve_fit(
            multiple_lorentzian, wavelengths, intensity,
            p0=p0,
            bounds=bounds_tuple
        )
    except RuntimeError:
        print(f"Fit failed for {label}.")
        return None

    # Generate fitted curve and individual Lorentzians
    fitted_curve = multiple_lorentzian(wavelengths, *popt)
    individual_peaks = [
        single_lorentzian(wavelengths, *popt[i * 3:(i + 1) * 3])
        for i in range(num_curves)
    ]

    # Plot the results
    plt.plot(wavelengths, intensity, color="blue", label="PL Spectrum")
    for i, peak in enumerate(individual_peaks, start=1):
        plt.plot(wavelengths, peak, label=f"Peak {i}")
    plt.plot(wavelengths, fitted_curve, color="black", linestyle="--", label="Total Fit")
    plt.title(title)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Counts (a.u.)")
    plt.legend(loc="upper right")
    plt.show()

    # Collect results: Peak centers and ratios (if applicable)
    results = {"Peak Centers": [popt[i * 3 + 1] for i in range(num_curves)]}
    if num_curves > 1:
        peak_ratios = {f"Ratio Peak {j+1}/{j}": max(individual_peaks[j]) / max(individual_peaks[j - 1])
                       for j in range(1, num_curves)}
        results.update(peak_ratios)

    return results


if __name__ == "__main__":
    # Prompt the user for folder path
    folder_path = get_folder_path()

    # Retrieve all .txt files from the folder
    file_list = [file for file in os.listdir(folder_path) if file.endswith(".txt")]

    # Prompt user for the number of Lorentzian peaks
    while True:
        try:
            num_curves = int(input("Enter the number of Lorentzian peaks to fit (1, 2, or 3): "))
            if num_curves in [1, 2, 3]:
                break
            else:
                print("Invalid choice. Please enter 1, 2, or 3.")
        except ValueError:
            print("Invalid input. Please enter a numeric value.")

    # Prompt user for bounds usage
    while True:
        use_bounds_input = input("Do you want to use bounds for the fitting? (yes/no): ").strip().lower()
        if use_bounds_input in ["yes", "no"]:
            use_bounds = use_bounds_input == "yes"
            break
        else:
            print("Invalid input. Please enter 'yes' or 'no'.")

    # If bounds are used, collect them from the user
    if use_bounds:
        bounds = get_bounds_from_user(num_curves)
    else:
        # If no bounds are used, initialize with default values
        bounds = {"initial": [1.0, 600.0, 10.0] * num_curves}

    # Dictionary to store results
    results_dict = {}

    for file in file_list:
        file_path = os.path.join(folder_path, file)
        data = np.loadtxt(file_path, delimiter="\t", skiprows=1)

        wavelengths = data[:, 0]
        intensity = data[:, 1]

        # Apply median filter to smooth the intensity
        intensity_filtered = medfilt(intensity, kernel_size=5)

        # Extract file info for labeling
        label = file.split("_")[0]
        title = "_".join(file.split("_")[:4])

        # Perform fitting
        results = fit_spectrum(wavelengths, intensity_filtered, num_curves, use_bounds, bounds, label, title)
        if results is not None:
            results_dict[label] = results

    # Save results to CSV
    result_df = pd.DataFrame.from_dict(results_dict, orient="index")
    result_df.to_csv(os.path.join(folder_path, "lorentzian_fit_results.csv"))

    print("Fitting completed. Results saved.")

