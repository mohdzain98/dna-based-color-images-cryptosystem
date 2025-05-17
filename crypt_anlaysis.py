import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from scipy.stats import chisquare
import math
from numba import njit


@njit
def fast_variance(image):
    mean_val = np.mean(image)
    diff_squared = (image - mean_val) ** 2
    var = np.sum(diff_squared) / (image.shape[0] * image.shape[1])
    return round(var, 3)


class experimentalAnalysis:
    def __init__(self, image, directory):
        self.m = image.shape[0]
        self.n = image.shape[1]
        self.directory = directory
        os.makedirs(self.directory, exist_ok=True)  # Ensure directory exists
        print("m,n:", self.m, self.n)

    def plotHistogram(self, red, green, blue, bins, name="histogram"):
        fig = plt.figure(figsize=(4, 3))
        values = red.flatten()
        plt.hist(
            values, bins=bins, color="red", alpha=0.6
        )  # Adjust the number of bins as needed
        values = green.flatten()
        plt.hist(values, bins=bins, color="green", alpha=0.6)
        values = blue.flatten()
        plt.hist(values, bins=bins, color="blue", alpha=0.6)
        plt.xlabel("Value")
        plt.ylabel("Frequency")
        title = (
            "Histogram of the Original Image"
            if name == "org_img_histogram"
            else "Histogram of the Cipher Image"
        )
        plt.title(title)
        file_path = os.path.join(self.directory, f"{name}.png")
        plt.savefig(file_path, dpi=300, bbox_inches="tight")
        print(f"Histogram saved successfully at {file_path}")
        plt.show()

    def calMean(self, image):
        intensity = np.mean(image)
        return round(intensity, 3)

    def variance(self, image):
        m = self.m
        n = self.n
        var = 0
        for i in range(0, m):
            for j in range(0, n):
                var = var + (pow((image[i][j] - self.calMean(image)), 2))
        var = var / (m * n)
        return round(var, 3)

    def variance_numba(self, image):
        return fast_variance(image)

    def chiSquare(self, image):
        hist, _ = np.histogram(image.flatten(), bins=256, range=(0, 256))
        expected_distribution = np.full(256, len(image) * len(image[0]) / 256)

        # Calculate the chi-squared statistic
        chi_squared, p_value = chisquare(f_obs=hist, f_exp=expected_distribution)
        return chi_squared

    def entropy(self, image):
        values = image.flatten()

        # Calculate the histogram of values
        hist, bins = np.histogram(
            values, bins=256
        )  # You can adjust the number of bins as needed

        # Convert the histogram to a probability distribution
        probabilities = hist / len(values)

        # Calculate the entropy
        entropy = -np.sum(probabilities * np.log2(probabilities + np.finfo(float).eps))
        return entropy

    def corAnalysisH(self, image):
        ans = 0
        for i in range(40):
            ans = ans + np.corrcoef(image[i], image[i + 1])[0][1]
        ans = ans / (40)
        return ans

    def corAnalysisV(self, image):
        ans = 0
        for j in range(40):
            ans = ans + np.corrcoef(image[:, j], image[:, j + 1])[0][1]
        ans = ans / 40
        return ans

    def corAnalysisD(self, matrix):
        ans = 0
        diagOne = []
        diagTwo = []
        for i in range(1, 50):
            for j in range(200):
                diagOne.append(matrix[j][i + j])
                diagTwo.append(matrix[j + i][j])
            ans = ans + np.corrcoef(diagOne, diagTwo)[0][1]
        ans = ans / 30
        return ans

    def noiceAttack(self, image):
        for i in range(0, 100):
            for j in range(0, 100):
                image[i][j] = 0
        return image

    def calculate_PSNR(self, original_image, cipher_image):
        MSE = 0
        m = n = 256
        for i in range(0, m):
            row = 0
            for j in range(0, n):
                row = row + ((original_image[i][j] - cipher_image[i][j]) ** 2)
            MSE = MSE + row
        MSE = MSE / (m * n)
        # Calculate PSNR in decibels
        PSNR = 10 * math.log10((255 * 255) / MSE)
        return PSNR
