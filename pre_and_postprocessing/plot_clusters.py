import matplotlib.pyplot as plt
import pandas as pd
import os

folder = os.getcwd()
file = folder + "/output_data.txt"
df = pd.read_csv(file, header=None, sep=r"\s+", skiprows=1)

# Plot the generated clusters
plt.scatter(df[0], df[1], c=df[2], cmap='viridis', alpha=0.6)
plt.title("Classified Cluster Points")
plt.xlabel("Feature 1")
plt.ylabel("Feature 2")
plt.show()