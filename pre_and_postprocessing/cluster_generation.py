import matplotlib.pyplot as plt
import pandas as pd
from sklearn.datasets import make_blobs

# Generate sample cluster points
n_samples = 80000  # Total number of points
n_features = 2   # Number of features (dimensions)
n_clusters = 8   # Number of clusters

X, y = make_blobs(n_samples=n_samples, centers=n_clusters, n_features=n_features, random_state=0x1337, center_box=(-20.0, 20.0))

#print("y = ", y)

#print("X[:,0], X[:,1] = ", X)

df = pd.DataFrame(
    {'Feature 1': X[:,0],
     'Feature 2': X[:,1],
     'Cluster': y
    })

df.to_csv('cluster_points_80k.csv', index=False)

# Plot the generated clusters
plt.scatter(X[:, 0], X[:, 1], c=y, cmap='viridis', alpha=0.6)
plt.title("Generated Cluster Points")
plt.xlabel("Feature 1")
plt.ylabel("Feature 2")
plt.show()