import numpy as np
import matplotlib.pyplot as plt

def random_dataset(d, num_clusters, avg_cluster_size, center_scale):
    while True:
        sizes = np.random.normal(avg_cluster_size, 100, num_clusters).round()
        if np.alltrue(sizes > 2):
            break
    centers = center_scale*np.random.normal(0, 5, num_clusters*d).reshape((-1,d))
    clusters = []
    for i in range(num_clusters):
        pts = np.random.normal(0, 1, int(sizes[i])*d).reshape((-1,d))
        pts += centers[i]
        clusters.append(pts)
    points = np.vstack(clusters)
    return points

def write_fvecs(points, fname):
    with open(fname, "wb") as f:
        for v in points:
            f.write(np.int32(len(v)).tobytes())
            f.write(v.tobytes())

points = random_dataset(8,30,20000,2)
plt.scatter(points[:,0], points[:,1], s=1)
plt.show()

points = points.astype('float32')
write_fvecs(points, "points")
