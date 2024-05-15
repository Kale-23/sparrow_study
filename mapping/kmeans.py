import numpy as np 
from sklearn.cluster import KMeans 

"""
############################################################################
"""
# randomly initializing K centroid by picking K samples from X
def initialize_random_centroids(K, X, useRand=True):
    if (useRand ) : 
        """Initializes and returns k random centroids"""
        m, n = np.shape(X)
        # a centroid should be of shape (1, n), so the centroids array will be of shape (K, n)
        centroids = np.empty((K, n))
        for i in range(K):
            # pick a random data point from X as the centroid
            centroids[i] =  X[np.random.choice(range(m))] 
        return centroids
    else :
        """ Assumes K=6  and  D=4 """
        centroids = np.array([[-0.18948127,  0.16649768,  0.02406444, -0.00601086],
          [ 0.53000236,  0.38294265, -0.15929458,  0.02184399],
          [ 0.1877262,  -0.22771385,  0.33630654,  0.08129628],
          [ 0.06614255, -0.18489067,  0.12179142,  0.06047589],
          [-0.16177389,  0.0932761,  -0.03489587,  0.08000099]]  )
        return centroids 

def euclidean_distance(x1, x2):
    """Calculates and returns the euclidean distance between two vectors x1 and x2"""
    return np.sqrt(np.sum(np.power(x1 - x2, 2)))

def closest_centroid(x, centroids, K):
    """Finds and returns the index of the closest centroid for a given vector x"""
    distances = np.empty(K)
    for i in range(K):
        distances[i] = euclidean_distance(centroids[i], x)
    return np.argmin(distances) # return the index of the lowest distance

def create_clusters(centroids, K, X):
    """Returns an array of cluster indices for all the data samples"""
    m, _ = np.shape(X)
    cluster_idx = np.empty(m)
    for i in range(m):
        cluster_idx[i] = closest_centroid(X[i], centroids, K)
    return cluster_idx

def compute_means(cluster_idx, K, X):
    """Computes and returns the new centroids of the clusters"""
    _, n = np.shape(X)
    centroids = np.empty((K, n))
    for i in range(K):
        points = X[cluster_idx == i] # gather points for the cluster i
        centroids[i] = np.mean(points, axis=0) # use axis=0 to compute means across points
    return centroids

def run_Kmeans(K, X, max_iterations=1900, useRand=True ):
    """Runs the K-means algorithm and computes the final clusters"""
    # initialize random centroids
    centroids = initialize_random_centroids(K, X, useRand=useRand)
    # loop till max_iterations or convergance
    print(f"initial centroids: {centroids}")
    kmeansMach = KMeans( init=centroids, n_clusters=K, random_state=0, n_init=1, max_iter=1900 ) 
    kmeansMach.fit( X ) 
    return(  kmeansMach.labels_ ) 

    for _ in range(max_iterations):
        # create clusters by assigning the samples to the closet centroids
        clusters = create_clusters(centroids, K, X)
        previous_centroids = centroids                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
        # compute means of the clusters and assign to centroids
        centroids = compute_means(clusters, K, X)
        # if the new_centroids are the same as the old centroids, return clusters
        diff = previous_centroids - centroids
        if not diff.any():
            return clusters
    return clusters

"""
############################################################################
"""

