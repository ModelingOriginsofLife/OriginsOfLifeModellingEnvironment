import numpy as np
from sklearn.decomposition import KernelPCA

data = np.loadtxt("heredity_PCA/rawdata.txt")

kpca = KernelPCA(kernel="rbf")
X_kpca = kpca.fit_transform(data)

np.savetxt("heredity_PCA/kpcascores.txt", X_kpca)
np.savetxt("heredity_PCA/kpcaeigs.txt", kpca.lambdas_)