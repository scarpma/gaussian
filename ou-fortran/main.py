import stats
import numpy as np


print('loading')
db = np.load('../../databases/velocities.npy')
print('velocities loaded')

print('python shape', db.shape)
sf = stats.compute_sf_mixed(db, npart=5000)

np.savetxt('prova.dat', sf)
