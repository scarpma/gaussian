from db_utils import *

run = 19
number = 1750

path = f'/scratch/scarpolini/gaussian/wgangp/runs/{run}/{number}_gen.h5'
print('Loading Model ...')
gen = load_model(path)


N = 10
bs = 50000
trajs = np.zeros(shape=(N*bs,2000,1))
print('Generating Trajectories ...')
for ii in range(N):
    print(ii)
    noise = np.random.normal(0, 1, size=(bs, 100))
    trajs[ii*bs:(ii+1)*bs,:,0:1] = gen.predict(noise, verbose=1, batch_size=bs)


try: os.mkdir(f'/storage/scarpolini/databases/gaussian/wgangp/runs/{run}')
except: print('Directory already exists')
print('Saving ...')
np.save(f'/storage/scarpolini/databases/gaussian/wgangp/runs/{run}/gen_trajs_{number}', trajs)
