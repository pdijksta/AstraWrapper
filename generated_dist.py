import numpy as np

def center_clock(filename):
    dist = np.loadtxt(filename)
    mask = dist[:,-1] == -1
    mask[0] = 0
    dist[1:,6] -= dist[mask,6].mean()
    np.savetxt(filename, dist)
    print('Centered t coordinate for %s' % filename)

