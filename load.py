import numpy as np
from scipy.constants import m_e, c, e
from ElegantWrapper.watcher import Watcher2

m_e_eV = m_e*c**2/e

def astra_dist_to_watcher(dist_file, allowed_status=(-1,5), remove_first=True):
    dist = np.loadtxt(dist_file)
    columns_dict = {}
    columns_dict['x'] = dist[:,0]
    columns_dict['y'] = dist[:,1]
    columns_dict['z'] = -dist[:,2]/c
    momentum_z = dist[:,5]
    momentum_z[1:] += dist[0,5]
    columns_dict['p'] = momentum_z / m_e_eV
    pis0 = momentum_z == 0
    columns_dict['xp'] = columns_dict['p'].copy()
    columns_dict['yp'] = columns_dict['p'].copy()
    columns_dict['xp'][pis0] = 0
    columns_dict['xp'][~pis0] = dist[:,3][~pis0]/momentum_z[~pis0]
    columns_dict['yp'][pis0] = 0
    columns_dict['yp'][~pis0] = dist[:,4][~pis0]/momentum_z[~pis0]
    columns_dict['t'] = columns_dict['clock'] = dist[:,6]*1e-9
    columns_dict['status'] = dist[:,9]
    if remove_first:
        for key, val in columns_dict.items():
            columns_dict[key] = val[1:]
    status = columns_dict['status']
    if allowed_status is not None:
        mask = np.zeros_like(status, bool)
        for allowed in allowed_status:
            mask = np.logical_or(mask, status==allowed)
        for key, val in columns_dict.items():
            columns_dict[key] = val[mask]
    return Watcher2({}, columns_dict)

def load_astra_temit(emit_file):
    data = np.loadtxt(emit_file)
    out = {
            'z': data[:,0],
            't': data[:,1]*1e-9,
            'pos_avr': data[:,2]*1e-3,
            'pos_rms': data[:,3]*1e-3,
            'ang_rms': data[:,4]*1e-3,
            'nemit': data[:,5]*1e-6,
            'posang_avr': data[:,6]*1e-6, # unsure
            }
    return out

def load_astra_zemit(emit_file):
    data = np.loadtxt(emit_file)
    out = {
            'z': data[:,0],
            't': data[:,1]*1e-9,
            'E': data[:,2]*1e6,
            'z_rms': data[:,3]*1e-3,
            'E_rms': data[:,4]*1e3,
            'nemit': data[:,5],
            'zE_avr': data[:,6], # what is this?
            }
    return out

