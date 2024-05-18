import numpy as np
from scipy.constants import c
from ElegantWrapper.watcher import Watcher2

def astra_dist_to_watcher(dist_file, allowed_status=(-1,5), remove_first=True):
    dist = np.loadtxt(dist_file)
    columns_dict = {}
    columns_dict['x'] = dist[:,0]
    columns_dict['y'] = dist[:,1]
    columns_dict['t'] = -dist[:,2]/c
    columns_dict['p'] = -dist[:,5]
    columns_dict['xp'] = dist[:,3]/columns_dict['p']
    columns_dict['yp'] = dist[:,4]/columns_dict['p']
    columns_dict['clock'] = dist[:,6]*1e-9
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



