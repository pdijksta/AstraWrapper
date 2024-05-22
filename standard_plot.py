import numpy as np
from scipy.constants import m_e, e, c

from PassiveWFMeasurement import myplotstyle as ms
from .load import astra_dist_to_watcher

m_e_eV = m_e*c**2/e

def plot(sim):
    infile = sim.infile
    charge = sim.charge

    fig = ms.figure('Standard plot for %s' % infile)
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    subplot = ms.subplot_factory(3, 3)
    sp_ctr = 1

    sp_energy = subplot(sp_ctr, title='Energy', xlabel='$s$ (m)', ylabel='$E$ (MeV)')
    sp_ctr += 1
    sp_espread = sp_energy.twinx()
    sp_espread.set_ylabel('$\sigma_E$ (keV)')

    sp_energy.plot(sim.emit_z['z'], sim.emit_z['E']/1e6, label='$E$')
    sp_espread.plot(sim.emit_z['z'], sim.emit_z['E_rms']/1e3, label='$\sigma_E$', color='tab:orange')
    ms.comb_legend(sp_energy, sp_espread)

    sp_beta = subplot(sp_ctr, title='Beta functions', xlabel='$s$ (m)', ylabel=r'$\beta$ (m)')
    sp_ctr += 1

    sp_current = subplot(sp_ctr, title='Current profiles', xlabel='$t$ (ps)', ylabel='$I$ (A)')
    sp_ctr += 1

    for dim, emit_dict in [
            ('X', sim.emit_x),
            ('Y', sim.emit_y),
            ]:
        emit = emit_dict['nemit']/(sim.emit_z['E']/m_e_eV)
        sp_beta.plot(emit_dict['z'], emit_dict['pos_rms']**2/emit, label=dim)

    sp_beta.legend()

    for n_dist, (dist_file, s) in enumerate(zip(sim.dist_files, sim.dist_files_s)):
        w = astra_dist_to_watcher(dist_file, allowed_status=(-1, 5))
        zz = w['clock']
        #zz_mean, zz_std = zz0.mean(), zz0.std()
        #min_z = zz_mean - n_sig_z*zz_std
        #max_z = zz_mean + n_sig_z*zz_std
        #mask = np.logical_and(zz0>min_z, zz0<max_z)
        #print('Dist file %i, consider %i/%.0f percent of particles for plot' % (n_dist, mask.size, mask.sum()/mask.size*100))
        #print(np.unique(w['status'], return_counts=True))
        #zz = zz0[mask]
        hist, curr_zz = np.histogram(zz-zz.mean(), bins=100)
        factor = charge/(np.diff(curr_zz)[0]*np.sum(hist))
        curr = np.concatenate(([0,],hist*factor))
        sp_current.step(curr_zz*1e12, curr, label='%.1f' % s)

        sp = subplot(sp_ctr, grid=False, title='Dist at s=%.1f m, %i particles' % (s, zz.size), xlabel='$t$ (ps)', ylabel='$\Delta E$ (MeV)')
        sp_ctr += 1


        hist2d, yedges, xedges = np.histogram2d(w['p'], w['clock'], bins=(100, 100))
        x_axis = (xedges[1:]+xedges[:-1])/2
        y_axis = (yedges[1:]+yedges[:-1])/2
        x_factor = 1e12
        y_factor = 1e-6
        extent = [x_axis[0]*x_factor, x_axis[-1]*x_factor, y_axis[0]*y_factor, y_axis[-1]*y_factor]
        sp.imshow(hist2d, aspect='auto', extent=extent, origin='lower', cmap=ms.plt.get_cmap('hot'))

    sp_current.legend()

