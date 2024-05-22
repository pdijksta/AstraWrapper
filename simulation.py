import os
import re
import glob
import numpy as np

from .load import load_astra_temit, load_astra_zemit

class AstraSimulation:
    def __init__(self, infile, distfilename, charge):
        self.distfilename = distfilename
        self.charge = charge

        self.infile = infile
        self.infile_base = os.path.basename(infile)
        self.dirname = os.path.dirname(infile)
        self.basename = self.infile_base.split('.in')[0]

        self.emit_x = load_astra_temit(os.path.join(self.dirname, '%s.Xemit.001' % self.basename))
        self.emit_y = load_astra_temit(os.path.join(self.dirname, '%s.Yemit.001' % self.basename))
        self.emit_z = load_astra_zemit(os.path.join(self.dirname, '%s.Zemit.001' % self.basename))

        dist_files0 = glob.glob(self.dirname+'/%s.*.001' % self.basename)
        re_dist = re.compile('%s\.(\d{4})\.001' % self.basename)
        self.dist_files = sorted(filter(lambda x: re_dist.match(os.path.basename(x)), dist_files0))
        self.dist_files_s = np.array([0] + [int(re_dist.match(os.path.basename(x)).group(1))/100 for x in self.dist_files])
        self.dist_files.insert(0, distfilename)

