
from numba import jit

from utils import RHO_WATER, PI

superdroplet_count = 0

class Superdroplet(object):

    def __init__(self, multi, rcubed, solute):
        global superdroplet_count
        superdroplet_count += 1

        self.multi  = multi
        self.rcubed = rcubed
        self.solute = solute
        self.density = RHO_WATER # kg/m3
        self.id = superdroplet_count

    def set_properties(self, multi, rcubed, solute):
        self.multi = multi
        self.rcubed = rcubed
        self.solute = solute

    @property
    def terminal_v(self):
        return self.calc_terminal_v(self.rcubed, self.mass)

    # @jit(nopython=True)
    @staticmethod
    def calc_terminal_v(rcubed, mass):
        ## BEARD, 1976
        r = rcubed**(1./3.)
        d = 2.*r*1e6 # diameter, m -> micron
        x = mass * 1e3 # convert kg -> g

        if d <= 134.43:
            alpha = 4.5795e5
            x_to_beta = x**(2./3.)
        elif 134.43 < d <= 1511.64:
            alpha = 4962.0
            x_to_beta = x**(1./3.)
        elif 1511.64 < d <= 3477.84:
            alpha = 1732.0
            x_to_beta = x**(1./6.)
        else:
            alpha = 917.0
            x_to_beta = 1.0

        return 1e-2 * alpha * x_to_beta # from cm/s -> m/s

    @property
    def volume(self):
        return self.calc_volume(self.rcubed)

    # @jit(nopython=True)
    @staticmethod
    def calc_volume(rcubed):
        return rcubed*4.*PI/3.

    @property
    def mass(self):
        return self.calc_mass(self.density, self.volume)

    # @jit(nopython=True)
    @staticmethod
    def calc_mass(density, volume):
        return density*volume

    def __repr__(self):
        return "%d" % self.id
