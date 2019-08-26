# Constants
PI = 3.1415926535897932384626433832
RHO_WATER = 1e3  # kg/m3
RHO_AIR = 1.0
THIRD = 1./3.
MULTI_THRESH = 1e4


# Global droplet counter
superdroplet_count = 1


class Superdroplet(object):
    """ Simple wrapper for a 'Superdroplet' object, to encapsulate
    functionality  """

    def __init__(self, multi, rcubed, solute):
        global superdroplet_count
        superdroplet_count += 1

        self.multi = multi
        self.rcubed = rcubed
        self.solute = solute
        self.density = RHO_WATER  # kg/m3
        self.id = superdroplet_count

    def set_properties(self, multi, rcubed, solute):
        """ Mutator - modifies state without returning a copy """
        self.multi = multi
        self.rcubed = rcubed
        self.solute = solute

    @property
    def terminal_v(self):
        return self._calc_terminal_v()

    def _calc_terminal_v(self):
        # BEARD, 1976
        r = self.rcubed**(1./3.)
        d = 2.*r*1e6  # diameter, m -> micron
        x = self.calc_mass() * 1e3  # convert kg -> g

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

        return 1e-2 * alpha * x_to_beta  # from cm/s -> m/s

    @property
    def volume(self):
        return self.rcubed*4.*PI/3.

    @property
    def mass(self):
        return self.density*self.calc_volume()

    def __repr__(self):
        return "%d" % self.id
