import os
from collections import OrderedDict
from pprint import pprint

import numpy as np
from scipy.interpolate import splrep, splev
import sncosmo


class SampledFunction(object):
    def __init__(self, x, y):
        self.x = np.asarray(x, dtype=np.float64)
        self.y = np.asarray(y, dtype=np.float64)
        self.xmin = x[0]
        self.xmax = x[-1]
        self._tck = splrep(self.x, self.y, k=1)

    def __call__(self, x):
        f = splev(x, self._tck, ext=1)


class AggregateBandpass(object):
    """Bandpass with spatially uniform transmission defined by multiple
    transmissions in series.

    Parameters
    ----------
    transmissions : list of SampledFunctions
    """
    def __init__(self, transmissions):
        self.transmissions

    def __repr__(self):
        return ("AggregateBandpass(" + repr(self.transmissions) + ")")

    def __call__(self, wave):
        t = self.transmissions[0](wave)
        for trans in self.transmissions[1:]:
            t *= trans(wave)
        return t


class RadialBandpass(object):
    """Bandpass where transmission depends on radius.

    Parameters
    ----------
    transmission : SampledFunction or list thereof
    dependent_transmissions : list of SampledFunction
    radii : list of float
    """
    def __init__(self, transmission, radial_transmissions, radii):
        self.transmissions = ([transmission]
                              if type(transmission) == SampledFunction
                              else transmission)
        self.radial_transmissions = radial_transmissions
        self.radii = radii


def simplify_sampled_functions(sfs):
    """Given a list of SampledFunctions and scalars (representing a constant
    sampled function), return a list of just SampledFunctions, with scalars
    multiplied into the functions."""
    pass

def _parse_value(s):
    try:
        x = int(s)
    except:
        try:
            x = float(s)
        except:
            x = s
    return x


def read_cards(fname):
    cards = OrderedDict()
    with open(fname, 'r') as f:
        for line in f:
            if line[0] != '@':
                continue

            words = line.split()
            key = words[0][1:]  # words is at least length 1
            tokens = words[1:]
            if len(tokens) == 0:
                value = None
            elif len(tokens) == 1:
                value = _parse_value(tokens[0])
            else:
                value = [_parse_value(v) for v in tokens]

            cards[key] = value

    return cards


def read_filterwheel(fname):
    
    # read filter filenames
    bandfnames = {}
    with open(fname, 'r') as f:
        for line in f:
            words = line.split()  # each line can have 2 or 3 words.
            band = words[0]
            bandfname = words[-1]
            bandfnames[band] = bandfname

    transmissions = {}
    for band, bandfname in bandfnames.items():
        x, y = np.loadtxt(bandfname, unpack=True)
        transmissions[band] = SampledFunction(x, y)

    return transmissions


def read_radial_filterwheel(fname):
    """Read radially variable filterwheel transmissions.

    Parameters
    ----------
    fname : str
        Name of filterwheel file, which contains a list of files defining the
        transmission for each bandpass.

    Returns
    -------
    dict of list
        Dictionary where keys are filter names and values are lists. Each
        item in the list is a two-tuple giving the radius and filter function
        (in the form of a SampledFunction).
    """
    
    # read filter filenames (multiple files per filter)
    bandfnames = {}
    with open(fname, 'r') as f:
        for line in f:
            band, _, bandfname = line.split()
            if band not in bandfnames:
                bandfnames[band] = []
            bandfnames[band].append(bandfnames)

    transmissions = {band: [] for band in bandfnames}
    for band, bandfnames in bandfnames.items():
        for bandfname in bandfnames:
            # read the filter function at a single radius.
            # TODO: re-organize the salt2-format reader.
            meta, data = sncosmo.io._read_salt2(bandfname)
            rad = meta["MEASUREMENT_RADIUS"].split()[0]  # in cm
            sf = SampledFunction(data['lambda'], data['tr'])
            transmissions[band].append((rad, sf))

        # sort transmissions (just in case they weren't in order)
        transmissions[band].sort(key=lambda x: x[0])

    return transmissions

    
def read_snfit_instrument(dirname, cardfile='instrument.cards'):
    """Read a set of bandpasses in the snfit format"""

    cards = read_cards(os.path.join(dirname, cardfile))

    transmissions = []  # scalars or `SampledFunction`s

    # required keys:
    for key in ("MIRROR_REFLECTIVITY",
                "OPTICS_TRANS",
                "QE",
                "ATMOSPHERIC_TRANS"):
        value = cards[key]
        if type(value) is int or type(value) is float:
            transmissions.append(value)
        else:
            fname = os.path.join(dirname, value)
            x, y = np.loadtxt(fname, unpack=True)
            transmissions.append(SampledFunction(x, y))

    # There can be *either* a `FILTERS` or a `RADIALLY_VARIABLE_FILTERS`
    # keyword.
    if "FILTERS" in cards:
        fname = os.path.join(dirname, cards["FILTERS"])
        filter_transmissions = read_filterwheel(fname)

    elif "RADIALLY_VARIABLE_FILTERS" in cards:  # -> filterwheel file
        fname = os.path.join(dirname, cards["RADIALLY_VARIABLE_FILTERS"])
        filter_transmissions = read_radial_filterwheel(fname)

    else:
        raise ValueError("'FILTERS' or 'RADIALLY_VARIABLE_FILTERS' keyword "
                         "not found")

    if "CHROMATIC_CORRECTIONS" in cards:
        fname = os.path.join(dirname, cards["CHROMATIC_CORRECTIONS"])
        corr_transmissions = read_filterwheel(fname)

        # check that keys are in filters
        for key in corr_transmissions:
            if key not in filter_transmissions:
                raise ValueError("band {!r} in {!r} not in filterwheel"
                                 .format(key, fname))
    else:
        corr_transmissions = 
        
    # construct filters
    if "FILTERS" in cards:
        pass

    else:  # radially variable filters
        

    print(transmissions)
        # how to translate to filter names?
    # A: Return dictionary of filters in filterwheel, translate externally

    # how does Product() function work in detail? (can we aggregate
    # everything into one 1-d function upon instantiation rather than
    # reflectivity, optics, etc?)

    # need a SampledFunction class to represent the above. See how Product()
    # works.

    # how will this be represented such that one can do:
    # sncosmo.get_bandpass("MEGACAMPSF::g", x, y)

    # A1: get_bandpass checks type of band... if it is variable, demand
    #     arguments
    # A2: maybe all bandpasses should have e.g., .at() attributes?
    #     get_bandpass always calls it, sometimes with no arguments.
    # A3: maybe need a new function get_transmission() or similar, used in
    #     integration functions, that returns transmission on a grid with
    #     set spacing. This would still call Bandpass.at() with 0 or more args.

# testing
read_snfit_instrument('snfit_data/Instruments/Megacam-PSF')
