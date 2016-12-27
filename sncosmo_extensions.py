from collections import OrderedDict
import copy
import os
from pprint import pprint

import numpy as np
from scipy.interpolate import splrep, splev
import sncosmo
import astropy.units as u


class SampledFunction(object):
    def __init__(self, x, y):
        self.x = np.asarray(x, dtype=np.float64)
        self.y = np.asarray(y, dtype=np.float64)
        self.xmin = x[0]
        self.xmax = x[-1]
        self._tck = splrep(self.x, self.y, k=1)

    def __call__(self, x):
        return splev(x, self._tck, ext=1)
        

class AggregateBandpass(object):
    """Bandpass with spatially uniform transmission defined by multiple
    transmissions in series.

    Parameters
    ----------
    transmissions : list of SampledFunctions
    """
    def __init__(self, transmissions, prefactor=1.0):
        if len(transmissions) < 1:
            raise ValueError("empty list of transmissions")
        self.transmissions = transmissions
        self.prefactor = prefactor
        
    def __repr__(self):
        return ("AggregateBandpass(" + repr(self.transmissions) +
                ", prefactor=" + repr(self.prefactor) + ")")

    def __call__(self, wave):
        t = self.transmissions[0](wave)
        for trans in self.transmissions[1:]:
            t *= trans(wave)
        t *= self.prefactor
        return t


class RadialBandpassSet(object):
    """A set of Bandpasses defined at different radii.

    Parameters
    ----------
    transmissions : list of SampledFunction
    radial_transmissions :  list of (float, SampledFunction)
    """
    def __init__(self, transmissions, radial_transmissions, prefactor=1.0):
        self.prefactor = prefactor
        self.trans = transmissions
        self.rtrans = radial_transmissions

    def at(self, r):
        """Return the bandpass at the given radius"""

        if r < 0.0:
            raise ValueError("negative radius")

        # interpolate radial transmission: find index
        i = 1
        while i < len(self.rtrans) and r > self.rtrans[i][0]:
            i += 1
        if i == len(self.rtrans):
            raise ValueError("radius greater than maximum radius of {:f}"
                             .format(self.rtrans[-1][0]))

        # linearly interpolate second transmission onto first
        weight = (r - self.rtrans[i-1][0]) / (self.rtrans[i][0] - self.rtrans[i-1][0])
        x = self.rtrans[i-1][1].x
        y = weight * self.rtrans[i-1][1].y + (1.0 - weight) * self.rtrans[i][1](x)

        trans = copy.copy(self.trans)
        trans.append(SampledFunction(x, y))

        return AggregateBandpass(trans, prefactor=self.prefactor)


def separate_sampled_functions(sfs):
    """Given a list of SampledFunctions and scalars (representing a constant
    sampled function), collect scalars into a single prefactor. Return the
    prefactor and a list of the SampledFunctions."""

    true_sfs = []
    prefactor = 1.0
    for sf in sfs:
        if type(sf) is int or type(sf) is float:
            prefactor *= sf
        else:
            true_sfs.append(sf)

    return prefactor, true_sfs
            

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


def read_filterwheel(dirname, fname, skiprows=0):
    
    # read filter filenames
    bandfnames = {}
    with open(os.path.join(dirname, fname), 'r') as f:
        for line in f:
            words = line.split()  # each line can have 2 or 3 words.
            band = words[0]
            bandfname = words[-1]
            bandfnames[band] = bandfname

    transmissions = {}
    for band, bandfname in bandfnames.items():
        x, y = np.loadtxt(os.path.join(dirname, bandfname), unpack=True,
                          skiprows=skiprows)
        transmissions[band] = SampledFunction(x, y)

    return transmissions


def read_radial_filterwheel(dirname, fname):
    """Read radially variable filterwheel transmissions.

    Parameters
    ----------
    fname : str
        Basename of filterwheel file, which contains a list of files defining the
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
    with open(os.path.join(dirname, fname), 'r') as f:
        for line in f:
            band, _, bandfname = line.split()
            if band not in bandfnames:
                bandfnames[band] = []
            bandfnames[band].append(bandfname)

    transmissions = {band: [] for band in bandfnames}
    for band, bandfname_list in bandfnames.items():
        for bandfname in bandfname_list:
            # read the filter function at a single radius.
            # TODO: re-organize the salt2-format reader.
            with open(os.path.join(dirname, bandfname), 'r') as f:
                meta, data = sncosmo.io._read_salt2(f)
            try:
                rad_str = meta["MEASUREMENT_RADIUS"]
            except KeyError:
                raise Exception("MEASUREMENT_RADIUS keyword not found in " +
                                os.path.join(dirname, bandfname))
            
            rad = float(rad_str.split()[0]) # parse string like '0 cm'
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
            x, y = np.loadtxt(os.path.join(dirname, value), unpack=True)
            transmissions.append(SampledFunction(x, y))

    # optional key:
    if "CHROMATIC_CORRECTIONS" in cards:
        corr_transmissions = read_filterwheel(dirname,
                                              cards["CHROMATIC_CORRECTIONS"],
                                              skiprows=3)
    else:
        corr_transmissions = None
            
    # There can be *either* a `FILTERS` or a `RADIALLY_VARIABLE_FILTERS`
    # keyword but not both.
    if "FILTERS" in cards:
        filter_transmissions = read_filterwheel(dirname, cards["FILTERS"])

        bands = {}
        for name in filter_transmissions:
            # get a list of all transmissions for this band
            trans_list = copy.copy(transmissions)
            trans_list.append(filter_transmissions[name])
            if corr_transmissions is not None:
                trans_list.append(corr_transmissions.get(name, 1.0))

            # simplify the list
            prefactor, trans_list = separate_sampled_functions(trans_list)

            # if there's only one sampled function, we can construct a normal
            # bandpass
            if len(trans_list) == 1:
                w = trans_list[0].x
                t = trans_list[0].y
                bands[name] = sncosmo.Bandpass(w, prefactor * t,
                                               wave_unit=u.AA)
            elif len(trans_list) > 1:
                bands[name] = AggregateBandpass(trans_list, prefactor=prefactor)
            else:
                raise Exception("Filter {} consists only of scalars!"
                                .format(name))

    elif "RADIALLY_VARIABLE_FILTERS" in cards:
        filter_transmissions = read_radial_filterwheel(
            dirname, cards["RADIALLY_VARIABLE_FILTERS"])

        bands = {}

        for name in filter_transmissions:
            # get a list of all transmissions for this band
            trans_list = copy.copy(transmissions)
            if corr_transmissions is not None:
                trans_list.append(corr_transmissions.get(name, 1.0))

            # simplify the list
            prefactor, trans_list = separate_sampled_functions(trans_list)

            bands[name] = RadialBandpass(trans_list, filter_transmissions[name],
                                         prefactor=prefactor)

    else:
        raise ValueError("'FILTERS' or 'RADIALLY_VARIABLE_FILTERS' keyword "
                         "not found")

    return bands



    # how will this be represented such that one can do:
    # sncosmo.get_bandpass("MEGACAMPSF::g", x, y)

    # A1: get_bandpass checks type of band... if it is variable, demand
    #     arguments
    # A2: maybe all bandpasses should have e.g., .at() attributes?
    #     get_bandpass always calls it, sometimes with no arguments.
    # A3: maybe need a new function get_transmission() or similar, used in
    #     integration functions, that returns transmission on a grid with
    #     set spacing. This would still call Bandpass.at() with 0 or more args.


_BANDSETS = None


def get_bandset(name):
    global _BANDSETS
    if _BANDSETS is None:
        bandsets = read_snfit_instrument('snfit_data/Instruments/Megacam-PSF')
        _BANDSETS = {'MEGACAMPSF::' + key: value
                     for key, value in bandsets.items()}
    return _BANDSETS[name]


def expand_snls_filters(data):
    """Return a new data table with snls megacam filters evaluted at the 
    radius given in the header"""

    r = math.sqrt(data.meta['X_FOCAL_PLANE']**2 +
                  data.meta['Y_FOCAL_PLANE']**2)

    bands = np.empty(len(data), dtype=np.object)
    
    for name in set(data['Filter']):
        mask = data['Filter'] == name
        
        # only treat things we know are 'bandsets'
        # TODO: replace this with something more general: lookup in dictionary
        # of known bandsets?
        if name.startswith("MEGACAMPSF"):
            bands[mask] = get_bandset(name).at(r)
        else:
            bands[mask] = sncosmo.get_bandpass(name)
        
    data.replace_column('Filter', bands)
