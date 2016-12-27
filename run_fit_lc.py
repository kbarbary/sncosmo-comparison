#!/usr/bin/env python

import os
import glob
import re
import time
from pprint import pprint

import sncosmo
import snfitio  # to read snfit results

REGEX = re.compile('jla_light_curves/lc-(.+)\.list')


def snname_from_fname(fname):
    return re.match(REGEX, fname).groups()[0]


if __name__ == "__main__":

    model = sncosmo.Model(source='salt2',
                          effects=[sncosmo.F99Dust()],
                          effect_names=['mw'],
                          effect_frames=['obs'])

    fnames = glob.glob("jla_light_curves/lc-SDSS19230.list")
    for fname in fnames[0:1]:
        snname = snname_from_fname(fname)

        data = sncosmo.read_lc(fname, format='salt2', read_covmat=True)

        model.set(mwebv=data.meta['MWEBV'], z=data.meta['Z_HELIO'])
        t0 = time.time()
        result, m = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'],
                                   modelcov=True, phase_range=(-15., 45.),
                                   wave_range=(3000., 7000.), verbose=True)
        print("time:", time.time() - t0, 's')
        print(result)

        snfit_result = snfitio.read_snfit_result('results_snfit/result-{}.dat'
                                                 .format(snname))
        print(snfit_result)
