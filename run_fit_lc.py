#!/usr/bin/env python

import os
import glob
import re
from pprint import pprint

import sncosmo

REGEX = re.compile('jla_light_curves/lc-(.+)\.list')


def snname_from_fname(fname):
    return re.match(REGEX, fname).groups()[0]


if __name__ == "__main__":

    fnames = glob.glob("jla_light_curves/lc-*.list")
    for fname in fnames[0:1]:
        snname = snname_from_fname(fname)

        print(snname)

        data = sncosmo.read_lc(fname, format='salt2')
        pprint(data.meta)
        print(data)

    
