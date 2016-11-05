"""Parse output from snfit"""
from collections import OrderedDict

import numpy as np
import six
from sncosmo.utils import Result


def _parse_numeric_value(s):
    try:
        x = int(s)
    except:
        x = float(s)
    return x


def read_snfit_result(name_or_obj):
    """Read a result file produced by snfit (SALT software).

    Parameters
    ----------
    name_or_obj : str or file-like object
        File name or file object.

    Returns
    -------
    result : Result
    
    """

    if isinstance(name_or_obj, six.string_types):
        f = open(name_or_obj, 'r')
    else:
        f = name_or_obj

    # first line has model name. E.g., Salt2Model
    model = f.readline().strip()

    meta = OrderedDict()
    fitparams = {}
    reading_fitparams = None  # which one are we currently reading?
    for line in f:
        words = line.split()

        # ignore blank lines
        if len(words) == 0:
            continue

        if reading_fitparams is None:

            # is this metadata?
            if words[0][0] == '@':
                meta[words[0][1:]] = _parse_numeric_value(words[1])

            # are we starting a fitparams section?
            elif words[0] == 'BEGIN_OF_FITPARAMS':
                fitparams[words[1]] = OrderedDict()
                reading_fitparams = words[1]

            else:
                raise Exception("cannot parse line: {!r}".format(line))

        # in a fit parameters section
        else:
            # are we at the end?
            if words[0] == 'END_OF_FITPARAMS':
                reading_fitparams = None

            else:
                key, val, err = words[0:3]
                fitparams[reading_fitparams][key] = (_parse_numeric_value(val),
                                                     _parse_numeric_value(err))

    # done reading file
    if isinstance(name_or_obj, six.string_types):
        f.close()

    # parse fitparams section for the model only (assumes this section exists)
    d = fitparams[model]

    # collect parameter names (everything that doesn't start with "Cov"):
    param_names = [key for key in d if not key.startswith("Cov")]
    parameters = np.array([d[name][0] for name in param_names])
    
    # "varied parameters are those that have errors greater than 0.
    vparam_names = [name for name in param_names if d[name][1] > 0.]
    errors = OrderedDict([(name, d[name][1]) for name in vparam_names])

    # construct covariance matrix
    n = len(vparam_names)
    cov = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i, n):
            key1 = "Cov" + vparam_names[i] + vparam_names[j]
            key2 = "Cov" + vparam_names[j] + vparam_names[i]
            if key1 in d:
                cov[i,j] = cov[j,i] = d[key1][0]
            elif key2 in d:
                cov[i,j] = cov[j,i] = d[key2][0]

    return Result(param_names=param_names,
                  parameters=parameters,
                  vparam_names=vparam_names,
                  errors=errors,
                  cov=cov,
                  meta=meta)
