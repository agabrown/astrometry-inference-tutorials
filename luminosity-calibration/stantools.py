"""
Utility functions for using PyStan.

Anthony Brown June 2017
"""
__all__ = ['load_stan_code', 'stan_cache']

import pystan
import pickle
from hashlib import md5

def load_stan_code(filename):
    """
    Load the Stan model code from the input file.

    Parameters
    ----------

    filename - Name of file with model code.

    Returns
    -------

    A string with the model code.
    """
    f = open(filename, 'r')
    contents = f.read()
    f.close()
    return contents

def stan_cache(model_code, model_name=None):
    """
    From the given model code create a pystan.StanModel instance by either compiling the model from
    scratch or by loading a cached version. This is to avoid repeated Stan code compilation in the use of
    PyStan. Code modified from: https://pystan.readthedocs.io/en/latest/avoiding_recompilation.html

    Parameters
    ----------

    model_code - String with Stan model code.

    Keywords
    --------

    model_name - Name of the model.

    Returns
    -------

    Instance of pystan.StanModel.
    """
    code_hash = md5(model_code.encode('ascii')).hexdigest()
    if model_name is None:
        cache_fn = 'cached-model-{}.pkl'.format(code_hash)
    else:
        cache_fn = 'cached-{}-{}.pkl'.format(model_name, code_hash)
    try:
        sm = pickle.load(open(cache_fn, 'rb'))
    except:
        sm = pystan.StanModel(model_code=model_code)
        with open(cache_fn, 'wb') as f:
            pickle.dump(sm, f)
    else:
        print("Using cached StanModel")
    return sm
