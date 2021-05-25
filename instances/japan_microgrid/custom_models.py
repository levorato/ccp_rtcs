import numpy as np
from ts_utils import *
from sklearn.base import RegressorMixin, BaseEstimator


class NaiveModel(BaseEstimator):
    '''
    Calcula e prediz a média calculada
    '''
    def __init__(self, param=None):
        super().__init__()
        self.param = param  # ghost parameter

    def fit(self, X, y):
        self._y = np.mean(y)
        return self

    def predict(self, X):
        return self._y * np.ones(len(X))
    
    # def get_params(self, deep=True):
    #     return {'param': self.param}

    # def set_params(self, **parameters):
    #     for parameter, value in parameters.items():
    #         setattr(self, parameter, value)


class RandomWalkModel(BaseEstimator):
    '''
    Generate a Random Walk process.
    '''
    def __init__(self, param=None):
        super().__init__()
        self.param = param  # ghost parameter

    def fit(self, X, y):
        self._y = y.copy()
        return self

    def predict(self, X):
        diff_ts = self._y - self._y.shift(1)
        rw = utils_generate_rw(y0=self._y[-1], n=len(X), sigma=diff_ts.std(), ymin=self._y.min(), ymax=self._y.max())
        return rw

    # def get_params(self, deep=True):
    #     return {'param': self.param}

    # def set_params(self, **parameters):
    #     for parameter, value in parameters.items():
    #         setattr(self, parameter, value)


