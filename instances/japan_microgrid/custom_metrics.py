# https://scikit-learn.org/stable/modules/model_evaluation.html#the-scoring-parameter-defining-model-evaluation-rules

import numpy as np
import heapq
from sklearn.metrics import make_scorer
from sklearn.metrics._regression import _check_reg_targets
from sklearn.utils.validation import (check_array, check_consistent_length,
                                _num_samples, _check_sample_weight)


def MAPE_eps(y_true, y_pred, sample_weight=None, multioutput='uniform_average'):
    '''
    Mean Absolute Percentual Error, com um pequeno epsilon no denominador
    '''
    eps = 0.001
    return np.mean(np.abs(y_true - y_pred) / (y_true + eps))


# make_scorer é usado para passar ao RandomizedSearchCV, se quiser. Mas algo não funciona direito.
# MAPE_eps_score = make_scorer(MAPE_eps, greater_is_better=False)

# Sugestão de mais uma métrica: desvio padrão do erro absoluto
# Favor fazer eventuais correções no código
def standard_absolute_error(y_true, y_pred, sample_weight=None, multioutput='uniform_average'):
    '''
    Desvio padrão do erro
    '''
    #return np.std(np.abs(y_true - y_pred))

    y_type, y_true, y_pred, multioutput = _check_reg_targets(
        y_true, y_pred, multioutput)
    std_errors = np.nanstd(np.abs(y_pred - y_true), axis=0)
    output_errors = np.average(std_errors,
                               weights=sample_weight, axis=0)

    if isinstance(multioutput, str):
        if multioutput == 'raw_values':
            return output_errors
        elif multioutput == 'uniform_average':
            # pass None as weights to np.average: uniform mean
            multioutput = None

    return np.average(output_errors, weights=multioutput)


std_score = make_scorer(standard_absolute_error, greater_is_better=False)

def mean_max10_error(y_true, y_pred, sample_weight=None, multioutput='uniform_average'):
    '''
    Calcula a media dos 10 maiores erros
    '''
    y_type, y_true, y_pred, multioutput = _check_reg_targets(
        y_true, y_pred, multioutput)
    max10 = heapq.nlargest(10, np.abs(y_true - y_pred))
    output_errors = np.average(max10,
                               weights=sample_weight, axis=0)

    if isinstance(multioutput, str):
        if multioutput == 'raw_values':
            return output_errors
        elif multioutput == 'uniform_average':
            # pass None as weights to np.average: uniform mean
            multioutput = None

    return np.average(output_errors, weights=multioutput)


def mean_absolute_percentage_error(y_true, y_pred,
                                   sample_weight=None,
                                   multioutput='uniform_average'):
    """Mean absolute percentage error regression loss.
    Note here that we do not represent the output as a percentage in range
    [0, 100]. Instead, we represent it in range [0, 1/eps]. Read more in the
    :ref:`User Guide <mean_absolute_percentage_error>`.
    .. versionadded:: 0.24
    Parameters
    ----------
    y_true : array-like of shape (n_samples,) or (n_samples, n_outputs)
        Ground truth (correct) target values.
    y_pred : array-like of shape (n_samples,) or (n_samples, n_outputs)
        Estimated target values.
    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights.
    multioutput : {'raw_values', 'uniform_average'} or array-like
        Defines aggregating of multiple output values.
        Array-like value defines weights used to average errors.
        If input is list then the shape must be (n_outputs,).
        'raw_values' :
            Returns a full set of errors in case of multioutput input.
        'uniform_average' :
            Errors of all outputs are averaged with uniform weight.
    Returns
    -------
    loss : float or ndarray of floats in the range [0, 1/eps]
        If multioutput is 'raw_values', then mean absolute percentage error
        is returned for each output separately.
        If multioutput is 'uniform_average' or an ndarray of weights, then the
        weighted average of all output errors is returned.
        MAPE output is non-negative floating point. The best value is 0.0.
        But note the fact that bad predictions can lead to arbitarily large
        MAPE values, especially if some y_true values are very close to zero.
        Note that we return a large value instead of `inf` when y_true is zero.
    """
    y_type, y_true, y_pred, multioutput = _check_reg_targets(
        y_true, y_pred, multioutput)
    check_consistent_length(y_true, y_pred, sample_weight)
    epsilon = np.finfo(np.float64).eps
    mape = np.abs(y_pred - y_true) / np.maximum(np.abs(y_true), epsilon)
    output_errors = np.average(mape,
                               weights=sample_weight, axis=0)
    if isinstance(multioutput, str):
        if multioutput == 'raw_values':
            return output_errors
        elif multioutput == 'uniform_average':
            # pass None as weights to np.average: uniform mean
            multioutput = None

    return np.average(output_errors, weights=multioutput)

