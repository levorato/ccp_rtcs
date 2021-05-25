import numpy as np
import click
import mlflow
import mlflow.sklearn
#import mlflow.statsmodels
import mlflow.keras
import pandas as pd

from pprint import pprint
import logging
import tempfile
import os
import glob


# métricas (funções: métrica(y, y_pred))
from sklearn.metrics import mean_absolute_error, mean_squared_error, max_error, median_absolute_error
#from sklearn.metrics import mean_absolute_percentage_error
from custom_metrics import MAPE_eps, standard_absolute_error, mean_max10_error, mean_absolute_percentage_error
# distribuições de probabilidade customizadas
#from .custom_distributions import random_power

from custom_models import *
from ts_utils import remove_outliers

# Finally, import function to make a machine learning pipeline
from sklearn.pipeline import make_pipeline
#from neuralprophet import NeuralProphet

# fbprophet, pmdarima, statsmodels, tbats
#from hcrystalball.wrappers import SarimaxWrapper
#from hcrystalball.wrappers import ProphetWrapper
#from hcrystalball.wrappers import SimpleSmoothingWrapper
#from hcrystalball.wrappers import HoltSmoothingWrapper
#from hcrystalball.wrappers import TBATSWrapper
#from hcrystalball.wrappers import BATSWrapper
from hcrystalball.wrappers import ExponentialSmoothingWrapper
from hcrystalball.wrappers import get_sklearn_wrapper

# modelos base
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor, AdaBoostRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import LinearRegression, Ridge, ElasticNet, TheilSenRegressor, HuberRegressor, Lasso
from sklearn.linear_model import RANSACRegressor
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_validate
from hcrystalball.model_selection import ModelSelector, FinerTimeSplit
from sklearn.model_selection import GridSearchCV
from hcrystalball.metrics import make_ts_scorer
from hcrystalball.compose import TSColumnTransformer
from sklearn.preprocessing import StandardScaler
from hcrystalball.preprocessing import TargetTransformer
from hcrystalball.feature_extraction import SeasonalityTransformer
from hcrystalball.feature_extraction import HolidayTransformer

from mlflow import log_metric, log_param, log_artifacts


import matplotlib.pyplot as plt

'''
Funções de utilidade geral, por assim dizer
'''

import os
import mlflow
import mlflow.pyfunc

def create_dir(dir):
    '''
    Cria diretório, caso não exista
    '''
    try:
        if not os.path.exists(dir):
            os.makedirs(dir)
    except FileExistsError:
        pass


class FbProphetWrapper(mlflow.pyfunc.PythonModel):
    def __init__(self, model):
        self.model = model
        super().__init__()

    def load_context(self, context):
        from fbprophet import Prophet

        return

    def predict(self, context, model_input):
        future = self.model.make_future_dataframe(periods=model_input["periods"][0])
        return self.model.predict(future)


def yield_artifacts(run_id, path=None):
    """Yield all artifacts in the specified run"""
    client = mlflow.tracking.MlflowClient()
    for item in client.list_artifacts(run_id, path):
        if item.is_dir:
            yield from yield_artifacts(run_id, item.path)
        else:
            yield item.path


def fetch_logged_data(run_id):
    """Fetch params, metrics, tags, and artifacts in the specified run"""
    client = mlflow.tracking.MlflowClient()
    data = client.get_run(run_id).data
    # Exclude system tags: https://www.mlflow.org/docs/latest/tracking.html#system-tags
    tags = {k: v for k, v in data.tags.items() if not k.startswith("mlflow.")}
    artifacts = list(yield_artifacts(run_id))
    return {
        "params": data.params,
        "metrics": data.metrics,
        "tags": tags,
        "artifacts": artifacts,
    }


def print_dict(d, metrics: bool):
    for key in sorted(d.keys()):
        if metrics:
            print('{}: {:.4f}'.format(key, d[key]))
        else:
            print('{}: {}'.format(key, d[key]))


def intersection(lst1, lst2): 
  
    # Use of hybrid method 
    temp = set(lst2) 
    lst3 = [value for value in lst1 if value in temp] 
    return lst3 


def one_step_ahead_test_reselect(best_model,dtf,split_num,target_col):
    lst=[]

    #df=dtf[target_col]
    #df=df.shift(periods=1)
    #df=df.fillna(0)
    #df=df.rename('x')
        
    #dtf=dtf.join(df)
    
    
    for s in range(split_num):
        
        split=split_num-s
        dtf_train = dtf[:-split]
        dtf_test = dtf[-split:]

        X_train, y_train = dtf_train[[col for col in dtf_train.columns if col != target_col]], dtf_train[target_col]
        X_test, y_test = dtf_test[[col for col in dtf_test.columns if col != target_col]], dtf_test[target_col]

        ms,best_model,result=seleciona_modelo_horizonte(dtf_train,target_col,52,1,X_train.columns.tolist(),lags=3)
        
        fitted_model = best_model.fit(X_train, y_train)
        
        predicted_y_test = fitted_model.predict(X_test)
        lst.append(predicted_y_test.iloc[0])
    pred_step_ahead=pd.DataFrame(lst)
    return pred_step_ahead



def one_step_ahead_test(best_model,dtf,split_num,target_col):
    lst=[]

    #df=dtf[target_col]
    #df=df.shift(periods=1)
    #df=df.fillna(0)
    #df=df.rename('x')
        
    #dtf=dtf.join(df)
    
    
    for s in range(split_num):
        
        split=split_num-s
        dtf_train = dtf[:-split]
        dtf_test = dtf[-split:]

        X_train, y_train = dtf_train[[col for col in dtf_train.columns if col != target_col]], dtf_train[target_col]
        X_test, y_test = dtf_test[[col for col in dtf_test.columns if col != target_col]], dtf_test[target_col]

        fitted_model = best_model.fit(X_train, y_train)
        predicted_y_test = fitted_model.predict(X_test)
        lst.append(predicted_y_test.iloc[0])
    pred_step_ahead=pd.DataFrame(lst)
    return pred_step_ahead


'''
    Calcula as metricas de erro (previsto x realizado) para cada uma das categorias de dados passadas como parametro.
'''
def calcula_erros_previsao(cat_list, X, y_true, y_pred,tag,tipo_previsao):
    dict_erro = dict()
    print('Calculando os erros de previsao para as colunas ' + str(cat_list))
    for cat in cat_list + ['Total']:
        #print('Filter on ' + cat)
        if cat == 'Total':
            df_filter = X['_holiday_JP'].notnull()
        elif cat == '_holiday_BR':
            df_filter = (X['_holiday_JP'].str.len() > 0)
        else:
            df_filter = (X[cat] == 1)
        if len(X[df_filter].index) == 0:
            #print('Pulando, devido a dataframe filtrado de tamanho zero.')
            continue
        # end if
        #print(X[df_filter].head(10))
        dict_erro[cat] = []
        dict_erro[cat].append(cat)
        dict_erro[cat].append(mean_absolute_percentage_error(y_true[df_filter], y_pred[df_filter]))
        dict_erro[cat].append(mean_absolute_error(y_true[df_filter], y_pred[df_filter]))
        dict_erro[cat].append(mean_squared_error(y_true[df_filter], y_pred[df_filter]))
        dict_erro[cat].append(median_absolute_error(y_true[df_filter], y_pred[df_filter]))
        dict_erro[cat].append(standard_absolute_error(y_true[df_filter], y_pred[df_filter]))
        dict_erro[cat].append(max_error(y_true[df_filter], y_pred[df_filter]))
        dict_erro[cat].append(mean_max10_error(y_true[df_filter], y_pred[df_filter]))
        dict_erro[cat].append(tag)
        dict_erro[cat].append(tipo_previsao)
    # end for
    return dict_erro, ['Periodo', 'MAPE', 'MAE', 'MSE', 'median_abs_error', 'std_abs_error', 'max_error',
                       'mean_max10_error','tag','tipo_previsao']





'''
    Cria uma tabela (CSV e HTML) com diversas metricas de erro de previsao e armazena no mlflow.  
'''
def loga_tabela_erros_por_categoria(exog_feature_list, X, y, y_pred, split_name,tag,tipo_previsao):
    erros_por_categoria_dict, column_names = calcula_erros_previsao(exog_feature_list, X, y, y_pred,tag,tipo_previsao)
    df_erros_por_categoria = pd.DataFrame.from_dict(erros_por_categoria_dict, orient='index',
                                                    columns=column_names)
    return df_erros_por_categoria

def computa_resultados(exog_features,y_test,y_train,target_col,X_test,regressor):
    prev=loga_tabela_erros_por_categoria(exog_features.columns.tolist(), X_test[X_test.index.isin(y_test.index)], y_test['Real'], y_test[regressor],'teste_previsao',target_col,'previsao')
    baseline=loga_tabela_erros_por_categoria(exog_features.columns.tolist(), X_test[X_test.index.isin(y_test.index)], y_test['Real'], y_test[target_col],'teste_baseline',target_col,'baseline')
    med=y_train.quantile(0.5)
    mediana=loga_tabela_erros_por_categoria(exog_features.columns.tolist(), X_test[X_test.index.isin(y_test.index)], y_test['Real'],np.ones(y_test.shape[0])*med,'teste_mediana',target_col,'mediana')
    listapd=[prev,baseline,mediana]
    resultados_df=pd.concat(listapd,axis=0)
    resultados_df['Algoritmo_vencedor']=regressor
    return resultados_df

#@click.command()
#@click.option("--dataset-tratado-dir", default=os.path.join(os.getcwd(), "datasets", "tratado"))
#@click.option("--source-type", default='PROJECT', type=str)
#@click.option("--source-name", default='Git', type=str)
#@click.option("--split-prop", default=0.8, type=float)
#@click.option("--rank", default=12, type=int)



def resultados_finais(resultados,text,mlflow):
    with tempfile.TemporaryDirectory() as temp_dir:  ## Write csv / html from stats dataframe
        filename = os.path.join(temp_dir+ '_'+text+'.csv')
        resultados.to_csv(filename)
        mlflow.log_artifact(filename)
        filename = os.path.join(temp_dir+ '_'+text+'.html')
        resultados.to_html(filename, index=False, header=True, float_format='%.2f')
        mlflow.log_artifact(filename)


def seleciona_modelo_horizonte(dtf_train='',target_col='',seed=42,horizonte=90,exog_features_list='',lags=3):
    #if horizonte<3:
        #lags=horizonte
    
    ms = ModelSelector(horizon=horizonte,
                       frequency='D',
                       country_code_column=None  # 'country' --> deixar None, se nao ocorre erro de execucao
                       )
    ms.create_gridsearch(sklearn_models=False,
                         n_splits=5,  # 10 cross-validation splits
                         between_split_lag=None,
                         sklearn_models_optimize_for_horizon=False,
                         autosarimax_models=False,  # Autosarimax agora esta funcionando, com pmdarima=1.5.3
                         prophet_models=False,  # Nao ativar, pois usaremos o NeuralProphet em seu lugar
                         tbats_models=False,  # TBATS funcionando OK (pip install tbats)
                         exp_smooth_models=False,  # exp_smooth funcionando OK
                         average_ensembles=False,  # average_ensembles, funcionando OK
                         stacking_ensembles=False,  # Nao vamos usar, demora muito e nao da bom resultado
                         exog_cols=exog_features_list,
                         # exog_cols=None,
                         #holidays_days_before=2,
                         #holidays_days_after=1,
                         #holidays_bridge_days=True,
                         )
    #ms.add_model_to_gridsearch(NeuralProphetWrapper(exog_cols=exog_features.columns.tolist()))
    use_scikit = True
    if use_scikit:
        #xgb_r = get_sklearn_wrapper(XGBRegressor,lags=lags)  # , random_state=seed))
        xgb_r = get_sklearn_wrapper(XGBRegressor,lags=lags,objective="reg:pseudohubererror")  # , random_state=seed))
        xgb_r.name = 'XGBRegressor_HuberLoss'
        ms.add_model_to_gridsearch(xgb_r)
        xgbs_r = get_sklearn_wrapper(XGBRegressor,lags=lags)  # , random_state=seed))
        xgbs_r.name = 'XGBRegressor_SquaredLoss'
        ms.add_model_to_gridsearch(xgbs_r)
        ridge_r = get_sklearn_wrapper(Ridge, lags=lags)
        ridge_r.name = 'RidgeRegressor'
        ms.add_model_to_gridsearch(ridge_r)
        huber_r = get_sklearn_wrapper(HuberRegressor)
        huber_r.name = 'Huber'
        ms.add_model_to_gridsearch(huber_r)
    # Method `select_model` is doing majority of the magic for you - it creates forecast for each combination of
    # columns specified in `partition_columns` and for each of the time series it will run grid_search mentioned
    # above. Optionally once can select list of columns over which the model selection will run in parallel using
    # prefect (`parallel_over_columns`).
    # Required format for data is Datetime index, unsuprisingly numerical column for `target_col_name` all other
    # columns except `partition_columns` will be used as exogenous variables - as additional features for modeling.
    ms.select_model(df=dtf_train,
                    target_col_name=target_col,
                    partition_columns=None,
                    #                 parallel_over_columns=['Assortment'],
                    #                 persist_model_selector_results=False,
                    #                 output_path='my_results',
                    #                 executor = LocalDaskExecutor(),
                    )

    ms.persist_results('results')
    #mlflow.log_metric("score", 0.75)

    # ============================  Train model  ==================================
    result = ms.results[0]
    print('Model selection result: \n', str(result))
    best_model = result.best_model
    return ms,best_model,result


def seleciona_modelo_horizonte2(dtf_train='',target_col='',seed=42,horizonte=90,exog_features_list='',lags=3):
    if horizonte<3:
        lags=horizonte
    
    ms = ModelSelector(horizon=horizonte,
                       frequency='D',
                       country_code_column=None  # 'country' --> deixar None, se nao ocorre erro de execucao
                       )
    ms.create_gridsearch(sklearn_models=False,
                         n_splits=5,  # 10 cross-validation splits
                         between_split_lag=None,
                         sklearn_models_optimize_for_horizon=False,
                         autosarimax_models=False,  # Autosarimax agora esta funcionando, com pmdarima=1.5.3
                         prophet_models=False,  # Nao ativar, pois usaremos o NeuralProphet em seu lugar
                         tbats_models=False,  # TBATS funcionando OK (pip install tbats)
                         exp_smooth_models=False,  # exp_smooth funcionando OK
                         average_ensembles=False,  # average_ensembles, funcionando OK
                         stacking_ensembles=False,  # Nao vamos usar, demora muito e nao da bom resultado
                         exog_cols=exog_features_list,
                         # exog_cols=None,
                         #holidays_days_before=2,
                         #holidays_days_after=1,
                         #holidays_bridge_days=True,
                         )
    #ms.add_model_to_gridsearch(NeuralProphetWrapper(exog_cols=exog_features.columns.tolist()))
    use_scikit = True
    if use_scikit:
        #xgb_r = get_sklearn_wrapper(XGBRegressor,lags=lags)  # , random_state=seed))
        #xgb_r = get_sklearn_wrapper(XGBRegressor,lags=lags,objective="reg:pseudohubererror")  # , random_state=seed))
        #xgb_r.name = 'XGBRegressor_HuberLoss'
        #ms.add_model_to_gridsearch(xgb_r)
        xgbs_r = get_sklearn_wrapper(XGBRegressor,lags=lags)  # , random_state=seed))
        xgbs_r.name = 'XGBRegressor_SquaredLoss'
        ms.add_model_to_gridsearch(xgbs_r)
        ridge_r = get_sklearn_wrapper(Ridge, lags=lags)
        ridge_r.name = 'RidgeRegressor'
        ms.add_model_to_gridsearch(ridge_r)
    # Method `select_model` is doing majority of the magic for you - it creates forecast for each combination of
    # columns specified in `partition_columns` and for each of the time series it will run grid_search mentioned
    # above. Optionally once can select list of columns over which the model selection will run in parallel using
    # prefect (`parallel_over_columns`).
    # Required format for data is Datetime index, unsuprisingly numerical column for `target_col_name` all other
    # columns except `partition_columns` will be used as exogenous variables - as additional features for modeling.
    ms.select_model(df=dtf_train,
                    target_col_name=target_col,
                    partition_columns=None,
                    #                 parallel_over_columns=['Assortment'],
                    #                 persist_model_selector_results=False,
                    #                 output_path='my_results',
                    #                 executor = LocalDaskExecutor(),
                    )

    ms.persist_results('results')
    #mlflow.log_metric("score", 0.75)

    # ============================  Train model  ==================================
    result = ms.results[0]
    print('Model selection result: \n', str(result))
    best_model = result.best_model
    return ms,best_model,result



def seleciona_modelo_horizonte3(dtf_train='',target_col='',seed=42,horizonte=1,exog_features_list='',lags=3,pack=''):
    #if horizonte<3:
     #   lags=horizonte
    
    ms = ModelSelector(horizon=horizonte,
                       frequency='D',
                       country_code_column=None  # 'country' --> deixar None, se nao ocorre erro de execucao
                       )
    ms.create_gridsearch(sklearn_models=False,
                         n_splits=2,  # 10 cross-validation splits
                         between_split_lag=None,
                         sklearn_models_optimize_for_horizon=False,
                         autosarimax_models=False,  # Autosarimax agora esta funcionando, com pmdarima=1.5.3
                         prophet_models=False,  # Nao ativar, pois usaremos o NeuralProphet em seu lugar
                         tbats_models=False,  # TBATS funcionando OK (pip install tbats)
                         exp_smooth_models=False,  # exp_smooth funcionando OK
                         average_ensembles=False,  # average_ensembles, funcionando OK
                         stacking_ensembles=False,  # Nao vamos usar, demora muito e nao da bom resultado
                         exog_cols=exog_features_list,
                         # exog_cols=None,
                         #holidays_days_before=2,
                         #holidays_days_after=1,
                         #holidays_bridge_days=True,
                         )
    #ms.add_model_to_gridsearch(NeuralProphetWrapper(exog_cols=exog_features.columns.tolist()))
    use_scikit = True
    
    huber,ridge,xgb_sq,xgb_hb=pack
    if use_scikit:
        if target_col in xgb_hb:
            xgb_r = get_sklearn_wrapper(XGBRegressor,lags=lags,objective="reg:pseudohubererror")  # , random_state=seed))
            #xgb_r = get_sklearn_wrapper(XGBRegressor,lags=lags,objective="reg:tweedie")  # , random_state=seed))
            xgb_r.name = 'XGBRegressor_Huber'
            ms.add_model_to_gridsearch(xgb_r)
        elif target_col in xgb_sq:
            xgb_r = get_sklearn_wrapper(XGBRegressor,lags=lags)  # , random_state=seed))
            #xgb_r = get_sklearn_wrapper(XGBRegressor,lags=lags,objective="reg:tweedie")  # , random_state=seed))
            xgb_r.name = 'XGBRegressor_Squared_Loss'
            ms.add_model_to_gridsearch(xgb_r)
        elif target_col in ridge:
            ridge_r = get_sklearn_wrapper(Ridge, random_state=seed,lags=lags)
            ridge_r.name = 'Ridge'
            ms.add_model_to_gridsearch(ridge_r)
        else:
            huber_r = get_sklearn_wrapper(HuberRegressor,max_iter=160)
            huber_r.name = 'Huber'
            ms.add_model_to_gridsearch(huber_r)
    # Method `select_model` is doing majority of the magic for you - it creates forecast for each combination of
    # columns specified in `partition_columns` and for each of the time series it will run grid_search mentioned
    # above. Optionally once can select list of columns over which the model selection will run in parallel using
    # prefect (`parallel_over_columns`).
    # Required format for data is Datetime index, unsuprisingly numerical column for `target_col_name` all other
    # columns except `partition_columns` will be used as exogenous variables - as additional features for modeling.
    ms.select_model(df=dtf_train,
                    target_col_name=target_col,
                    partition_columns=None,
                    #                 parallel_over_columns=['Assortment'],
                    #                 persist_model_selector_results=False,
                    #                 output_path='my_results',
                    #                 executor = LocalDaskExecutor(),
                    )

    ms.persist_results('results')
    #mlflow.log_metric("score", 0.75)

    # ============================  Train model  ==================================
    result = ms.results[0]
    print('Model selection result: \n', str(result))
    best_model = result.best_model
    return ms,best_model,result
    