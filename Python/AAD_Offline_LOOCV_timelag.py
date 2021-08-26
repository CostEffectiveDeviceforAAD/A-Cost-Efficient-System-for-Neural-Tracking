""""""""""""""""""""""""""""""""""""""""""
#        OpenBCI - Brainflow           #
#           For Offline AAD              #
""""""""""""""""""""""""""""""""""""""""""

###### Imports #####
import librosa, warnings, random, time, os, sys, serial, logging, argparse, mne, scipy.io, math
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from pylsl import StreamInlet, resolve_stream, StreamInfo
# from OpenBCI_lsl import *
from scipy import signal, io
from scipy.signal import butter, lfilter, resample, filtfilt
from helper import *
from pymtrf import *
from psychopy import visual, core, event
from preprocessing_ha import *
from Comments import *
from Direction import *
from scipy.stats import pearsonr

# -------------------------------- SETTING ---------------------------------#


#loc = 'kist'
loc = 'hyu'

if loc == 'kist':
    path = 'C:/Users/LeeJiWon/Desktop/OpenBCI'
elif loc == 'hyu':
    path = 'C:/Users/user/Desktop/hy-kist/OpenBCI'


# ----------------------------- Load Speech segment data ------------------------------#

# Load All speech
allspeech = np.load(path + '/AAD/Python/Allspeech.npy')
# 60 by 3840  /1-30 : left by time / 31-60 right by time // srat : 64

stim_R = allspeech[:30, :]  # 30 by 3840  Twenty   // trial by time
stim_L = allspeech[30:, :]  # 30 by 3840  Journey  // trial by time


##############################################

def mtrf_train(stim, resp, fs, mapping_direction, tmin, tmax, reg_lambda, idx):
    if stim.shape[0] < stim.shape[1]:
        warnings.warn(f'stim: more features {stim.shape[0]} ' +
                      f'than samples {stim.shape[1]}, check input dimensions!')
    if resp.shape[0] < resp.shape[1]:
        warnings.warn(f'resp: more features {resp.shape[0]} ' +
                      f'than samples {resp.shape[1]}, check input dimensions!')

    assert tmin < tmax, 'tmin has to be smaller than tmax'
    assert reg_lambda > 0, 'reg_lambda has to be positive and larger than 0!'

    x, y, tmin, tmax = stimulus_mapping(mapping_direction, stim, resp, tmin, tmax)

    t_min = np.floor(tmin / 1e3 * fs * mapping_direction).astype(int)
    t_max = np.ceil(tmax / 1e3 * fs * mapping_direction).astype(int)

    ##
    lags_range = lag_builder(t_min, t_max)
    lags = lags_range[32 - idx:33 - idx]
    #lags = lag_builder(t_min, t_max)
    ##

    lag_x = lag_gen(x, lags)
    x_input = np.hstack([np.ones(x.shape), lag_x])
    n_feat = x_input.shape[1]

    if x.shape[1] == 1:
        reg_matrix = quadratic_regularization(n_feat)
    else:
        reg_matrix = np.eye(n_feat)

    coefficients = regularized_regression_fit(x_input, y, reg_matrix, reg_lambda)
    model, inter = coefficient_to_model(coefficients, x.shape[1],
                                        lags.shape[0], y.shape[1])
    time_lags = lags / fs * 1e3

    return model, time_lags, inter


def mtrf_predict(stim, resp, model, fs, mapping_direction, tmin, tmax, constant, idx):
    # Define x and y
    assert tmin < tmax, 'Value of tmin must be < tmax'

    if constant is None:
        constant = np.zeros((model.shape[0], model.shape[2]))
    else:
        assert np.all(constant.shape == np.array([model.shape[0],
                                                  model.shape[2]]))

    x, y, tmin, tmax = stimulus_mapping(mapping_direction, stim, resp, tmin, tmax)

    t_min = np.floor(tmin / 1e3 * fs * mapping_direction).astype(int)
    t_max = np.ceil(tmax / 1e3 * fs * mapping_direction).astype(int)

    #lags = lag_builder(t_min, t_max)
    lags_range = lag_builder(t_min, t_max)
    lags = lags_range[32 - idx:33 - idx]  # [26-idx:27]
    ##

    x_lag = np.hstack([np.ones(x.shape), lag_gen(x, lags)])

    model = model_to_coefficients(model, constant)

    pred = regularized_regression_predict(x_lag, model)

    # Calculate accuracy
    if y is not None:
        r = np.zeros((1, y.shape[1]))
        p = np.zeros((1, y.shape[1]))
        mse = np.zeros((1, y.shape[1]))
        for i in range(y.shape[1]):
            r[:, i], p[:, i] = pearsonr(y[:, i], pred[:, i])
            mse[:, i] = np.mean((y[:, i] - pred[:, i]) ** 2)
    else:
        r = None
        p = None
        mse = None

    return pred, r, p, mse

############################################

####

for s in range(8,12):

    sub = str(s)
    # Load data

    raw_mat = io.loadmat(path + '/Recording data/'+ sub + '/RAW_' + sub + '.mat')
    raw = raw_mat['RAW']        # channel by time
    raw = np.concatenate((raw, np.ones([16,100])), axis=1)  # for final trial (lack of time)
    tri_mat = io.loadmat(path + '/Recording data/' + sub + '/TRIGGER_' + sub + '.mat')
    tri = tri_mat['TRIGGER']    # 3 by time

    # ----------------------------- Parameter Setting -----------------------------#
    #srate = 1000
    srate = 125
    fs = 64
    tmin = -250
    tmax = 250
    Dir = -1
    reg_lambda = 10

    ##############################################
    # Set int

    r_L = []
    r_R = []
    Acc = []
    ACC = []
    Accuracy = []
    onset = []
    tr = 0
    train = 0

    model = np.zeros([16, 17, 1])
    inter = np.zeros([16, 1])


    # ==================================================================================================#
    # -------------------------------------- START EXPERIMENT ------------------------------------------#
    # ==================================================================================================#


    # ---------- Start 30 trial ----------#

    # Trigger detection
    dif = np.diff(tri[1,:])
    ind = np.where(dif > 0)
    ind = ind[0] + 1

    for i in range(0, len(ind) - 1):
        if ind[i + 1] - ind[i] != 1:
            onset.append(ind[i])
    onset.append(ind[len(ind) - 1])


    for idx in range(0,33):

        for train in range(0,30):  # 30 / select train set

            # Format per trial

            model_w = np.zeros([15,1,1])
            inter_w = np.zeros([15,1])
            # ------------------------------- Train set -------------------------------#
            for tr in range(0,30):

                speech = onset[tr] + (srate * 3) + 1  # 3초 후 부터가 speech 이기에 +1
                win = raw[:, speech: speech + srate * 60]
                win = np.delete(win, 7, axis=0)
                win = Preproccessing(win, srate, 0.5, 8, 601)

                if tr != train:

                    ## mTRF train function ##
                    model, tlag, inter = mtrf_train(stim_L[tr:tr+1,:].T, win.T, fs, Dir,
                                                    tmin, tmax, reg_lambda, idx)

                    'model - (channel,17,1)  / tlag - (17,1) / inter - (16,1)'

                    model_w = np.add(model_w, model)
                    inter_w = np.add(inter_w, inter)

                    print("Train - {0}".format(tr))

            model = model_w / (29)
            inter = inter_w / (29)


            #===========================================================================
            # ------------------------------- Test set -------------------------------#
            speech = onset[train] + (srate * 3) + 1  # 3초 후 부터가 speech 이기에 +1
            win = raw[:, speech: speech + srate * 60]
            win = np.delete(win, 7, axis=0)
            win = Preproccessing(win, srate, 0.5, 8, 601)

            ## Calculate Predicted signal ##
            pred_l, r_l, p, mse = mtrf_predict(stim_L[train:train+1,:].T, win.T, model, fs,
                                             Dir, tmin, tmax, inter, idx)
            pred_r, r_r, p, mse = mtrf_predict(stim_R[train:train+1,:].T, win.T, model, fs,
                                             Dir, tmin, tmax, inter, idx)

            print("Test - {0}".format(train))

            ######  Estimate accuracy  #####
            if r_l > r_r:
                acc = 1
            else:
                acc = 0

            # Save acc for entire Accuracy
            Acc = np.append(Acc, acc)   # 하나의 corr 씩 30개

            print("done one index - {0}".format(idx))


        Accuracy.append(mean(Acc))    # 한 time-lag 지점 마다의 acc(30개) 쌓기   # 33개의 individual timelag point acc 가 쌓임
        Acc = []

    Accuracy = np.asarray(Accuracy)
    scipy.io.savemat(path + '/save_data/Accuracy' + '_timelag_off_' + sub + '.mat', {'Acc': Accuracy})
    Accuracy = []


print("The End")



