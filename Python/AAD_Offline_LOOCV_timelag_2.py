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

def mtrf_crossval(stim, stim2, resp, fs, mapping_direction, tmin, tmax, reg_lambda, idx):

    assert tmin < tmax, 'Value of tmin must be < tmax'

    x, y_un, tmin, tmax = stimulus_mapping(mapping_direction, stim2, resp, tmin, tmax)
    x, y, tmin, tmax = stimulus_mapping(mapping_direction, stim, resp, tmin, tmax)


    #n_trials, n_feat = test_input_dimensions(x_un)
    n_trials, n_feat = test_input_dimensions(x_at)

    n_trl_y, n_targets = test_input_dimensions(y)

    assert n_trials == n_trl_y, 'stim and resp should have the same no of trials!'

    reg_lambda = test_reg_lambda(reg_lambda)

    n_lambda = len(reg_lambda)

    t_min = np.floor(tmin / 1e3 * fs * mapping_direction).astype(int)
    t_max = np.ceil(tmax / 1e3 * fs * mapping_direction).astype(int)

    ##
    lags_range = lag_builder(t_min, t_max)
    lags = lags_range[32 - idx:33 - idx]
    #lags = lag_builder(t_min, t_max)
    ##

    # Set up regularisation
    dim1 = n_feat * lags.shape[0] + n_feat
    model = np.zeros((n_trials, n_lambda, dim1, n_targets))

    if n_feat == 1:
        reg_matrix = quadratic_regularization(dim1)
    else:
        reg_matrix = np.eye(dim1)

    # Training
    x_input = []

    for c_trials in range(n_trials):
        # Generate lag matrix
        x_input.append(np.hstack([np.ones(x[c_trials].shape), lag_gen(x[c_trials], lags)]))
        # Calculate model for each lambda value
        for c_lambda in range(n_lambda):
            temp = regularized_regression_fit(x_input[c_trials],
                                              y[c_trials], reg_matrix, reg_lambda[c_lambda])
            model[c_trials, c_lambda, :, :] = temp

    r = np.zeros((n_trials, n_lambda, n_targets))
    r_un = np.zeros((n_trials, n_lambda, n_targets))
    p = np.zeros(r.shape)
    mse = np.zeros(r.shape)
    pred = []

    for trial in range(n_trials):
        pred.append(np.zeros((n_lambda, y[trial].shape[0], n_targets)))

        # Perform cross-validation for each lambda value
        for c_lambda in range(n_lambda):
            # Calculate prediction
            cv_coef = np.mean(model[np.arange(n_trials) != trial, c_lambda, :, :], 0, keepdims=False)
            pred[trial][c_lambda, :, :] = regularized_regression_predict(x_input[trial], cv_coef)

            # Calculate accuracy
            for k in range(n_targets):
                temp_pred = np.squeeze(pred[trial][c_lambda, :, k]).T
                r[trial, c_lambda, k], p[trial, c_lambda, k] = pearsonr(y[trial][:, k], temp_pred)
                r_un[trial, c_lambda, k], p[trial, c_lambda, k] = pearsonr(y_un[trial][:, k], temp_pred)

                mse[trial, c_lambda, k] = np.mean((y[trial][:, k] - temp_pred) ** 2)

    return r, r_un, p, mse, pred, model

############################################

####

for s in range(2,12):

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
    reg_lambda = [10]

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

    ### stim /resp
    resp = []
    stim_at = []
    stim_un = []
    for tr in range(0, 30):  # 30 trial eeg data 쌓기

        speech = onset[tr] + (srate * 3) + 1  # 3초 후 부터가 speech 이기에 +1
        win = raw[:, speech: speech + srate * 60]
        win = np.delete(win, 7, axis=0)
        win = Preproccessing(win, srate, 0.5, 8, 601)

        resp.append(win.T)  # 3840 by 15 가 30개

    resp = np.asarray(resp)  # 30 by 3840 by 15

    stim_at.append(stim_L.T)
    stim_at = np.asarray(stim_at).T
    stim_un.append(stim_R.T)
    stim_un = np.asarray(stim_un).T

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

        # ------------------------------- Train set -------------------------------#

        ## mTRF train function ##
        # stim / resp = [n_trials, n_samples(time), n_features]
        r, r_un, p, mse, pred, model = mtrf_crossval(stim_at, stim_un, resp, fs, Dir,
                                                tmin, tmax, reg_lambda, idx)

        ######  Estimate accuracy  #####

        for i in range(0,30):

            if r[i] > r_un[i]:
                acc = 1
            else:
                acc = 0
            # Save acc for entire Accuracy
            Acc = np.append(Acc, acc)   # 하나의 corr 씩 30개

        print("done one index - {0}".format(idx))

        Accuracy.append(mean(Acc))    # 한 time-lag 지점 마다의 acc(30개) 쌓기   # 33개의 individual timelag point acc 가 쌓임
        Acc = []

    Accuracy = np.asarray(Accuracy)
    scipy.io.savemat(path + '/save_data/Accuracy' + '_timelag_off_2_' + sub + '.mat', {'Acc': Accuracy})
    Accuracy = []


print("The End")



