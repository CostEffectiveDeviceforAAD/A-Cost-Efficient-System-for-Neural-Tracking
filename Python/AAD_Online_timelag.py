""""""""""""""""""""""""""""""""""""""""""
#        OpenBCI - Brainflow           #
#           For Online AAD              #
""""""""""""""""""""""""""""""""""""""""""

###### Imports #####
import librosa, warnings, random, time, os, sys, serial, logging, argparse, mne, scipy.io, math
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from pylsl import StreamInlet, resolve_stream, StreamInfo
# from OpenBCI_lsl import *
from scipy import signal, io, linalg
from scipy.signal import butter, lfilter, resample, filtfilt
from pymtrf import *
from psychopy import visual, core, event
from preprocessing_ha import *
from Comments import *
from Direction import *
from helper import *
from scipy.stats import pearsonr
import warnings


# -------------------------------- SETTING ---------------------------------#

#########
subject = '_timelag'
###########

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



for sub in range(2,12):
    # Load data
    #sub = '0806_LJH'
    sub = str(sub)
    raw_mat = io.loadmat(path + '/Recording data/'+ sub + '/RAW_' + sub + '.mat')
    raw = raw_mat['RAW']        # channel by time
    raw = np.concatenate((raw, np.ones([16,100])), axis=1)  # for final trial (lack of time)
    tri_mat = io.loadmat(path + '/Recording data/' + sub + '/TRIGGER_' + sub + '.mat')
    tri = tri_mat['TRIGGER']    # 3 by time

    # ----------------------------- Parameter Setting -----------------------------#
    #srate = 1000
    srate = 125
    fs = 64
    tmin = 0
    tmax = 400
    Dir = -1
    reg_lambda = 10
    train = 14
    ##############################################
    # Set int

    r_L = []
    r_R = []
    Acc = []
    ACC = []
    model_w = []
    inter_w = []
    entr_L = []
    entr_R = []
    Accuracy = []

    # for trial array / for entire array

    onset = []
    predic_l = []
    predic_r = []
    Pre_L = []
    Pre_R = []
    w = 1  # To avoid repeat when detect trigger
    j = 0  # Question number
    tr = 0  # trial
    ##############################################


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


    for idx in range(0,27):

        for tr in range(0,30):  # 30

            # ----------------------------- Trigger detection -----------------------------#
            # per trial

            # Format per trial
            i = 0  # Window number
            Acc = []
            predic_l = []
            predic_r = []

            #raw = eeg[:60,:,tr]
            # ----------------------------- Working while 60s -----------------------------#

            # onset 부터 3초 지나고 원하는 시간(한 trial) 동안 돌아가도록
            speech = onset[tr] + (srate * 3) + 1  # 3초 후 부터가 speech 이기에 +1
            #speech = 0
            while i != 46:  # 46 번째 window 까지 (0~45)

                # RAW data

                win = raw[:, speech+srate*(i) : speech+srate*(15 + i)]

                # ----------------------------- Pre-processing -----------------------------#
                # preprocessing_ha.py

                win = np.delete(win, 7, axis=0)
                win = Preproccessing(win, srate, 0.5, 8, 601)  # data, sampling rate, low-cut, high-cut, filter order

                # ------------------------------- Train set -------------------------------#
                if tr < train:  # int train
                    state = "Train set"

                    #################################
                    def mtrf_train(stim, resp, fs, mapping_direction, tmin, tmax, reg_lambda,idx):

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
                        lags = lags_range[26-idx:27]
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

                ##################

                    ## mTRF train function ##
                    model, tlag, inter = mtrf_train(stim_L[tr:tr + 1, 64 * (i) : 64 * (15 + i)].T, win.T, fs, Dir,
                                                    tmin, tmax, reg_lambda, idx)

                    'model - (16,17,1)  / tlag - (17,1) / inter - (16,1)'   #channel by tau by 1

                    #===========================================================================

                    # Sum w - window
                    if i == 0:
                        model_w = model
                        inter_w = inter
                    else:  # i > 0 - 45까지
                        model_w = np.add(model_w, model)
                        inter_w = np.add(inter_w, inter)

                    i = i + 1

                # ------------------------------- Test set -------------------------------#
                else:
                    state = "Test set"

                ##################################
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
                        lags = lags_range[26-idx:27]
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

                    #########################

                    ## Calculate Predicted signal ##
                    pred_l, r_l, p, mse = mtrf_predict(stim_L[tr:tr+1, 64 * (i) : 64 * (15 + i)].T, win.T, model, fs,
                                                     Dir, tmin, tmax, inter, idx)
                    pred_r, r_r, p, mse = mtrf_predict(stim_R[tr:tr+1, 64 * (i) : 64 * (15 + i)].T, win.T, model, fs,
                                                     Dir, tmin, tmax, inter, idx)

                    print("Test")

                    # Stock correlation value per window(i)
                    r_L = np.append(r_L, r_l)
                    r_R = np.append(r_R, r_r)

                    ######  Estimate accuracy  #####
                    if r_l > r_r:
                        acc = 1
                    else:
                        acc = 0

                    print("======= acc : {0} ".format(acc))

                    # Save acc for entire Accuracy
                    Acc = np.append(Acc, acc)

                    # Up window number
                    i = i + 1

                # -- End one windo --w

            # ------------------------ End 60s - one trial ------------------------#

            ###### The things that have to calculate per trial ######
            ## Add model_w case train
            if state == "Train set":
                if tr == 0:
                    model_wt = model_w
                    inter_wt = inter_w
                elif tr > 0:
                    model_wt = np.add(model_wt, model_w)
                    inter_wt = np.add(inter_wt, inter_w)

                # Average at last train trial
                if tr == 13:
                    model_wm = model_wt / (46 * 14)
                    inter_wm = inter_wt / (46 * 14)
                    model = model_wm
                    inter = inter_wm

                #=============================================================

                print("Train_check_trial")  # 당연히 높게 나와야함
                ACC = np.append(ACC, np.mean(Acc))
                Pre_L.append(predic_l)
                Pre_R.append(predic_r)
                print("\n==================================\n")
                print("Present Accuracy = {0}%".format(ACC[-1] * 100))
                print("\n==================================\n")

                #=============================================================
            elif state == "Test set":
                # Stack correlation value collected during one trial
                entr_L.append(r_L)
                entr_R.append(r_R)

                r_L = []
                r_R = []

                # Collect Accuracy per trial
                ACC = np.append(ACC, np.mean(Acc))
                Pre_L.append(predic_l)
                Pre_R.append(predic_r)
                print("\n==================================\n")
                print("Present Accuracy = {0}%".format(ACC[-1] * 100))
                print("\n==================================\n")

            print(tr)

            w = 1

        Accuracy.append(ACC)
        print("done one trial")
        #print("Accuracy = {0}%".format(mean(ACC[0,14:]) * 100))
        print("Present index = {0}".format(idx))
        ACC = []

    Accuracy2 = np.asarray(Accuracy)
    scipy.io.savemat(path + '/save_data/Accuracy' + subject + '_' + sub + '.mat', {'Acc': Accuracy2})
# ----------------------------- 30 trial End -----------------------------#


print("The End")

#### save ####
# mat save
# Save per trial // eeg, trigger, accuracy ,behavior
Accuracy2 = np.asarray(Accuracy)

scipy.io.savemat(path + '/save_data/Accuracy' + subject + '_' + sub + '.mat', {'Acc': Accuracy2})



