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
from scipy import signal, io
from scipy.signal import butter, lfilter, resample, filtfilt
# from helper import *
from pymtrf import *
from psychopy import visual, core, event
from PreProcessing import *
from Comments import *
from Direction import *

#################################################################################################
#---------------------------------- Experimental SETTING ---------------------------------------#
#################################################################################################
# Load data / mat file
# set info
path = 'C:/Users/LeeJiWon/Desktop/OpenBCI'          # Base Directory
sub = '6'
raw_mat = io.loadmat(path + '/Recording data/'+ sub + '/RAW_' + sub + '.mat')
raw = raw_mat['RAW']        # channel by time
raw = np.concatenate((raw, np.ones([16,100])), axis=1)  # for final trial (lack of time)
tri_mat = io.loadmat(path + '/Recording data/' + sub + '/TRIGGER_' + sub + '.mat')
tri = tri_mat['TRIGGER']    # 3 by time

# Load All speech
allspeech = np.load(path + '/AAD/Python/Allspeech.npy')
# 60 by 3840  /1-30 : left by time / 31-60 right by time // srat : 64

stim_T = allspeech[:30, :]  # 30 by 3840  Twenty   // trial by time
stim_J = allspeech[30:, :]  # 30 by 3840  Journey  // trial by time

# ----------------------------- Parameter Setting -----------------------------#
srate = 125
fs = 64
tmin = 0
tmax = 250
Dir = -1
reg_lambda = 10
train = 14
##############################################
# Set int

corr_J = []
corr_T = []
Acc = []
ACC = []
model_w = []
inter_w = []
entr_J = []
entr_T = []
Accuracy = []
# for trial array / for entire array
onset = []
tr = 0  # trial

# ==================================================================================================#
# -------------------------------------- START EXPERIMENT ------------------------------------------#
# ==================================================================================================#
# ---------- Start 30 trial ----------#
# Trigger detection
dif = np.diff(tri[1,:])
ind = np.where(dif > 0)
ind = ind[0] + 1
for i in range(0, len(ind) - 2):
    try:
        if ind[i + 1] - ind[i] == 1:
           ind = np.delete(ind,i+1)
    except:
        pass

#onset.append(ind[len(ind) - 1])

while tr < 30:  # 30
    # Format per trial
    i = 0  # Window number
    Acc = []
    # ----------------------------- Working while 60s -----------------------------#
    speech = ind[tr] + (srate * 3) + 1  # 3초 후 부터가 speech 이기에 +1
    while i != 46:  # window number 0~45
        # RAW data
        win = raw[:, speech+srate*(i) : speech+srate*(15 + i)]

        # ----------------------------- Pre-processing -----------------------------#
        win = np.delete(win, 7, axis=0)
        win = Preproccessing(win, srate, 0.5, 8, 601)  # data, sampling rate, low-cut, high-cut, filter order
        # ------------------------------- Train set -------------------------------#
        if tr < train:  # int train
            state = "Train set"
            # Train decoder model
            model, tlag, inter = mtrf_train(stim_J[tr:tr + 1, 64 * (i) : 64 * (15 + i)].T, win.T, fs, Dir,
                                            tmin, tmax, reg_lambda)
            # Sum w - window
            if i == 0:
                model_w = model
                inter_w = inter
            else:  # i > 0 - 45까지
                model_w = np.add(model_w, model)
                inter_w = np.add(inter_w, inter)

        # ------------------------------- Test set -------------------------------#
        else:
            state = "Test set"
            ## Calculate Predicted signal ##
            pred_j, corr_j, p, mse = mtrf_predict(stim_J[tr:tr+1, 64 * (i) : 64 * (15 + i)].T, win.T, model, fs,
                                             Dir, tmin, tmax, inter)
            pred_t, corr_t, p, mse = mtrf_predict(stim_T[tr:tr+1, 64 * (i) : 64 * (15 + i)].T, win.T, model, fs,
                                             Dir, tmin, tmax, inter)
            # Stock correlation value per window(i)
            corr_J = np.append(corr_J, corr_j)
            corr_T = np.append(corr_T, corr_t)

            ######  Estimate accuracy  #####
            if corr_j > corr_t:
                acc = 1
            else:
                acc = 0

            # Save acc for entire Accuracy
            Acc = np.append(Acc, acc)

        # Plus window number
        i = i + 1

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
        print("Train set {0}".format(tr))

        #=============================================================
    elif state == "Test set":
        # Stack correlation value collected during one trial
        entr_J.append(corr_J)
        entr_T.append(corr_T)
        corr_J = []
        corr_T = []

        # Collect Accuracy per trial
        ACC = np.append(ACC, np.mean(Acc))
        print("\n==================================\n")
        print("Test set {0}".format(tr+1))
        print("Present Accuracy = {0}%".format(ACC[-1] * 100))
        print("\n==================================\n")

    tr = tr + 1

# ----------------------------- 30 trial End -----------------------------#
print("Total Accuracy = {0}%".format(mean(ACC[14:30])*100))
Accuracy.append(ACC)
print("The End")

#### save ####
# mat save
# Save per trial // eeg, trigger, accuracy ,behavior
Accuracy = np.asarray(Accuracy)
np.save(path + '/save_data/Accuracy_EMA_' + sub, Accuracy)
np.save(path + '/save_data/AllcorrJ_EMA_' + sub, entr_J)
np.save(path + '/save_data/AllcorrT_EMA_' + sub, entr_T)



