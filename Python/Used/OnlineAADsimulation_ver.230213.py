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
from preprocessing_ha import *


# Direction
path = 'C:/Users/LeeJiWon/Desktop/hykist/인수인계/Recording data_CEDA/'

# Load All speech
ALLSPEECH = np.load(path + 'Allspeech.npy')     # 60 by 3840
# 60 by 3840  /1-30 : Jouney by time / 31-60 Twenty by time // srat : 64
stim_T = ALLSPEECH[:30, :]
stim_J = ALLSPEECH[30:, :]

# ----------------------------- Parameter Setting -----------------------------#
srate = 125
fs = 64
tmin = 0
tmax = 250
Dir = -1
reg_lambda = 10
train = 14
width = 15

# Set int
Indi_Acc = []

# for trial array / for entire array
onset = []

#==============================================================================#
# loop for subjects

for Sub in range(1,10):
    # recreate array per subject
    model_w = np.zeros((15, 17, 1))
    inter_w = np.zeros((15, 1))
    model_train = np.zeros((15, 17, 1))
    inter_train = np.zeros((15, 1))
    acc = []
    dif = []
    ind = []
    acc = np.zeros((16,46))

    # Load RAW and AUX data
    RAW = np.load(path + str(Sub) + '/RAW_' + str(Sub) +'.npy')  # channel * time
    AUX = np.load(path + str(Sub) + '/AUX_' + str(Sub) +'.npy')  # 3 * time

    # Trigger detection
    dif = np.diff(AUX[1, :])
    ind = np.where(dif > 0)
    ind = ind[0] + 1
    for i in range(0, len(ind) - 2):
        try:
            if ind[i + 1] - ind[i] == 1:
                ind = np.delete(ind, i + 1)
        except:
            pass

    # Training
    for tr in range(train):

        speech = ind[tr] + (srate * 3) + 1

        for i in range(46):  # window number 0~45
            # RAW data
            win = RAW[:, speech + (srate * (i)): speech + (srate * (width + i))]

            # ----------------------------- Pre-processing -----------------------------#
            win = np.delete(win, 7, axis=0)
            win = Preproccessing(win, srate, 0.5, 8, 601)  # data, sampling rate, low-cut, high-cut, filter order
            # ------------------------------- Train set -------------------------------#
            # Train decoder model
            model, tlag, inter = mtrf_train(stim_J[tr:tr + 1, fs * (i): fs * (width + i)].T, win.T, fs, Dir,
                                            tmin, tmax, reg_lambda)
            # Sum w - window
            model_w += model
            inter_w += inter

        # Sum tr
        model_train += model_w
        inter_train += inter_w
        print("Train set {0}".format(tr))
    # Average after train trial
    model_train /= (46 * train)
    inter_train /= (46 * train)

    # Testing
    for tr in range(train,30):
        speech = ind[tr] + (srate * 3) + 1  # 3초 후 부터가 speech 이기에 +1

        for i in range(46):  # window number 0~45
            # RAW data
            win = RAW[:, speech + (srate * (i)): speech + (srate * (width + i))]

            # ----------------------------- Pre-processing -----------------------------#
            win = np.delete(win, 7, axis=0)
            win = Preproccessing(win, srate, 0.5, 8, 601)  # data, sampling rate, low-cut, high-cut, filter order
            # ------------------------------- Test set -------------------------------#
            ## Calculate Predicted signal ##
            pred, corr_att, p, mse = mtrf_predict(stim_J[tr:tr + 1, fs * (i): fs * (width + i)].T, win.T, model_train,
                                                 fs, Dir, tmin, tmax, inter_train)
            pred, corr_utt, p, mse = mtrf_predict(stim_T[tr:tr + 1, fs * (i): fs * (width + i)].T, win.T, model_train,
                                                 fs, Dir, tmin, tmax, inter_train)
            # Stock correlation value per window(i)
            # Compare with both correlation values
            if corr_att > corr_utt:
                acc[tr - train, i] = 1
                print("acc=1")
            else:
                acc[tr - train, i] = 0
                print("acc=0")

        print("Test set {0}".format(tr))


    # Save acc for entire Accuracy
    Indi_Acc.append(np.mean(acc) * 100)
    print("Done Subject " + str(Sub))
    print("Accuracy =  " + str(np.mean(acc) * 100))

print("Done")