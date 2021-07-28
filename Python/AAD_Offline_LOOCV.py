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
# from helper import *
from pymtrf import *
from psychopy import visual, core, event
from preprocessing_ha import *
from Comments import *
from Direction import *


# -------------------------------- SETTING ---------------------------------#

#########
subject = '_loocv'
###########

loc = 'kist'
#loc = 'hyu'

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

# Load data

raw_mat = io.loadmat(path + '/Recording data/0714_LJW/RAW_01LJW.mat')
raw = raw_mat['RAW']        # channel by time
raw = np.concatenate((raw, np.ones([16,100])), axis=1)  # for final trial (lack of time)
tri_mat = io.loadmat(path + '/Recording data/0714_LJW/TRIGGER_01LJW.mat')
tri = tri_mat['TRIGGER']    # 3 by time

ch = 2
'''
raw_mat = io.loadmat('C:/Users/user/Desktop/hy-kist/OpenBCI/AAK/Seg.mat')
eeg = raw_mat['eeg']    # channel by time by trial

'''
# ----------------------------- Parameter Setting -----------------------------#
#srate = 1000
srate = 125
fs = 64
tmin = 0
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

#for ch in [5,11]:

while train < 30:  # 30

    # Format per trial
    Acc = []
    model_w = np.zeros([16,17,1])
    inter_w = np.zeros([16,1])
    # ------------------------------- Train set -------------------------------#
    for tr in range(0,30):

        speech = onset[tr] + (srate * 3) + 1  # 3초 후 부터가 speech 이기에 +1
        win = raw[:, speech: speech + srate * 60]
        #win = np.delete(win, 7, axis=0)
        win = Preproccessing(win, srate, 0.5, 8, 601)

        if tr != train:
            ## mTRF train function ##
            model, tlag, inter = mtrf_train(stim_L[tr:tr+1,:].T, win.T, fs, Dir,
                                            tmin, tmax, reg_lambda)

            'model - (16,17,1)  / tlag - (17,1) / inter - (16,1)'

        model_w = np.add(model_w, model)
        inter_w = np.add(inter_w, inter)

    model = model_w / (29)
    inter = inter_w / (29)


    #===========================================================================
    # ------------------------------- Test set -------------------------------#
    speech = onset[train] + (srate * 3) + 1  # 3초 후 부터가 speech 이기에 +1
    win = raw[:, speech: speech + srate * 60]
    #win = np.delete(win, 7, axis=0)
    win = Preproccessing(win, srate, 0.5, 8, 601)

    ## Calculate Predicted signal ##
    pred_l, r_l, p, mse = mtrf_predict(stim_L[train:train+1,:].T, win.T, model, fs,
                                     Dir, tmin, tmax, inter)
    pred_r, r_r, p, mse = mtrf_predict(stim_R[train:train+1,:].T, win.T, model, fs,
                                     Dir, tmin, tmax, inter)

    print("Test")

    ######  Estimate accuracy  #####
    if r_l > r_r:
        acc = 1
    else:
        acc = 0

    print("======= acc : {0} ".format(acc))

    # Save acc for entire Accuracy
    Acc = np.append(Acc, acc)

    ACC.append(Acc)
    print("\n==================================\n")
    #print("Present Accuracy = {0}%".format(ACC[-1] * 100))
    print("Present trial = {0}".format(train))
    print("\n==================================\n")

    train = train + 1


Accuracy.append(ACC)

print("The End")

#### save ####
# mat save
# Save per trial // eeg, trigger, accuracy ,behavior
Accuracy = np.asarray(Accuracy)
scipy.io.savemat(path + '/save_data/Accuracy' + subject + '.mat', {'Acc': Accuracy})




