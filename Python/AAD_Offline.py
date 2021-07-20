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
from preprocessing_ha import *
from Comments import *
from Direction import *


# -------------------------------- SETTING ---------------------------------#

#########
subject = ''
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

# Load data

raw_mat = io.loadmat(path + '/0714_LJW/RAW_01LJW.mat')
raw = raw_mat['RAW']        # channel by time
raw = np.concatenate((raw, np.ones([16,100])), axis=1)  # for final trial (lack of time)
tri_mat = io.loadmat(path + '/0714_LJW/TRIGGER_01LJW.mat')
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

#for ch in [5,11]:

while tr < 30:  # 30

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

        win = Preproccessing(win, srate, 0.5, 8, 601)  # data, sampling rate, low-cut, high-cut, filter order

        # Select channel
        win = win[[0,1,2,3,4,5,8,9,10,11,12,13,14,15],:]

        # ------------------------------- Train set -------------------------------#
        if tr < train:  # int train
            state = "Train set"

            ## mTRF train function ##
            model, tlag, inter = mtrf_train(stim_L[tr:tr + 1, 64 * (i) : 64 * (15 + i)].T, win.T, fs, Dir,
                                            tmin, tmax, reg_lambda)

            'model - (16,17,1)  / tlag - (17,1) / inter - (16,1)'

            #================================================================

            pred_l, r_l, p, mse = mtrf_predict(stim_L[tr:tr+1, 64 * (i) : 64 * (15 + i)].T, win.T, model, fs,
                                             Dir, tmin, tmax, inter)
            pred_r, r_r, p, mse = mtrf_predict(stim_R[tr:tr+1, 64 * (i) : 64 * (15 + i)].T, win.T, model, fs,
                                             Dir, tmin, tmax, inter)

            predic_l.append(pred_l)
            predic_r.append(pred_r)

            print("Train")
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

            ## Calculate Predicted signal ##
            pred_l, r_l, p, mse = mtrf_predict(stim_L[tr:tr+1, 64 * (i) : 64 * (15 + i)].T, win.T, model, fs,
                                             Dir, tmin, tmax, inter)
            pred_r, r_r, p, mse = mtrf_predict(stim_R[tr:tr+1, 64 * (i) : 64 * (15 + i)].T, win.T, model, fs,
                                             Dir, tmin, tmax, inter)

            predic_l.append(pred_l)
            predic_r.append(pred_r)

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
        print("Present Channel = {0}".format(ch))
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
        print("Present Channel = {0}".format(ch))
        print("\n==================================\n")

    print(tr)
    tr = tr + 1
    w = 1


# ----------------------------- 30 trial End -----------------------------#
######################################################################################
for tr in range(0,train):
    Acc = []
    speech = onset[tr] + (srate * 3) + 1

    #raw = eeg[:60, :, tr]

    for i in range(0,46):

        win = raw[:, speech + srate * (i): speech + srate * (15 + i)]
        win = Preproccessing(win, srate, 0.5, 8, 601)

        # Select channel
        win = win[[0,1,2,3,4,5,8,9,10,11,12,13,14,15],:]  #exclude o12, fp1

        ## Calculate Predicted signal ##
        pred_l, r_l, p, mse = mtrf_predict(stim_L[tr:tr + 1, 64 * (i): 64 * (15 + i)].T, win.T, model, fs,
                                         Dir, tmin, tmax, inter)
        pred_r, r_r, p, mse = mtrf_predict(stim_R[tr:tr + 1, 64 * (i): 64 * (15 + i)].T, win.T, model, fs,
                                         Dir, tmin, tmax, inter)

        predic_l.append(pred_l)
        predic_r.append(pred_r)

        print("Train check")
        # Stock correlation value per window(i)
        r_L = np.append(r_L, r_l)
        r_R = np.append(r_R, r_r)

        ###### Estimate accuracy #####
        if r_l > r_r:
            acc = 1
        else:
            acc = 0

        print("======= acc : {0} ".format(acc))

        # Save acc for entire Accuracy
        Acc = np.append(Acc, acc)

        # Up window number
        i = i + 1

    ACC = np.append(ACC, np.mean(Acc))
    Pre_L.append(predic_l)
    Pre_R.append(predic_r)
    print("\n==================================\n")
    print("Present Accuracy = {0}%".format(ACC[-1] * 100))
    print("Present Channel = {0}".format(ch))
    print("\n==================================\n")


    tr = tr + 1

print("Total Accuracy = {0}%".format(mean(ACC[14:29])*100))
Accuracy.append(ACC)
#ACC = []
#tr = 0


print("The End")

#### save ####
# mat save
# Save per trial // eeg, trigger, accuracy ,behavior
Accuracy = np.asarray(Accuracy)
scipy.io.savemat(path + '/save_data/Accuracy' + subject + '.mat', {'Acc': Accuracy})
scipy.io.savemat(path + '/save_data/Predict_L' + subject + '.mat', {'Pre_L': Pre_L})
scipy.io.savemat(path + '/save_data/Predict_R' + subject + '.mat', {'Pre_R': Pre_R})




