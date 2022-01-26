""""""""""""""""""""""""""""""""""""""""""
#        EMA-applied Accuracy           #
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

#################################################################################################
#---------------------------------- Experimental SETTING ---------------------------------------#
#################################################################################################
# Load data / mat file
# set info
path = 'C:/Users/LeeJiWon/Desktop/OpenBCI'          # Base Directory
sub = '5'

Corr_j = np.load(path + '/Recording data/'+ sub + '/All_correlation_att_' + sub + '.npy')
Corr_t = np.load(path + '/Recording data/'+ sub + '/All_correlation_unatt_' + sub + '.npy')      # 16 by 46
EmaCorr_j = np.zeros((16,46))
EmaCorr_t = np.zeros((16,46))
acc = np.zeros((16,46))
ACC = []

# Load All speech
allspeech = np.load(path + '/AAD/Python/Allspeech.npy')
# 60 by 3840  /1-30 : left by time / 31-60 right by time // srat : 64

stim_T = allspeech[:30, :]  # 30 by 3840  Twenty   // trial by time
stim_J = allspeech[30:, :]  # 30 by 3840  Journey  // trial by time

# ----------------------------- Parameter Setting -----------------------------#

for tr in range(0,16):
    for i in range(0,46):
        # Exponential Moving Average
        if i == 0:
            EmaCorr_j[tr, i] = Corr_j[tr,i]
            EmaCorr_t[tr , i] = Corr_t[tr,i]
        else:
            weight = 2 / (i + 2)
            EmaCorr_j[tr, i] = np.add((weight * Corr_j[tr,i]), ((1 - weight) * EmaCorr_j[tr, i - 1]))
            EmaCorr_t[tr, i] = np.add((weight * Corr_t[tr,i]), ((1 - weight) * EmaCorr_t[tr, i - 1]))

        ######  Estimate accuracy  #####
        if EmaCorr_j[tr, i] > EmaCorr_t[tr, i]:
            acc[tr, i] = 1
        else:
            acc[tr, i] = 0

    ACC = np.append(ACC, mean(acc[tr:tr+1,:]))

# ----------------------------- 30 trial End -----------------------------#
print("Total Accuracy = {0}%".format(mean(ACC)*100))
print("Fix Accuracy = {0}%".format(mean(ACC[:12]*100)))
print("Switch Accuracy = {0}%".format(mean(ACC[12:]*100)))
print("The End")

#### save ####
np.save(path + '/save_data/Accuracy_EMA_' + sub, ACC)
np.save(path + '/save_data/AllcorrJ_EMA_' + sub, EmaCorr_j)
np.save(path + '/save_data/AllcorrT_EMA_' + sub, EmaCorr_t)



