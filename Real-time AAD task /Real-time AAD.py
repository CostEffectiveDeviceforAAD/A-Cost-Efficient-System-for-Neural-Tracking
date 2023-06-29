""""""""""""""""""""""""""""""""""""""""""
#        A cost-effective system         #
#         For rea-time AAD task          #
""""""""""""""""""""""""""""""""""""""""""

###### Imports #####
import librosa, warnings, time, serial, scipy.io
import numpy as np
import pandas as pd
from pymtrf import *
from PreProcessing import *
from Brainflow_stream import *
from EMA import *

#################################################################################################
#---------------------------------- Experimental SETTING ---------------------------------------#
#################################################################################################
# set info
subject = '001'             # Subject number
arduino = 'COM3'            # Arduino serial port number
cyton = 'COM11'             # OpenBCI board Bluetooth port number
path = 'C:/Users/'          # Base Directory

# Connection
port = serial.Serial(arduino, 9600)     # Connect to port of arduino
board, args = Brainflow_stream(cyton)   # Connect to Cyton with Brainflow network
# Set channels number & sampling rate
eeg_channels = board.get_eeg_channels(args.board_id)
aux_channels = board.get_analog_channels(args.board_id)
srate = board.get_sampling_rate(args.board_id)
# ----------------------------------- Load Speech segment data -----------------------------------#
# Load Stimulus data
allspeech = np.load(path + '/YourStimuli.npy')
stim_att = allspeech[:] # attended speech
stim_utt = allspeech[:] # unattended speech
# ---------------------------------------- Parameter Setting -------------------------------------#
tmin = 0
tmax = 250          # Time-lag
Dir = -1            # Backward
reg_lambda = 10     # Lambda value
fs = 64             # post sampling rate
timelag = len(range(int(floor(tmin/1000*fs)), int(ceil(tmax/1000*fs)+1)))

# Ex_int
Ttrial = 30
Train_tr = 14
Winwidth = 15
window_Tnum = 60-(width-1)
ChLen = len(eeg_channels)

# Array
#corr_att = np.zeros((ChLen,window_Tnum))
#corr_utt = np.zeros((ChLen,window_Tnum))
EmaCorr_att = np.zeros((ChLen,window_Tnum))
EmaCorr_utt = np.zeros((ChLen,window_Tnum))
acc = np.zeros((ChLen, window_Tnum))
start = []
end = []

# ==================================================================================================#
# -------------------------------------- START EXPERIMENT ------------------------------------------#
# ==================================================================================================#
# ----------------------------------------- Start 30 trial -----------------------------------------#
tr = 0  # trial number
z = 1   # To avoid repeat when detect trigger
# Throw data which don't need
input = board.get_board_data()
while tr < Ttrial:
    if z == 1:
        port.write(b'1')        # To sent Signal for trial start to arduino
        z = 0

    # For Next trial, reset the data format
    eeg_record = np.zeros((ChLen, 0))
    aux_record = np.zeros((3, 0))

    # Data acquisition
    input = board.get_board_data()
    eeg_record = np.concatenate((eeg_record, -input[eeg_channels, :]), axis=1)   #If you connected with cathode electrode on openBCI board, you should transmit cathode of EEG data
    aux_record = np.concatenate((aux_record, input[aux_channels, :]), axis=1)

    # ----------------------------- Trigger detection -----------------------------#
    if 1 in input[aux_channels, :][1, :]:   # if the trigger is detected at 12 pin, start process. (include beep sound)
        print("Input Trigger {0}".format(tr + 1))
        print("Start Sound")
        # Find onset point
        index = np.where(aux_record[1, :] != 0)
        onset = index[0][0]
        # Format per trial
        i = 0        # Window number
        work = 0     # Time count
        check = -3   # attention cue sound 3 second

        # ----------------------------- Working while 60s = one trial-----------------------------#
        # Find speech onset point exclude attention cue sound
        speech = onset + (srate * 3) + 1
        while i != window_Tnum:          # During 46 windows

            # If the processing time exceeds 1s, no time sleep
            if work > 1:
                work = 1
            # Wait 1 second to be stacked EEG data and update EEG data per 1 second
            time.sleep(1 - work)
            # Visualize Time
            check = check + 1
            # Time counting
            start = time.perf_counter()

            # Acquire data
            input = board.get_board_data()
            eeg_record = np.concatenate((eeg_record, -input[eeg_channels, :]), axis=1)  # channel by time
            aux_record = np.concatenate((aux_record, input[aux_channels, :]), axis=1)

            # Work time
            end = time.perf_counter()
            work = end - start
            # Stack data until 15s and window sliding per 1s
            # Go to next step after 15s.
            if check >= Winwidth:
                # Adjust data as acquired from that time.
                win = eeg_record[:, speech + srate * (i):]      # channel by time

                if len(win.T) > srate * (Winwidth):   # when exceed minimum length of window
                    win = eeg_record[:, speech + srate * (i): speech + srate * (Winwidth + i)]      # 16 by 1875

                # ----------------------------- Pre-processing -----------------------------#
                win = np.delete(win, 7, axis=0)                 # delete 7 row ( 8 channel/Fp1 )
                win = Preproccessing(win, srate, 0.5, 8, 601)   # data, sampling rate, low-cut, high-cut, filter order
                data_l = len(win.T)                             # To check the length of inputted data
                # ============================== Train set ==================================#
                if tr < Train_tr:  # int train
                    # Train decode model
                    model, tlag, inter = mtrf_train(stim_att[tr:tr + 1, fs * (i):fs * (i) + data_l].T, win.T, fs, Dir,
                                                            tmin, tmax, reg_lambda)
                    model_w += model
                    inter_w += inter

                    # Average at last train trial
                    if tr == Train_tr-1:
                        Train_model = model_w / (window_Tnum * Train_tr)  # number of window * of trial
                        Train_inter = inter_w / (window_Tnum * Train_tr)

                # ============================== Test set ===================================#
                else:
                    # Reconstruct speech
                    pred, corr_att, p, mse = mtrf_predict(stim_att[tr:tr + 1, fs * (i):fs * (i) + data_l].T, win.T, Train_model, fs,
                                                     Dir, tmin, tmax, Train_inter)
                    pred, corr_utt, p, mse = mtrf_predict(stim_utt[tr:tr + 1, fs * (i):fs * (i) + data_l].T, win.T, Train_model, fs,
                                                     Dir, tmin, tmax, Train_inter)
                    ''' # without EMA
                    # Compare with both correlation values
                    if corr_att > corr_utt:
                        acc[tr-14,i] = 1
                    else:
                        acc[tr-14,i] = 0

                    # Stock correlation value per window(i)
                    corr_att[tr-14,i] = np.array(corr_att)
                    corr_utt[tr-14,i] = np.array(corr_utt)
                    '''
                    # Exponential Moving Average
                    EmaCorr_att, EmaCorr_utt = EMA(corr_att, corr_utt, EmaCorr_att, EmaCorr_utt, i, tr)

                    # Compare with both correlation values
                    if EmaCorr_att[tr - Train_tr, i] > EmaCorr_utt[tr - Train_tr, i]:
                        acc[tr-Train_tr, i] = 1
                    else:
                        acc[tr-Train_tr, i] = 0

                # Plus window number
                i += 1

                # Time count
                end = time.perf_counter()
                work = end - start
        # ------------------------ End one trial ------------------------#
        # Format current trial
        tr += 1
        z = 1

#################################################################################################
#                                       EXPERIMENT FINISH                                       #
#################################################################################################
port.close()
screen.close()
board.stop_stream()
board.release_session()
