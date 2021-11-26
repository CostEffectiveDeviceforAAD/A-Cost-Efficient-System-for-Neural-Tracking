""""""""""""""""""""""""""""""""""""""""""
#        OpenBCI - Brainflow            #
#           For Online AAD              #
""""""""""""""""""""""""""""""""""""""""""

###### Imports #####
import librosa, warnings, random, time, os, sys, serial, logging, argparse, mne, scipy.io, math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal
from scipy.signal import butter, lfilter, resample, filtfilt
# from helper import *
from pymtrf import *
from psychopy import visual, core, event
from preprocessing_ha import *
from Comments import *
from Direction import *
from brainflow.board_shim import BoardShim, BrainFlowInputParams, BoardIds, BrainFlowError, LogLevels
from brainflow.data_filter import DataFilter, FilterTypes, AggOperations, WindowFunctions, DetrendOperations
from Brainflow_stream import *

#################################################################################################
# ---------------------------------- Experimental SETTING --------------------------------------#
#################################################################################################

# set info
subject = '001'             # Subject number
original = 'R'              # First attention direction
opposite = 'L'              # opposite direction
arduino = 'COM8'            # Arduino serial port number
cyton = 'COM7'              # OpenBCI board Bluetooth port number
path = 'C:/Users/'          # Base Directory

# Connect to port of arduino
port = serial.Serial(arduino, 9600)
# Connect Cyton with Brainflow network
board, args = Brainflow_stream(cyton)
# Set channels number & sampling rate
eeg_channels = board.get_eeg_channels(args.board_id)
aux_channels = board.get_analog_channels(args.board_id)
srate = board.get_sampling_rate(args.board_id)

# ----------------------------------- Load Speech segment data -----------------------------------#
# Load Stimulus data
allspeech = np.load(path + '/Allspeech.npy')
# 60 by 3840 / trial by time / 1-30 : Twenty / 31-60 Journey // sampling rate : 64 Hz

stim_T = allspeech[:30, :]
stim_J = allspeech[30:, :]

# ---------------------------------------- Parameter Setting -------------------------------------#
tmin = 0
tmax = 250          # Time-lag
Dir = -1            # Backward
reg_lambda = 10     # Lambda value
fs = 64             # Hope Sampling rate
train = 14          # The number of train trial




# Set int
r_J = []
r_T = []

ACC = []
model_w = []
inter_w = []
entr_J = []
entr_T = []
EEG = []
AUX = []
Correct = []
Answer = []
start = []
end = []

# for trial array / for entire array
EEG_all = np.array([])
AUX_all = np.array([])
answer_all = np.array([])
correct_all = np.array([])



w = 1  # To avoid repeat when detect trigger
j = 0  # Question number

##############################################



# Make the window for visual presentation
screen = visual.Window([960, 900], screen=1, pos=[600, 0], fullscr=True,
                       winType='pyglet', allowGUI=False, allowStencil=False,
                       monitor='testMonitor', color=[-1, -1, -1], blendMode='avg',
                       units='pix')

# ------------------------------------------ Warm up -------------------------------------------------#
# Warm up the device before experiment
event.waitKeys(keyList=['space'], clearEvents=True)
print("Warming up")

tic = time.perf_counter()
toc = time.perf_counter()

port.write(b'2')
while toc - tic < 10:  # During 10s

    input = board.get_board_data()
    print(input[aux_channels, :])  # Check input Trigger

    time.sleep(1)
    toc = time.perf_counter()

print("Warming up Done")
# If you have some problem with experiment, do experiment again.


# ------------------------------------------- Intro -------------------------------------------------#
# Load Intro command file
file = pd.read_excel(path + "/AAD/Python/question.xlsx")
file_3 = pd.read_excel(path + "/AAD/Python/prePractice.xlsx")

event.waitKeys(keyList=['space'], clearEvents=True)
Comments('intro', path, screen, original)


# -------------------------------------------- Practice --------------------------------------------#
# Presentation Command for practice before experiment
text = visual.TextStim(screen, text=file_3.comment[0], height=50, color=[1, 1, 1], wrapWidth=2000)
text.draw()
screen.flip()
event.waitKeys(keyList=['space'], clearEvents=True)

for p in range(0, 2):
    port.write(b'3')            # For practice speech
    practice(p, path, screen)


# ==================================================================================================#
# -------------------------------------- START EXPERIMENT ------------------------------------------#
# ==================================================================================================#
RAWData = np.zeros((16, 1))    # EEG 16 channel by 1
AUXData = np.zeros((3, 1))     # Trigger 3 channel by 1
tr = 0
w = 1
# Comment before first session.
Comments(tr, path, screen, original)

# ----------------------------------------- Start 30 trial ------------------------------------------#

# Throw data which don't need
input = board.get_board_data()

while tr < 30:  # 30

    if w == 1:
        Direction(tr, original, opposite, screen, port)
        port.write(b'1')        # Signal for trial onset to arduino
        w = 0

    # For Next trial, reset the data format
    eeg_record = np.zeros((16, 1))
    aux_record = np.zeros((3, 1))

    # Data acquisition
    input = board.get_board_data()
    eeg_record = np.concatenate((eeg_record, -input[eeg_channels, :]), axis=1)
    # If you connected with cathode electrode in openBCI board, you should transmit cathode of EEG data
    aux_record = np.concatenate((aux_record, input[aux_channels, :]), axis=1)


    # ----------------------------- Trigger detection -----------------------------#

    if 1 in input[aux_channels, :][1, :]:  # if the trigger is detected at 12 pin, start process. (include beep sound)
        print("Input Trigger {0}".format(tr + 1))
        print("Start Speech")

        # Find onset point
        index = np.where(aux_record[1, :] != 0)  # Find Onset index
        onset = index[0][0]
        # Format per trial
        i = 0       # Window number
        work = 0    # Time count
        check = -3  # attention cue sound 3 second
        Acc = []

        # ----------------------------- Working while 60s = one trial-----------------------------#

        # Find speech onset point exclude attention cue sound
        speech = onset + (srate * 3) + 1
        while i != 46:          # During 46 windows

            # If the processing time exceeds 1s, no time sleep
            if work > 1:
                work = 1
            # Wait 1 second to be stacked EEG data and update EEG data per 1 second
            time.sleep(1 - work)

            # Visualize Time
            check = check + 1
            print("Time sleep : {0}".format(check))
            # Time count
            start = time.perf_counter()

            ### Receive sample ###
            input = board.get_board_data()
            # Stack data
            eeg_record = np.concatenate((eeg_record, -input[eeg_channels, :]), axis=1)  # channel by time
            aux_record = np.concatenate((aux_record, input[aux_channels, :]), axis=1)

            # Count time
            end = time.perf_counter()
            work = end - start

            # Stack data until 15s and window sliding per 1s
            # After 15s, rely on time sleep.
            if check >= 15:

                # Adjust data as acquired from that time.
                win = eeg_record[:, speech + srate * (i):]      # channel by time

                if len(win.T) > srate * (15):       # when exceed long of 15 second
                    win = eeg_record[:, speech + srate * (i): speech + srate * (15 + i)]        # 15 by 1875

                # Check print
                print("Window number : {0}".format(i + 1))
                print("Time Check : {0}s".format(check))

                # ----------------------------- Pre-processing -----------------------------#
                win = np.delete(win, 7, axis=0)                 # delete 7 row ( 8 channel/Fp1 )
                win = Preproccessing(win, srate, 0.5, 8, 601)   # data, sampling rate, low-cut, high-cut, filter order
                data_l = len(win.T)                             # To check the length of inputted data

                # ============================== Train set ==================================#
                if tr < train:  # int train
                    state = "Train set"
                    # Train decode model
                    model, tlag, inter = mtrf_train(stim_J[tr:tr + 1, 64 * (i):64 * (i) + data_l].T, win.T, fs, Dir,
                                                            tmin, tmax, reg_lambda)

                    # Sum w - window
                    if i == 0:
                        model_w = model
                        inter_w = inter
                    else:  # i > 0 - 45까지
                        model_w = np.add(model_w, model)
                        inter_w = np.add(inter_w, inter)
                    i = i + 1
                    end = time.perf_counter()

                # ============================== Test set ===================================#
                else:
                    state = "Test set"

                    # Reconstruct speech
                    pred, corr_j, p, mse = mtrf_predict(stim_J[tr:tr + 1, 64 * (i):64 * (i) + data_l].T, win.T, model, fs,
                                                     Dir, tmin, tmax, inter)
                    pred, corr_t, p, mse = mtrf_predict(stim_T[tr:tr + 1, 64 * (i):64 * (i) + data_l].T, win.T, model, fs,
                                                     Dir, tmin, tmax, inter)

                    # Stock correlation value per window(i)
                    corr_J = np.append(corr_J, corr_j)
                    corr_T = np.append(corr_T, corr_t)

                    # Compare with both correlation values
                    if corr_j > corr_t:
                        acc = 1
                    else:
                        acc = 0
                    print("======= acc : {0} =======".format(acc))

                    # Save acc for entire Accuracy
                    Acc = np.append(Acc, acc)
                    i = i + 1

                # ===== End one window =====#

                # In trial 27~30, switch direction
                if tr + 1 >= 27 and check > 25:
                    switching(tr, check, original, opposite, screen)

                # Time count
                end = time.perf_counter()
                work = end - start
                print("working time = {0}s".format(work))

        # ------------------------ End 60s - one trial ------------------------#
        # Calculate per trial
        if state == "Train set":
            # Sum decoder model to average
            if tr == 0:
                model_wt = model_w
                inter_wt = inter_w
            elif tr > 0:
                model_wt = np.add(model_wt, model_w)
                inter_wt = np.add(inter_wt, inter_w)

            # Average at last train trial
            if tr == 13:
                model = model_wt / (i * (tr + 1))
                inter = inter_wt / (i * (tr + 1))

        elif state == "Test set":
            # Stack correlation value collected during one trial
            Allcorr_J.append(corr_J)
            Allcorr_T.append(corr_T)
            corr_T = []
            corr_J = []

            # Collect Accuracy per trial
            ACC = np.append(ACC, np.mean(Acc))
            print("\n==================================\n")
            print("Present Accuracy = {0}%".format(ACC[-1] * 100))
            print("\n==================================\n")

        # --------------------------- Questions --------------------------- #
        try:
            print("Question Time")
            correct, answer = Question(j, path, screen)
            Correct.append(correct)
            Answer.append(answer)
            j = j + 1
        except KeyError:  # for error of last question
            pass

        #=======  Data acquisition for rest  =======#
        input = board.get_board_data()
        eeg_record = np.concatenate((eeg_record, -input[eeg_channels, :]), axis=1)  # channel by time
        aux_record = np.concatenate((aux_record, input[aux_channels, :]), axis=1)

        #===== Stack eeg_record per trial & Save =====#
        RAWData = np.concatenate((RAWData, eeg_record), axis=1)
        AUXData = np.concatenate((AUXData, aux_record), axis=1)

        # =================================== SAVE DATA ===================================== #
        # Save per trial - RAW data, AUX data, accuracy ,behavior
        # numpy file
        np.save(path + '/save_data/RAW_' + subject, RAWData)
        np.save(path + '/save_data/AUX_' + subject, AUXData)

        # Mat file
        scipy.io.savemat(path + '/save_data/RAW_' + subject + '.mat', {'RAW': RAWData})
        scipy.io.savemat(path + '/save_data/TRIGGER_' + subject + '.mat', {'TRIGGER': AUXData})
        scipy.io.savemat(path + '/save_data/Accuracy_' + subject + '.mat', {'Acc': ACC})
        correct_all = np.asarray(Correct)
        scipy.io.savemat(path + '/save_data/Behavior_' + subject + '.mat', {'Behavior': correct_all})

        # Format current trial
        tr = tr + 1
        w = 1

        # ------------------ comment about session ---------------------#
        Comments(tr, path, screen, original)


#################################################################################################
#                                       EXPERIMENT DONE                                         #
#################################################################################################

# END
print("The End")
final = visual.TextStim(screen, text="실험이 끝났습니다. \n\n 수고하셨습니다.", height=50, color=[1, 1, 1], wrapWidth=2000)
final.draw()
screen.flip()
time.sleep(3)

# Total accuracy
print("\n===================================\n")
print("=== Total Accuracy = {0}% ===".format(mean(ACC)*100))
print("\n===================================\n")

port.close()
screen.close()
board.stop_stream()
board.release_session()

# Save Detail
# mat save
answer_all = np.asarray(Answer)
scipy.io.savemat(path + '/save_data/Answer_' + subject + '.mat', {'Answer': answer_all})

# np save
entr_J = np.asarray(entr_J)
entr_T = np.asarray(entr_T)
np.save(path + '/save_data/All_correlation_att_' + subject, entr_J)
np.save(path + '/save_data/All_correlation_unatt_' + subject, entr_T)
scipy.io.savemat(path + '/save_data/Correlation_att_' + subject + '.mat', {'corr_att': entr_J})
scipy.io.savemat(path + '/save_data/Correlation_unatt_' + subject + '.mat', {'corr_unatt': entr_T})