""""""""""""""""""""""""""""""""""""""""""
#        OpenBCI - Brainflow            #
#           For Online AAD              #
""""""""""""""""""""""""""""""""""""""""""

###### Imports #####
import librosa, warnings, random, time, os, sys, serial, logging, argparse, mne, scipy.io, math
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from pylsl import StreamInlet, resolve_stream, StreamInfo
# from OpenBCI_lsl import *5
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

# -------------------------------- SETTING ---------------------------------#

#########

subject = ''
#subject = '0726_KKM'

###########

loc = 'kist'
#loc = 'hyu'

if loc == 'kist':
    arduino = "COM8"
    path = 'C:/Users/LeeJiWon/Desktop/OpenBCI'
    cyton = 'COM7'
elif loc == 'hyu':
    arduino = "COM4"  # casing COM4 // File box COM10
    path = 'C:/Users/user/Desktop/hy-kist/OpenBCI'
    cyton = 'COM3'  # casing COM3 // File box COM14

# Connect to port of arduino
port = serial.Serial(arduino, 9600)

# Connect Cyton with Brainflow network
board, args = Brainflow_stream(cyton)

##################
original = 'L'
opposite = 'R'

# ----------------------------- Load Speech segment data ------------------------------#

# Load All speech
allspeech = np.load(path + '/AAD/Python/Allspeech.npy')
# 60 by 3840  /1-30 : Twenty by time / 31-60 Journey by time // srat : 64

stim_T = allspeech[:30, :]  # 30 by 3840   // trial by time
stim_J = allspeech[30:, :]  # 30 by 3840   // trial by time

# ----------------------------- Parameter Setting -----------------------------#

fs = 64
tmin = 0
tmax = 250
Dir = -1
reg_lambda = 10
train = 14

# Set channels number & sampling rate
eeg_channels = board.get_eeg_channels(args.board_id)
aux_channels = board.get_analog_channels(args.board_id)
srate = board.get_sampling_rate(args.board_id)

##############################################
# Set int

r_J = []
r_T = []
Acc = []
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
eeg_record = np.zeros((16, 1))
raw_data = np.zeros((16, 1))
aux_record = np.zeros((3, 1))
tri_data = np.zeros((3, 1))
EEG_all = np.array([])
AUX_all = np.array([])
answer_all = np.array([])
correct_all = np.array([])

w = 1  # To avoid repeat when detect trigger
j = 0  # Question number
tr = 0  # trial
##############################################

# ----------------------------- Make the window for Psychopy -----------------------------#

screen = visual.Window([960, 900], screen=0, pos=[600, 0], fullscr=False,
                       winType='pyglet', allowGUI=False, allowStencil=False,
                       monitor='testMonitor', color=[-1, -1, -1], blendMode='avg',
                       units='pix')

# ------------------------------------- Warm up --------------------------------------------#
event.waitKeys(keyList=['space'], clearEvents=True)
tic = time.perf_counter()
toc = time.perf_counter()

print("Warming up")
port.write(b'2')
while toc - tic < 10:  # During 20s

    input = board.get_board_data()
    eeg_data = input[eeg_channels, :]
    aux_data = input[aux_channels, :]
    print(aux_data)  # Check input Trigger

    # If Nan value is entered, restart
    if True in np.isnan(np.asarray(eeg_data)):
        print("******************")
        print("Input NAN!")
        print("******************")
        break
    key = event.getKeys()
    if key == ["escape"]:
        core.quit()

    time.sleep(1)
    toc = time.perf_counter()

print("Warming up End")

# -------------------------------------- Intro ---------------------------------------------#
file = pd.read_excel(path + "/AAD/Python/question.xlsx")
#file_2 = pd.read_excel(path + "/AAD/Python/Comments.xlsx")
file_3 = pd.read_excel(path + "/AAD/Python/prePractice.xlsx")

event.waitKeys(keyList=['space'], clearEvents=True)
Comments('intro', path, screen,original)

# ---------------------------------------- Practice ----------------------------------------#
# Comment
text = visual.TextStim(screen, text=file_3.comment[0], height=50, color=[1, 1, 1], wrapWidth=2000)
text.draw()
screen.flip()
event.waitKeys(keyList=['space'], clearEvents=True)

for p in range(0, 2):
    port.write(b'3')  # practice speech
    practice(p, path, screen)

# ==================================================================================================#
# -------------------------------------- START EXPERIMENT ------------------------------------------#
# ==================================================================================================#

Comments(tr, path, screen,original)

# ---------- Start 30 trial ----------#

# Throw data which dont need
input = board.get_board_data()

while tr < 30:  # 30

    # ----------------------------- Psychopy Window & Serial Write ------------------------------#
    if w == 1:
        Direction(tr, original, opposite, screen, port)
        port.write(b'1')
        w = 0

    # Data acquisition
    input = board.get_board_data()
    eeg_data = input[eeg_channels, :]
    aux_data = input[aux_channels, :]
    eeg_record = np.concatenate((eeg_record, -eeg_data), axis=1)
    aux_record = np.concatenate((aux_record, aux_data), axis=1)
    raw_data = np.concatenate((raw_data, -eeg_data), axis=1)
    tri_data = np.concatenate((tri_data, aux_data), axis=1)
    print(aux_data)

    # ----------------------------- Trigger detection -----------------------------#
    # per trial
    if 1 in aux_data[1, :]:  # if the trigger is entered at 12 pin, start process. (include beep sound)

        print("Input Trigger {0}".format(tr + 1))

        print("Start Speech")

        # Find onset point
        index = np.where(aux_record[1, :] != 0)  # Find Onset index
        onset = index[0][0]

        # Format per trial
        i = 0  # Window number
        work = 0  # Time count
        check = -3
        Acc = []

        # ----------------------------- Working while 60s -----------------------------#

        # onset 부터 3초 지나고 원하는 시간(한 trial) 동안 돌아가도록
        speech = onset + (srate * 3) + 1  # 3초 후 부터가 speech 이기에 +1
        while i != 46:  # 46 번째 window 까지 (0~45)

            # if the processing time exceeds 1s, no time sleep
            if work > 1:
                work = 1

            # Wait 1s
            time.sleep(1 - work)

            # predicted time count
            check = check + 1
            print("--------------------------------------")
            print("Time sleep : {0}".format(check))

            # Time count
            start = time.perf_counter()

            ### Receive sample ###
            input = board.get_board_data()

            # Separate EEG, AUX
            eeg_data = input[eeg_channels, :]
            aux_data = input[aux_channels, :]  # 11,12,13 / 0 or 1
            print("INPUT : {0}".format(len(eeg_data.T)))

            # Stack data
            eeg_record = np.concatenate((eeg_record, -eeg_data), axis=1)  # channel by time
            aux_record = np.concatenate((aux_record, aux_data), axis=1)
            raw_data = np.concatenate((raw_data, -eeg_data), axis=1)
            tri_data = np.concatenate((tri_data, aux_data), axis=1)

            # Count time
            end = time.perf_counter()
            work = end - start

            # Stack samples until 15s and window sliding per 1s
            # After 15s, rely on time sleep.
            if check >= 15:

                # Adjust data as acquired from that time.
                win = eeg_record[:, speech + srate * (i):]      # channel by time

                if len(win.T) > srate * (15):
                    win = eeg_record[:, speech + srate * (i): speech + srate * (15 + i)]        # 15 by 1875

                # Check print
                print("Window number : {0}".format(i + 1))
                print("Time Check : {0}s".format(len(eeg_record[:, speech:].T) / srate))

                # ----------------------------- Pre-processing -----------------------------#
                # preprocessing_ha.py
                # Delete No.8 channel (not use)
                win = np.delete(win, 7, axis=0)     # delete 7 row ( 8 channel )

                win = Preproccessing(win, srate, 0.5, 8, 601)  # data, sampling rate, low-cut, high-cut, filter order

                data_l = len(win.T)     # To check the length of inputted data

                # ------------------------------- Train set -------------------------------#
                if tr < train:  # int train
                    state = "Train set"

                    ## mTRF train function ##
                    model, tlag, inter = mtrf_train(stim_J[tr:tr + 1, 64 * (i):64 * (i) + data_l].T, win.T, fs, Dir,
                                                    tmin, tmax, reg_lambda)

                    'model - (16,17,1)  / tlag - (17,1) / inter - (16,1)'

                    # Sum w - window
                    if i == 0:
                        model_w = model
                        inter_w = inter
                    else:  # i > 0 - 45까지
                        model_w = np.add(model_w, model)
                        inter_w = np.add(inter_w, inter)

                    i = i + 1
                    end = time.perf_counter()

                # ------------------------------- Test set -------------------------------#
                else:
                    state = "Test set"

                    ## Calculate Predicted signal ##
                    pred, r_j, p, mse = mtrf_predict(stim_J[tr:tr + 1, 64 * (i):64 * (i) + data_l].T, win.T, model, fs,
                                                     Dir, tmin, tmax, inter)
                    pred, r_t, p, mse = mtrf_predict(stim_T[tr:tr + 1, 64 * (i):64 * (i) + data_l].T, win.T, model, fs,
                                                     Dir, tmin, tmax, inter)

                    # Stock correlation value per window(i)
                    r_J = np.append(r_J, r_j)
                    r_T = np.append(r_T, r_t)

                    ###### Estimate accuracy #####
                    if r_j > r_t:
                        acc = 1
                    else:
                        acc = 0

                    print("======= acc : {0} ".format(acc))

                    # Save acc for entire Accuracy
                    Acc = np.append(Acc, acc)

                    # Up window number
                    i = i + 1

                # -- End one window --#

                # In trial 27~30, switch direction
                if tr + 1 >= 27 and check > 25:
                    switching(tr, check, original, opposite, screen)

                # Time count
                end = time.perf_counter()
                work = end - start
                print("working time = {0}s".format(work))

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
                model_wm = model_wt / (i * (tr + 1))
                inter_wm = inter_wt / (i * (tr + 1))
                model = model_wm
                inter = inter_wm

        elif state == "Test set":
            # Stack correlation value collected during one trial
            entr_J.append(r_J)
            entr_T.append(r_T)

            r_T = []
            r_J = []

            # Collect Accuracy per trial
            ACC = np.append(ACC, np.mean(Acc))
            print("\n==================================\n")
            print("Present Accuracy = {0}%".format(ACC[-1] * 100))
            print("\n==================================\n")


        #####----- Question -----#####
        try:
            print("Question Time")

            correct, answer = Question(j, path, screen)

            Correct.append(correct)
            Answer.append(answer)
            # Correct.append(correct)
            Answer.append(answer)
            j = j + 1
        except KeyError:  # 마지막 질문 후 에러나서
            pass

        #=======  Data acquisition for rest  =======#
        input = board.get_board_data()

        # Separate EEG, AUX
        eeg_data = input[eeg_channels, :]
        aux_data = input[aux_channels, :]  # 11,12,13 / 0 or 1

        # Stack data
        eeg_record = np.concatenate((eeg_record, -eeg_data), axis=1)  # channel by time
        aux_record = np.concatenate((aux_record, aux_data), axis=1)
        raw_data = np.concatenate((raw_data, -eeg_data), axis=1)
        tri_data = np.concatenate((tri_data, aux_data), axis=1)

        #===== Stack eeg_record per trial & Save =====#
        EEG.append(eeg_record.T)
        AUX.append(aux_record.T)
        # For Next trial
        eeg_record = np.zeros((16, 1))
        aux_record = np.zeros((3, 1))

        #--------------------------------------------------------------------------------#
        # Save per trial // eeg, trigger, accuracy ,behavior
        EEG_all = np.asarray(EEG)
        AUX_all = np.asarray(AUX)
        scipy.io.savemat(path + '/save_data/E_' + subject + '.mat', {'EEG': EEG_all})
        scipy.io.savemat(path + '/save_data/A_' + subject + '.mat', {'AUX': AUX_all})
        scipy.io.savemat(path + '/save_data/RAW_' + subject + '.mat', {'RAW': raw_data})
        scipy.io.savemat(path + '/save_data/TRIGGER_' + subject + '.mat', {'TRIGGER': tri_data})
        scipy.io.savemat(path + '/save_data/Accuracy_' + subject + '.mat', {'Acc': ACC})
        correct_all = np.asarray(Correct)
        scipy.io.savemat(path + '/save_data/Behavior_' + subject + '.mat', {'Behavior': correct_all})

        tr = tr + 1
        w = 1

        # ------------------ comment about session ---------------------#
        Comments(tr, path, screen, original)

# ----------------------------- 30 trial End -----------------------------#

# Represent Total accuracy
print("\n===================================\n")
print("=== Total Accuracy = {0}% ===".format(mean(ACC)*100))
print("\n===================================\n")

# END
print("The End")
final = visual.TextStim(screen, text="실험이 끝났습니다. \n\n 수고하셨습니다.", height=50, color=[1, 1, 1], wrapWidth=2000)
final.draw()
screen.flip()
time.sleep(3)

port.close()
screen.close()
board.stop_stream()
board.release_session()

#### save ####
# mat save
answer_all = np.asarray(Answer)
scipy.io.savemat(path + '/save_data/Answer_' + subject + '.mat', {'Answer': answer_all})

# np save
np.save(path + '/save_data/EEG_' + subject, EEG_all)
np.save(path + '/save_data/A_' + subject, AUX_all)
np.save(path + '/save_data/RAW_' + subject, raw_data)
np.save(path + '/save_data/EEG_' + subject, tri_data)
np.save(path + '/save_data/All_Accuracy_' + subject, ACC)
entr_J = np.asarray(entr_J)
entr_T = np.asarray(entr_T)
np.save(path + '/save_data/All_correlation_att_' + subject, entr_J)
np.save(path + '/save_data/All_correlation_unatt_' + subject, entr_T)