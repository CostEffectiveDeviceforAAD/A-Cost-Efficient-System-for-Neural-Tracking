
""""""""""""""""""""""""""""""""""""""""""
 #        OpenBCI - Brainflow           #
 #           For Online AAD              #
""""""""""""""""""""""""""""""""""""""""""


###### Imports #####
import librosa, warnings, random, time, os, sys, serial, logging, argparse, mne, scipy.io, math, datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal
from scipy.signal import butter, lfilter, resample, filtfilt
#from helper import *
from pymtrf import *
from psychopy import visual, core, event
from preprocessing_ha import *
from brainflow.board_shim import BoardShim, BrainFlowInputParams, BoardIds, BrainFlowError, LogLevels
from brainflow.data_filter import DataFilter, FilterTypes, AggOperations, WindowFunctions, DetrendOperations
from Brainflow_stream import *



#-------------------------------- SETTING ---------------------------------#
#loc = 'kist'
loc = 'hyu'

if loc == 'kist':
    arduino = "COM8"
    path = 'C:/Users/LeeJiWon/Desktop/OpenBCI'
    cyton = 'COM7'
elif loc == 'hyu':
    arduino = "COM10"
    path = 'C:/Users/user/Desktop/hy-kist/OpenBCI'
    cyton = 'COM15'

# Connect to port of arduino
port = serial.Serial(arduino, 9600)

# Connect Cyton with Brainflow network
board, args = Brainflow_stream(cyton)       # kist : COM7 / hy: COM15

#----------------------------- Load Speech segment data ------------------------------#

# Load All speech
allspeech = np.load(path + '/AAD/Python/Allspeech.npy')
# 60 by 3840  /1-30 : left by time / 31-60 righy by time // srat : 64

stim_L = allspeech[:30, :]       # 30 by 3840   // trial by time
stim_R = allspeech[30:, :]       # 30 by 3840   // trial by time

#------------------------------------ Question ------------------------------------------------#

file = pd.read_excel(path + "/AAD/Python/question.xlsx")

def Question(j, file):

    try :
        # Question 1
        text3 = visual.TextStim(screen, text = file.tweenty_Q1[j], height=50, color=[1, 1, 1], wrapWidth=2000)
        text3.draw()
        screen.flip()

        key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

        answer.append(key)
        if file.tweenty_A1[j] == int(key[0]):
            correct.append("T")
        else:
            correct.append("F")

        # Question 2
        text3 = visual.TextStim(screen, text = file.tweenty_Q2[j], height=50, color=[1, 1, 1], wrapWidth=2000)
        text3.draw()
        screen.flip()

        key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

        answer.append(key)
        if file.tweenty_A2[j] == int(key[0]):
            correct.append("T")
        else:
            correct.append("F")

        # Question 3
        text3 = visual.TextStim(screen, text = file.journey_Q1[j], height=50, color=[1, 1, 1], wrapWidth=2000)
        text3.draw()
        screen.flip()

        key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

        answer.append(key)
        if file.journey_A1[j] == int(key[0]):
            correct.append("T")
        else:
            correct.append("F")

        # Question 4
        text3 = visual.TextStim(screen, text = file.journey_Q2[j], height=50, color=[1, 1, 1], wrapWidth=2000)
        text3.draw()
        screen.flip()

        key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

        answer.append(key)
        if file.journey_A2[j] == int(key[0]):
            correct.append("T")
        else:
            correct.append("F")

    except:
        correct.append("N")

    return correct, answer

#----------------------------- Parameter Setting -----------------------------#

############## Exp parameter ################

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
r_L = r_R = Acc = ACC = []
model_w = inter_w = entr_L = entr_R = []
EEG = AUX = Correct = Answer = []
start = end = []

# for trial array / for entire array
eeg_record = raw_data= np.zeros((16,1))
aux_record = tri_data = np.zeros((3,1))
EEG_all = AUX_all = np.array([])
answer_all = correct_all = np.array([])

w = 1           # To avoid repeat when detect trigger
j = 0           # Question number
tr = 0          # trial

#----------------------------- Make the window for Psychopy -----------------------------#

screen = visual.Window([960, 900], screen = 0, pos = [600,0], fullscr = False,
                       winType = 'pyglet', allowGUI = False, allowStencil = False,
                       monitor ='testMonitor', color = [-1,-1,-1], blendMode = 'avg',
                       units = 'pix')

text_p  = visual.TextStim(screen, text = " + ", height = 100, color = [1, 1, 1])
text_d = visual.TextStim(screen, text = ">>>>", height = 100, color = [1, 1, 1])

# Draw window
screen.flip()

#--------------------------------------- Warm up ------------------------------------------#
tic = time.perf_counter()
toc = time.perf_counter()

print("Warming up")
port.write(b'2')
while toc-tic < 30:         # During 60s

    input = board.get_board_data()
    eeg_data = input[eeg_channels, :]
    aux_data = input[aux_channels, :]
    print(aux_data)

    # If Nan value is entered, restart
    if True in np.isnan(eeg_data):
        print("Input NAN")
        break
    key = event.getKeys()
    if key == ["escape"]:
        core.quit()

    time.sleep(1)
    toc = time.perf_counter()

print("Warming up End")
#-------------------------------------- Intro ---------------------------------------------#
file_2 = pd.read_excel(path + "/AAD/Python/intro.xlsx")

for i in range(0,9):

    text = visual.TextStim(screen, text=file_2.coment[i], height=50, color=[1, 1, 1],wrapWidth=2000)

    key = event.waitKeys(keyList=["space", "escape"], clearEvents=True)
    if key == ["escape"]:
        core.quit()

    text.draw()
    screen.flip()

#==================================================================================================#
#-------------------------------------- START EXPERIMENT ------------------------------------------#
#==================================================================================================#


######      Start 30 trial      ######

while tr < 30:   # 30

    #----------------------------- Psychopy Window & Serial Write ------------------------------#

    # Press Button for start
    key2 = event.getKeys()
    if key2 == ["space"] and tr == 0:

        time.sleep(3)                       # First Interval
        port.write(b'1')                    # Send signal to arduino for start sound
        input = board.get_board_data()      # Do not stack previous data
        raw_data = np.concatenate((raw_data, input[eeg_channels, :]), axis=1)
        tri_data = np.concatenate((tri_data, input[aux_channels, :]), axis=1)

        # Draw Text of direction
        text_d.draw()
        screen.flip()

    elif tr > 0 and w == 1 :                # After First trial, continues

        # Draw Text
        text_p.draw()
        screen.flip()
        time.sleep(3)                       # Interval
        port.write(b'1')                    # Send signal to arduino for start sound

        input = board.get_board_data()      # Do not stack previous data
        raw_data = np.concatenate((raw_data, input[eeg_channels, :]), axis=1)
        tri_data = np.concatenate((tri_data, input[aux_channels, :]), axis=1)
        # Draw Text of direction
        text_d.draw()
        screen.flip()
        # Avoid repeat
        w = 0

    # Trigger detection
    input = board.get_board_data()
    aux_data = input[aux_channels, :]
    eeg_record = np.concatenate((eeg_record, input[eeg_channels, :]), axis=1)
    aux_record = np.concatenate((aux_record, input[aux_channels, :]), axis=1)
    raw_data = np.concatenate((raw_data, input[eeg_channels, :]), axis=1)
    tri_data = np.concatenate((tri_data, input[aux_channels, :]), axis=1)
    print(aux_data)

    #----------------------------- Trigger detection -----------------------------#
    # per trial
    if 1 in aux_data[1,:]:      # if the trigger is entered at 12 pin, start process. (include beep sound)

        print("Input Trigger {0}".format(tr+1))

        print("Start Speech")

        # Find onset point
        index = np.where(aux_record[1,:] != 0)     # Find Onset index
        onset = index[0][0]

        # Format per trial
        i = 0           # Window number
        work = 0        # Time count


    #----------------------------- Working while 60s -----------------------------#

        # onset 부터 3초 지나고 원하는 시간(한 trial) 동안 돌아가도록
        speech = onset + (srate*3) + 1   # 3초 후 부터가 speech 이기에 +1
        while i != 46:                   # 46 번째 window 까지 (0~45)

            # if the processing time exceeds 1s, no time sleep
            if work > 1:
                work = 1

            time.sleep(1 - work)            # Wait 1s
            start = time.perf_counter()     # Time count

            #####  Receive sample  #####
            input = board.get_board_data()
            print(len(eeg_data.T))
            # Stack data
            eeg_record = np.concatenate((eeg_record, input[eeg_channels, :]), axis=1)     # channel by time
            aux_record = np.concatenate((aux_record, input[aux_channels, :]), axis=1)
            raw_data = np.concatenate((raw_data, input[eeg_channels, :]), axis=1)
            tri_data = np.concatenate((tri_data, input[aux_channels, :]), axis=1)

            # Count time
            end = time.perf_counter()
            work = end - start

            # Stack samples until 15s and window sliding per 1s
            # After 15s, rely on time sleep.
            if len(eeg_record.T) >= (speech + (srate * 15)):

                # Adjust data as acquired from that time
                win = eeg_record[:, speech + srate*(i) :]

                if len(win.T) > srate*(15):
                    win = eeg_record[:, speech + srate*(i) : speech + srate*(15+i)]
                    print("over")

                # Check print
                print("Window number : {0}".format(i+1))
                print("Time Check : {0}s".format(len(eeg_record[:, speech:].T) / srate))


                #----------------------------- Pre-processing -----------------------------#
                # preprocessing_ha.py

                win = Preproccessing(win, srate, 0.5, 8, 750)    # data, sampling rate, low-cut, high-cut, filter order
                data_l = len(win.T)

                #------------------------------- Train set -------------------------------#
                if tr < train:  #int train
                    state = "Train set"

                    ###   mTRF train function   ###
                    model, tlag, inter = mtrf_train(stim_R[tr:tr+1, 64*(i):64*(i)+data_l].T, win.T, fs, Dir, tmin, tmax, reg_lambda)
                    
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

            #------------------------------- Test set -------------------------------#
                else:
                    state = "Test set"

                    ## Calculate Predicted signal ##
                    pred, r_l, p, mse = mtrf_predict(stim_L[tr:tr+1, 64*(i):64*(i)+data_l].T, win.T, model, fs, Dir, tmin, tmax, inter)
                    pred, r_r, p, mse = mtrf_predict(stim_R[tr:tr+1, 64*(i):64*(i)+data_l].T, win.T, model, fs, Dir, tmin, tmax, inter)

                    # Stock correlation value per window(i)
                    r_L = np.append(r_L, r_l)
                    r_R = np.append(r_R, r_r)

                    ###### Estimate accuracy #####
                    if r_r > r_l:
                        acc = 1
                    else:
                        acc = 0

                    print("-------------------------------")
                    print("acc : {0}".format(acc))
                    print("-------------------------------")

                    # Save acc for entire Accuracy
                    Acc = np.append(Acc, acc)
1
                    # Up window number
                    i = i + 1
                    # Time count
                    end = time.perf_counter()

                ## End one window

                #end = time.process_time()
                work = end - start
                print("working time = {0}s".format(work))

        #------------------------ End 60s - one trial ------------------------#
        ##### Question #####
        try :
            if tr+1 == file.TrNum[j] :
                print("Question Time")
                correct = []
                answer = []

                correct, answer = Question(j, file)

                Correct.append(correct)
                Answer.append(answer)
                j = j+1
        except KeyError:                # 마지막 질문 후 에러나서
            pass

        # Stack eeg_record per trial & Save
        EEG.append(eeg_record.T)
        AUX.append(aux_record.T)
        # Reset
        eeg_record = np.zeros((16,1))
        aux_record = np.zeros((3,1))

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
                model_wm = model_wt/(i*(tr+1))
                inter_wm = inter_wt/(i*(tr+1))
                model = model_wm
                inter = inter_wm

        elif state == "Test set":
            # Stack correlation value collected during one trial
            entr_L.append(r_L)
            entr_R.append(r_R)

            r_L = []
            r_R = []
            plt.close()

            # Collect Accuracy per trial
            ACC = np.append(ACC,np.mean(Acc))
            print("\n==================================\n")
            print("Present Accuracy = {0}%".format(ACC[-1]*100))
            print("\n==================================\n")

        # Save per trial // eeg, trigger, accuracy ,behavior
        EEG_all = np.asarray(EEG)
        AUX_all = np.asarray(AUX)

        scipy.io.savemat(path + '/save_data/E.mat', {'EEG': EEG_all})
        scipy.io.savemat(path + '/save_data/A.mat', {'AUX': AUX_all})
        scipy.io.savemat(path + '/save_data/RAW.mat', {'RAW': raw_data})
        scipy.io.savemat(path + '/save_data/TRIGGER.mat', {'TRIGGER': tri_data})
        scipy.io.savemat(path + '/save_data/Accuracy.mat', {'Acc': ACC})
        correct_all = np.asarray(Correct)
        scipy.io.savemat(path + '/save_data/Behavior.mat', {'Behavior': correct_all})

        # For Next trial
        tr = tr+1
        w = 1

#----------------------------- 30 trial End -----------------------------#
final = visual.TextStim(screen, text="수고하셨습니다.", height=70, color=[1, 1, 1], wrapWidth=2000)

text.draw()
screen.flip()
time.sleep(3)

port.close()
screen.close()
board.stop_stream()
board.release_session()
print("The End")


#### save ####
# mat save
answer_all = np.asarray(Answer)
scipy.io.savemat(path + '/save_data/Answer.mat', {'Answer': answer_all})

# np save
np.save(path+'/save_data/EEG', EEG_all)
np.save(path+'/save_data/A', AUX_all)
np.save(path+'/save_data/All_Accuracy', ACC)
entr_L = np.asarray(entr_L)
entr_R = np.asarray(entr_R)
np.save(path+'/save_data/All_correlation_right', entr_R)
np.save(path+'/save_data/All_correlation_left', entr_L )

                            
