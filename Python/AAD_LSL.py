
""""""""""""""""""""""""""""""""""""""""""
 #        OpenBCI - Python LSL           #
 #           For Online AAD              #
""""""""""""""""""""""""""""""""""""""""""


#================================== SET EXPERIMENT ================================================#

###### Imports #####
import librosa, warnings, random, time, os, sys, serial
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pylsl import StreamInlet, resolve_stream, StreamInfo
#from OpenBCI_lsl import *
from scipy import signal
from scipy.signal import butter, lfilter, resample, filtfilt
#from helper import *
from pymtrf import *
from psychopy import visual, core, event
from preprocessing_ha import *


#----------------------------- Load Speech segment data ------------------------------#

stim_L = np.load('C:/Users/LeeJiWon/Desktop/OpenBCI/AAD/AAK/ORIGINAL_SPEECH/Stim_seg_L.npy')   ## 30*46*960 numpy array
stim_R = np.load('C:/Users/LeeJiWon/Desktop/OpenBCI/AAD/AAK/ORIGINAL_SPEECH/Stim_seg_R.npy')   ## 30*46*960 numpy array


# kist : 'C:/Users/LeeJiWon/Desktop/OpenBCI/AAD/AAK/ORIGINAL_SPEECH/'
# hyu : 'C:/Users/user/Desktop/hy-kist/OpenBCI/save_data/'

#----------------------------- Connect to port of arduino ------------------------------#

port = serial.Serial("COM8", 9600)

# kist - COM8
#----------------------------- Open LSL network -----------------------------#

# first resolve an EEG stream on the lab network
print("looking for an EEG stream...")
streams_eeg = resolve_stream('type', 'EEG')
streams_aux = resolve_stream('type', 'AUX')


# create a new inlet to read from the stream
print("StreamInlet")
inlet_eeg = StreamInlet(streams_eeg[0], 360, 125)    # channel
inlet_aux = StreamInlet(streams_aux[0])            # aux


#----------------------------- Parameter Setting -----------------------------#

srate = 125
fs = 64
tmin = 0
tmax = 250
Dir = -1
reg_lambda = 10
train = 14

# Set int
r_L = []
r_R = []
Acc = []
ACC = []
model_w = []
inter_w = []
entr_L =[]
entr_R = []
EEG = []
tr = 0
start = []
end = []
w = 1
s = 0

#----------------------------- Make the window for Psychopy -----------------------------#

screen = visual.Window([960, 900],
    screen = 0,
    pos = [600,0],
    fullscr = False,
    winType = 'pyglet',
    allowGUI = False,
    allowStencil = False,
    monitor ='testMonitor',
    color = [-1,-1,-1],
    blendMode = 'avg',
    units = 'pix'
    #pos = [100,0]
                    )

# Set Text_1 - Start
text = visual.TextStim(screen, text=" + ", height=100, color=[1, 1, 1])

# Draw text
text.draw()
screen.flip()

#==================================================================================================#
#-------------------------------------- START EXPERIMENT ------------------------------------------#
#==================================================================================================#

##### Start 30 trial #####
while tr < 30:   # 30

    #inlet_eeg = StreamInlet(streams_eeg[0], 360, 125)  # channel
    #inlet_aux = StreamInlet(streams_aux[0])

    # Trigger detection
    [sample_aux, ts_aux] = inlet_aux.pull_sample()
    print("{0}".format(sample_aux))

#----------------------------- Psychopy Window & Serial Write ------------------------------#

    # Press Button for start
    key = event.getKeys()
    if key == ["space"] and tr == 0:

        # Send signal to arduino for start sound
        port.write(b'1')

        # Set Text_2
        text2 = visual.TextStim(screen, text="<<<<", height=80, color=[1, 1, 1])
        # Draw text
        text2.draw()
        screen.flip()
        s = 1;

    elif tr > 0 and w == 1 :

        # Set Text_1 - Start
        text = visual.TextStim(screen, text=" + ", height=100, color=[1, 1, 1])

        # Draw text
        text.draw()
        screen.flip()

        [sample_aux, ts_aux] = inlet_aux.pull_sample()
        print("wait")
        time.sleep(3)

        # Send signal to arduino for start sound
        port.write(b'1')

        # Set Text_2
        text2 = visual.TextStim(screen, text="<<<<", height=80, color=[1, 1, 1])
        # Draw text
        text2.draw()
        screen.flip()
        w = 0

    [sample_aux, ts_aux] = inlet_aux.pull_sample()
    [sample_eeg, ts_eeg] = inlet_eeg.pull_chunk(0, 125)
    print("{0}".format(sample_aux))

#----------------------------- Trigger detection -----------------------------#

    if sample_aux[1] > 0 and s == 1:

        print("Input Trigger {0}".format(tr))

        # Wait beep sound
        time.sleep(3)

        print("Start Speech")
        [sample_eeg, ts_eeg] = inlet_eeg.pull_chunk(0, 125)
        # Format per trial
        eeg_record = np.array([])
        i = 0
        work = 1
        eeg = []


#----------------------------- Working while 60s -----------------------------#

        while len(eeg_record.T) != 7500:

            if work > 1:
                work = 1

            time.sleep(1 - work)

            start = time.perf_counter()

            # receive sample per 1s

            [sample_eeg, ts_eeg] = inlet_eeg.pull_chunk(0, 125)
            [sample_aux, ts_aux] = inlet_aux.pull_sample()

            # To prevent to precess with empty samples
            if sample_eeg:

                eeg = eeg + sample_eeg   # list - a number of samples

                if len(sample_eeg) < 125:
                    sample = sample.expend(sample, 0, (125-len(sample_eeg)))

            end = time.perf_counter()
            work = end - start


            # Stack samples until 15s and window sliding per 1s
            if len(eeg) >= srate * (i+15):

                eeg_record = np.asarray(eeg).T    # channel by sample

                # zero padding for lack of sample
                if len(sample_eeg) < 125:

                    eeg_record = np.pad(eeg_record,((0,0),(0,(i+15)*srate-len(eeg_record.T))), 'constant', constant_values = 0)  # 16 by 1875

                # Adjust the window length to match the seconds
                win = eeg_record[:, round(i * srate):(i+15)*srate]

                # Check print
                print("Window number : {0}".format(i))
                print("Time Check : {0}s".format(len(eeg_record.T) / srate))


#----------------------------- Pre-processing -----------------------------#
                # preprocessing_ha.py

                win = Preproccessing(win, srate, 0.5, 8, 3)  # data, sampling rate, low-cut, high-cut, filt order



#------------------------------- Train set -------------------------------#
                if tr < 14:  #int train
                    state = "Train set"

                    ## mTRF train function ##
                    model, tlag, inter = mtrf_train(stim_L[tr,i:i+1].T, win.T, fs, Dir, tmin, tmax, reg_lambda)
                    
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
                    pred, r_l, p, mse = mtrf_predict(stim_L[tr, i:i+1].T, win.T, model, fs, Dir, tmin, tmax, inter)
                    pred, r_r, p, mse = mtrf_predict(stim_R[tr, i:i+1].T, win.T, model, fs, Dir, tmin, tmax, inter)

                    # Stock correlation value per window(i)
                    r_L = np.append(r_L, r_l)
                    r_R = np.append(r_R, r_r)

                    ## Real-time Plotting ##
                    plt.clf()
                    if i == 0:
                        plt.ion()
                        fig, ax1 = plt.subplots()
                        ax2 = ax1.twiny()

                    x = np.arange(14, i + 15)

                    plt.plot(x, r_L, 'ob-', label='Left')
                    plt.plot(x, r_R, 'or-', label='Right')

                    # trial labeling
                    plt.ylabel("Correlation")
                    plt.xlabel("Time")
                    plt.grid(True)
                    plt.legend()
                    plt.xlim(0, 60)
                    plt.ylim(-0.3, 0.3)
                    fig.canvas.draw()
                    fig.canvas.flush_events()
                    plt.draw()

                    ###### Estimate accuracy #####
                    if r_l > r_r:

                        acc = 1

                    else:

                        acc = 0

                    #print("acc : {0}".format(acc))

                    # acc save for entire Accuracy
                    Acc = np.append(Acc, acc)

                    i = i + 1
                    end = time.perf_counter()
                # End one window

                #end = time.process_time()
                work = end - start
                print("working time = {0}s".format(work))

#----------------------------- End 60s - one trial -----------------------------#

        # Stack eeg_record per trial
        EEG.append(eeg_record)

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
                model_wm = model_wt/(i*tr)
                inter_wm = inter_wt/(i*tr)
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
            ACC = np.append(ACC,mean(Acc))
            print("\n==================================\n")
            print("Present Accuracy = {0}%".format(ACC[-1]*100))
            print("\n==================================\n")

        # Next trial
        tr = tr+1
        w = 1

#----------------------------- 30 trial End -----------------------------#
port.close()
screen.close()
print("The End")

#### save ####
EEG = np.asarray(EEG)
entr_L = np.asarray(entr_L)
entr_R = np.asarray(entr_R)
np.save('C:/Users/user/Desktop/hy-kist/OpenBCI/save_data/EEG_record', EEG)
np.save('C:/Users/user/Desktop/hy-kist/OpenBCI/save_data/All_Accuracy', ACC)
np.save('C:/Users/user/Desktop/hy-kist/OpenBCI/save_data/All_correlation_right', entr_R)
np.save('C:/Users/user/Desktop/hy-kist/OpenBCI/save_data/All_correlation_left', entr_L )

                            

                    



