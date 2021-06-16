
""""""""""""""""""""""""""""""""""""""""""
 #        OpenBCI - Python LSL           #
 #           For Online AAD              #
""""""""""""""""""""""""""""""""""""""""""


#================================== SET EXPERIMENT ================================================#

###### Imports #####
import librosa, warnings, random, time, os, sys, serial, logging, argparse, mne
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
from brainflow.board_shim import BoardShim, BrainFlowInputParams, BoardIds, BrainFlowError
from brainflow.data_filter import DataFilter, FilterTypes, AggOperations, WindowFunctions, DetrendOperations


#----------------------------- Load Speech segment data ------------------------------#

stim_L = np.load('C:/Users/LeeJiWon/Desktop/OpenBCI/AAD/AAK/ORIGINAL_SPEECH/Stim_seg_L.npy')   ## 30*46*960 numpy array
stim_R = np.load('C:/Users/LeeJiWon/Desktop/OpenBCI/AAD/AAK/ORIGINAL_SPEECH/Stim_seg_R.npy')   ## 30*46*960 numpy array


# kist : 'C:/Users/LeeJiWon/Desktop/OpenBCI/AAD/AAK/ORIGINAL_SPEECH/'
# hyu : 'C:/Users/user/Desktop/hy-kist/OpenBCI/save_data/'

#----------------------------- Connect to port of arduino ------------------------------#

port = serial.Serial("COM8", 9600)

# kist - COM8
#----------------------------- Open Brainflow network -----------------------------#

BoardShim.enable_dev_board_logger()

parser = argparse.ArgumentParser()
# use docs to check which parameters are required for specific board, e.g. for Cyton - set serial port
parser.add_argument('--timeout', type=int, help='timeout for device discovery or connection', required=False,
                    default=0)
parser.add_argument('--ip-port', type=int, help='ip port', required=False, default=0)
parser.add_argument('--ip-protocol', type=int, help='ip protocol, check IpProtocolType enum', required=False,
                    default=0)
parser.add_argument('--ip-address', type=str, help='ip address', required=False, default='')
parser.add_argument('--serial-port', type=str, help='serial port', required=False, default='COM7')      # kist : COM7 / hy: COM15
parser.add_argument('--mac-address', type=str, help='mac address', required=False, default='')
parser.add_argument('--other-info', type=str, help='other info', required=False, default='')
parser.add_argument('--streamer-params', type=str, help='streamer params', required=False, default='')
parser.add_argument('--serial-number', type=str, help='serial number', required=False, default='')
parser.add_argument('--board-id', type=int, help='board id, check docs to get a list of supported boards',
                    required=False, default='2')
parser.add_argument('--file', type=str, help='file', required=False, default='')
args = parser.parse_args()

params = BrainFlowInputParams()
BoardShim.enable_dev_board_logger()
params.ip_port = args.ip_port
params.serial_port = args.serial_port
params.mac_address = args.mac_address
params.other_info = args.other_info
params.serial_number = args.serial_number
params.ip_address = args.ip_address
params.ip_protocol = args.ip_protocol
params.timeout = args.timeout
params.file = args.file

board = BoardShim(args.board_id, params)
board.prepare_session()

# board.start_stream () # use this for default options
board.start_stream(45000, args.streamer_params)
time.sleep(10)
#data = board.get_current_board_data (125) # get latest 256 packages or less, doesnt remove them from internal buffer
input = board.get_board_data()

# Set channels number
eeg_channels = board.get_eeg_channels(args.board_id)
aux_channels = board.get_analog_channels(args.board_id)

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
AUX = []
tr = 0
start = []
end = []
w = 1
s = 0
temp=[]

eeg_record = np.zeros((16,1))
aux_record = np.zeros((3,1))
eeg_data = np.array([])
temp=np.zeros((16,1))
w = 1
num = 1
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
while tr < 1:   # 30

    #inlet_eeg = StreamInlet(streams_eeg[0], 360, 125)  # channel
    #inlet_aux = StreamInlet(streams_aux[0])

#------------ Each trial

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

        #reset
        eeg_record = np.zeros((16,1))
        aux_record = np.zeros((3,1))

    # Trigger detection
    input = board.get_board_data()
    eeg_data = input[eeg_channels, :]
    aux_data = input[aux_channels, :]
    eeg_record = np.concatenate((eeg_record, eeg_data), axis=1)
    aux_record = np.concatenate((aux_record, aux_data), axis=1)
    print(aux_data)

#----------------------------- Trigger detection -----------------------------#
# per trial
    if 1 in aux_data[1,:]:

        print("Input Trigger {0}".format(tr))

        print("Start Speech")

        # Find onset point
        index = np.where(aux_record[1,:]==1)     # Onset 지점!
        onset = index[0][0]

        # Format per trial
        #eeg_record = np.array([])
        i = 0
        work = 1


#----------------------------- Working while 60s -----------------------------#
        # onset 부터 3초 지나고 원하는 시간(한 trial) 동안 돌아가도록
        speech = onset + (srate*3)
        while len(eeg_record.T) < (speech + srate*60):

            if work > 1:
                work = 1

            time.sleep(1 - work)

            start = time.perf_counter()

            # Receive sample
            input = board.get_board_data()
            eeg_data = input[eeg_channels, :]
            aux_data = input[aux_channels, :]                 # 11,12,13 / 0 or 1
            eeg_record = np.concatenate((eeg_record, eeg_data), axis=1)     # channel by time
            aux_record = np.concatenate((aux_record, aux_data), axis=1)


            end = time.perf_counter()
            work = end - start


            # Stack samples until 15s and window sliding per 1s
            if len(eeg_record.T) >= (speech + srate * (i+15)):

                # Adjust the window length to match the seconds
                win = eeg_record[:, speech + srate*(i) : speech + srate*(i+15)]
                trig = aux_record[:, speech + srate*(i) : speech + srate*(i+15)]


                # Check print
                print("Window number : {0}".format(i))
                print("Time Check : {0}s".format(len((eeg_record.T)-speech) / srate))


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
        AUX.append(aux_record)

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

                            

                    



