
""""""""""""""""""""""""""""""""""""""""""
 #        OpenBCI - Python LSL           #
 #           For Online AAD              #
""""""""""""""""""""""""""""""""""""""""""


#================================== SET EXPERIMENT ================================================#

###### Imports #####
import librosa, warnings, random, time, os, sys, serial, logging, argparse, mne, scipy.io
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
from brainflow.board_shim import BoardShim, BrainFlowInputParams, BoardIds, BrainFlowError, LogLevels
from brainflow.data_filter import DataFilter, FilterTypes, AggOperations, WindowFunctions, DetrendOperations


#----------------------------- Load Speech segment data ------------------------------#

speech = np.load('C:/Users/LeeJiWon/Desktop/OpenBCI/AAD/AAK/Allspeech.npy')     # 60 by 3840  /1-30 : left by time / 31-60 righy by time // srat : 64

stim_L = speech[:30, :]       # 30 by 3840   // trial by time
stim_R = speech[30:, :]       # 30 by 3840   // trial by time


# kist : 'C:/Users/LeeJiWon/Desktop/OpenBCI/AAD/AAK/Allspeech.npy'
# hyu : 'C:/Users/user/Desktop/hy-kist/OpenBCI/Sound data/AAK/ORIGINAL_SPEECH'

######################

attend = "L"

######################

#----------------------------- Connect to port of arduino ------------------------------#

port = serial.Serial("COM8", 9600)

# kist = COM8
# hyu = COM10
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
parser.add_argument('--serial-port', type=str, help='serial port', required=False, default='COM9')      # kist : COM7 / hy: COM15
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
BoardShim.log_message(LogLevels.LEVEL_INFO.value, 'start sleeping in the main thread')
time.sleep(10)
#data = board.get_current_board_data (125) # get latest 256 packages or less, doesnt remove them from internal buffer
input = board.get_board_data()      # get all data and remove it from internal buffer

# Set channels number
eeg_channels = board.get_eeg_channels(args.board_id)
aux_channels = board.get_analog_channels(args.board_id)

srate = board.get_sampling_rate(args.board_id)
#print(srate)


#------------------------------------ Question ------------------------------------------------#

q = "C:/Users/LeeJiWon/Desktop/OpenBCI/AAD/Python/question.xlsx"
file = pd.read_excel(q)


def Question(j, file):

    # Question 1
    text3 = visual.TextStim(screen, text = file.tweenty_Q1[j], height=50, color=[1, 1, 1], wrapWidth=2000)
    text3.draw()
    screen.flip()

    key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

    answer.append(key)
    if file.tweenty_A1[j] == int(key[0]):
        correct.append("True")
    else:
        correct.append("False")

    # Question 2
    text3 = visual.TextStim(screen, text = file.tweenty_Q2[j], height=50, color=[1, 1, 1], wrapWidth=2000)
    text3.draw()
    screen.flip()

    key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

    answer.append(key)
    if file.tweenty_A1[j] == int(key[0]):
        correct.append("True")
    else:
        correct.append("False")

    # Question 3
    text3 = visual.TextStim(screen, text = file.journey_Q1[j], height=50, color=[1, 1, 1], wrapWidth=2000)
    text3.draw()
    screen.flip()

    key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

    answer.append(key)
    if file.tweenty_A1[j] == int(key[0]):
        correct.append("True")
    else:
        correct.append("False")

    # Question 4
    text3 = visual.TextStim(screen, text = file.journey_Q2[j], height=50, color=[1, 1, 1], wrapWidth=2000)
    text3.draw()
    screen.flip()

    key = event.waitKeys(keyList=['1', '2', '3', '4'], clearEvents=True)

    answer.append(key)
    if file.tweenty_A1[j] == int(key[0]):
        correct.append("True")
    else:
        correct.append("False")



#----------------------------- Parameter Setting -----------------------------#

#srate = 125
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
leng= []
answer = []
correct = []


check = np.zeros(((16,1)))
eeg_record = np.zeros((16,1))
aux_record = np.zeros((3,1))
eeg_data = np.array([])
temp=np.zeros((16,1))
w = 1
num = 1
j = 0
#----------------------------- Make the window for Psychopy -----------------------------#

screen = visual.Window([960, 900],
    screen = 0,
    pos = [600,0],
    fullscr = True,
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
while tr < 5:   # 30

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

        print("wait")


        # Send signal to arduino for start sound
        time.sleep(3)
        port.write(b'1')
        input = board.get_board_data()

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
    if 1 in aux_data[1,:]:      # 12 pin 에 trigger 가 들어오면 시작.

        print("Input Trigger {0}".format(tr))

        print("Start Speech")

        # Find onset point
        index = np.where(aux_record[1,:] != 0)     # Onset 지점찾기
        onset = index[0][0]

        # Format per trial
        #eeg_record = np.array([])
        i = 0
        work = 0


#----------------------------- Working while 60s -----------------------------#

        # onset 부터 3초 지나고 원하는 시간(한 trial) 동안 돌아가도록
        speech = onset + (srate*3) + 1   # 3초 후 부터가 speech 이기에 +1
        while i != 46:  # 46 번째 window 까지

            # process 가 1 초 넘었다면 time.sleep 없이 pass
            if work > 1:
                work = 1

            # Wait 1s
            time.sleep(1 - work)

            # time count
            start = time.perf_counter()


            ### Receive sample ###
            input = board.get_board_data()

            # Separate EEG, AUX
            eeg_data = input[eeg_channels, :]
            aux_data = input[aux_channels, :]                 # 11,12,13 / 0 or 1
            print(len(eeg_data.T))

            # Stack data
            eeg_record = np.concatenate((eeg_record, eeg_data), axis=1)     # channel by time
            aux_record = np.concatenate((aux_record, aux_data), axis=1)

            # Count time
            end = time.perf_counter()
            work = end - start


            # Stack samples until 15s and window sliding per 1s
            # 15s 부터는 time.sleep(1s) 기준으로
            if len(eeg_record.T) >= (speech + (srate * 15)):


                # Adjust the window length to match the seconds
                win = eeg_record[:, speech + srate*(i) :]
                trg = aux_record[:, speech + srate*(i) :]



                if len(win.T) > srate*(15):
                    win = eeg_record[:, speech + srate*(i) : speech + srate*(15+i)]
                    print("over")

                # Check print
                print("Window number : {0}".format(i))
                print("Time Check : {0}s".format(len(eeg_record[:, speech:].T) / srate))


#----------------------------- Pre-processing -----------------------------#
                # preprocessing_ha.py

                win = Preproccessing(win, srate, 0.5, 8, 3)  # data, sampling rate, low-cut, high-cut, filt order
                data_l = len(win.T)

#------------------------------- Train set -------------------------------#
                if tr < 14:  #int train
                    state = "Train set"


                    ## mTRF train function ##
                    if attend == "L":
                        model, tlag, inter = mtrf_train(stim_L[tr:tr+1, 64*(i):64*(i)+data_l].T, win.T, fs, Dir, tmin, tmax, reg_lambda)
                    elif attend == "R":
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


        # Question
        if tr+1 == file.TrNum[j] :

            Question(j,file)
            j = j+1



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
            ACC = np.append(ACC,np.mean(Acc))
            print("\n==================================\n")
            print("Present Accuracy = {0}%".format(ACC[-1]*100))
            print("\n==================================\n")

        # Next trial
        tr = tr+1
        w = 1

#----------------------------- 30 trial End -----------------------------#
port.close()
screen.close()
board.stop_stream()
board.release_session()
print("The End")

path = 'C:/Users/LeeJiWon/Desktop/OpenBCI/save_data/'
scipy.io.savemat(path+'E.mat', {'EEG':EEG})
scipy.io.savemat(path+'A.mat', {'AUX':AUX})

# kist : C:/Users/LeeJiWon/Desktop/OpenBCI/save_data/
# hyu : 'C:/Users/user/Desktop/hy-kist/OpenBCI/Test/'

#### save ####
EEG = np.asarray(EEG)
entr_L = np.asarray(entr_L)
entr_R = np.asarray(entr_R)
np.save(path+'EEG_record', EEG)
np.save(path+'All_Accuracy', ACC)
np.save(path+'All_correlation_right', entr_R)
np.save(path+'All_correlation_left', entr_L )

                            
