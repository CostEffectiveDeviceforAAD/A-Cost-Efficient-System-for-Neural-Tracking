import pandas as pd
import librosa
#import librosa.display
#import IPython.display
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as mpl
#import matplotlib.font_manager as fm

# C:\Users\user\Desktop\hy-kist\aak\ORIGINAL_SPEECH\L_Twenty

####### Left

speech = np.random.random((1,1323000)) # 틀만들기

for i in range(1,31):
    audio = 'hy-kist/aak/ORIGINAL_SPEECH/L_Twenty/L'+str(i)+'.wav'
    y, sr = librosa.load(audio)       #sr = 22050
    speech = np.vstack((speech, y))

speech = np.delete(speech,0,0)        #30 by 132300  틀 없앰
    

from scipy.signal import hilbert, chirp
from scipy.stats import stats


for i in range(0,30):

    y = speech[i]
    y = hilbert(y)
    y = np.abs(y)
    y = librosa.resample(y, sr, 64)
    y = stats.zscore(y)
    if i == 0:
        Allspeech = pd.DataFrame(y)   # 3840 by 30
    else:
        Allspeech[i] = pd.DataFrame(y)   


####### Right


speech = np.random.random((1,1323000)) # 틀만들기

for i in range(1,31):
    audio = 'hy-kist/aak/ORIGINAL_SPEECH/R_Journey/R'+str(i)+'.wav'
    y, sr = librosa.load(audio)       #sr = 22050
    speech = np.vstack((speech, y))

speech = np.delete(speech,0,0)        #30 by 132300  틀 없앰
    

from scipy.signal import hilbert, chirp
from scipy.stats import stats


for i in range(0,30):

    y = speech[i]
    y = hilbert(y)
    y = np.abs(y)
    y = librosa.resample(y, sr, 64)
    y = stats.zscore(y)

    Allspeech[i+30] = pd.DataFrame(y)

Allspeech = Allspeech.T
np.save('hy-kist/Allspeech',Allspeech)


'''''''''''''''''''''''''''''''''
### save speech segment ###
'''''''''''''''''''''''''''''''''

speech = np.load('C:/Users/user/Desktop/hy-kist/OpenBCI/Python_code/Allspeech.npy')  # openbci2 파일과 동일한 파일에 위치
speech_L = speech[0:30]
speech_R = speech[30:61]

### Speech window segment

# for tr in range(0,30)   tr포함 2차원으로 쌓으려면 이렇게 for 이중 그럼 30by46 쌓일듯
tr = 0
fs = 64
win_off = []
win_on = []
stim_L = []
stim_R = []
Stim_L = []
Stim_R = []

for tr in range(0,30):

    for i in range(0, 46):  # window j

        win_on.append((i) * 1)
        win_off.append((i) * 1 + 15)

        if i == 0:
            stim_L = speech_L[tr, win_on[i] * fs:win_off[i] * fs]  # 960by1 46개 > 960by46
            stim_R = speech_R[tr, win_on[i] * fs:win_off[i] * fs]

        else:
            stim_L = np.vstack([stim_L, speech_L[tr, win_on[i] * fs:win_off[i] * fs]])  # 960by1 46개 > 960by46
            stim_R = np.vstack([stim_R, speech_R[tr, win_on[i] * fs:win_off[i] * fs]])

    Stim_L.append(stim_L)   # list 30 ( trial)
    Stim_R.append(stim_R)   # list 30 ( trial)

    Stim_seg_L = np.asarray(Stim_L)
    Stim_seg_R = np.asarray(Stim_R)


np.save('C:/Users/user/Desktop/hy-kist/OpenBCI/save_data/Stim_seg_L', Stim_seg_L)
np.save('C:/Users/user/Desktop/hy-kist/OpenBCI/save_data/Stim_seg_R', Stim_seg_R)

