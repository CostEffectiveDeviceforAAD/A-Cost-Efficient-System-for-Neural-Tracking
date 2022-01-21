import pandas as pd
import librosa
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert, chirp
from scipy.stats import stats

# C:\Users\user\Desktop\hy-kist\aak\ORIGINAL_SPEECH\L_Twenty

####### Twenty
# make sample
speech = np.random.random((1,1323000))

for i in range(1,31):
    audio = 'hy-kist/aak/ORIGINAL_SPEECH/L_Twenty/L'+str(i)+'.wav'
    y, sr = librosa.load(audio)       #sr = 22050
    speech = np.vstack((speech, y))

speech = np.delete(speech,0,0)        # delete 30 by 132300

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

####### Journey
speech = np.random.random((1,1323000))

for i in range(1,31):
    audio = 'hy-kist/aak/ORIGINAL_SPEECH/R_Journey/R'+str(i)+'.wav'
    y, sr = librosa.load(audio)       #sr = 22050
    speech = np.vstack((speech, y))

speech = np.delete(speech,0,0)

for i in range(0,30):
    y = speech[i]
    y = hilbert(y)
    y = np.abs(y)
    y = librosa.resample(y, sr, 64)
    y = stats.zscore(y)

    Allspeech[i+30] = pd.DataFrame(y)

np.save('hy-kist/Allspeech',Allspeech.T)


