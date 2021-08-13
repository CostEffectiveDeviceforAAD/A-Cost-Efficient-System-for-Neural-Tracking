# -------------------------------------- Preprocessing --------------------------------------#
import librosa
from scipy import signal
from scipy.signal import butter, lfilter, resample, filtfilt
from scipy.stats import stats
#from OpenBCI_lsl import *
import numpy as np
from pylab import *

def butter_bandpass_filter(data, lowcut, highcut, fs, order):   # IIR
    nyq = 0.5 * fs
    low = lowcut/nyq
    high = highcut/nyq
    b, a = signal.butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y

def FIR_filter(data, lowcut, highcut, fs, order):   # FIR
    nyq = 0.5 * fs
    #low = lowcut/nyq
    #high = highcut/nyq
    b = signal.firwin(order, [lowcut, highcut], window='hamming', pass_zero='bandpass', nyq = nyq)   # filter design
    #y = signal.lfilter(b, 1.0, data)
    y = filtfilt(b, 1.0, data)
    return y

#4 Filter order = /Normalised width of transition band

#print("a")
def Preproccessing(win,srate, low, high, order):

    #### Re-reference ####
    for t in range(0, len(win.T)):  # win : channel by samples
        re = np.mean(win[:, t])
        win[:, t] = win[:, t] - re

    '''
    #### Re-reference ####
    for c in range(0,15):
        m = np.mean(win[c,:])
        win[c,:] = win[c,:] - m


    for t in range(0, len(win.T)):  # win : channel by samples
        re_1 = np.mean(win[:7, t])
        re_2 = np.mean(win[7:, t])
        win[:7, t] = win[:7, t] - re_1
        win[7:, t] = win[7:, t] - re_2
    '''

    #### Filtering ####

    try:
        # FIR
        win = FIR_filter(win, low, high, srate, order)  # order 조정필요.
        # IIR
        #win = butter_bandpass_filter(win, low, high, srate, order)

    except ValueError:  # To extend length of data by using zeropad, when data is too short.
        win = np.concatenate((win, np.zeros((len(win), 50))), axis = 1)
        print("***** zero padding *****")
        win = FIR_filter(win, low, high, srate, order)  # order 조정필요.


    #### Resampling ####

    win = librosa.resample(win, srate, 64)

    ## Z-scoring

    win = (win - win.mean()) / win.std()
    #win = stats.zscore(win)

    return win

#Plot frequency and phase response
def mfreqz(b,a=1):
    w,h = signal.freqz(b,a)
    h_dB = 20 * log10 (abs(h))
    subplot(211)
    plot(w/max(w),h_dB)
    ylim(-150, 5)
    ylabel('Magnitude (db)')
    xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    title(r'Frequency response')
    subplot(212)
    h_Phase = unwrap(arctan2(imag(h),real(h)))
    plot(w/max(w),h_Phase)
    ylabel('Phase (radians)')
    xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    title(r'Phase response')
    subplots_adjust(hspace=0.5)

