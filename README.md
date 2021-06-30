# Auditory BCI for Online AAD

This is Auditory BCI for Real-time Auditory Attention Detection(AAD).

It can not only collect EEG data from human brain, but it can play sounds and record triggers corresponding to sound.


## Requirements

1. OpenBCI (for EEG data acquisition)

    : OpenBCI cyton + daisy

    : EEG Electrode Cap

    https://openbci.com/

2. Arduino UNO (for Trigger)

    https://www.arduino.cc/

3. WAV Trigger (for sound play)

    : SD card

    https://www.sparkfun.com/products/13660

---------------------------

## OpenBCI 

For OpenBCI board running with Arduino IDE, see the OpenBCI Tutorial and Library.

https://docs.openbci.com/docs/02Cyton/CytonProgram

https://github.com/OpenBCI/OpenBCI_Cyton_Library

------------------------------

## Python LSL

- Lab Streaming Layer (LSL) with python

https://docs.openbci.com/docs/06Software/02-CompatibleThirdPartySoftware/LSL

+ We use advanced LSL code ( Not OpenBCI_LSL library )

https://docs.openbci.com/docs/09Deprecated/Python

-------------------------------

## WAV Trigger

WAV Trigger Tutorial & Library

http://robertsonics.com/2015/04/25/arduino-serial-control-tutorial/

https://github.com/robertsonics/WAV-Trigger-Arduino-Serial-Library

--------------------
## Device Pradigm

1. Start streaming with Python. (LSL) > EEG acquisition start 
2. Start Experiment with Psychopy.
3. Sends the signal of sound onset to Arduino via the Serial Port. [ 1. line ]
4. Play sound from WAV Trigger via Arduino ( generated jitter ) and send Triggers on sound onset. [ 2. line ]
5. Sound flows through earphone to subject and analog signals of sound are sent to Arduino at the same time. [ 3. line ]
6. Detect sound onset timing in Arduino by analog signals of sound. [ 4. line ]
7. Sends the tigger for the sound onset to Cyton + daisy ( No jitter ). [ 5. line ]
8. Sends EEG data and Trigger to the desktop, received and recorded them by python. [ 6. line ]



![제목 없는 프레젠테이션](https://user-images.githubusercontent.com/85104167/120435956-4b9eda80-c3b9-11eb-9f6c-0167660a8a2e.jpg)
