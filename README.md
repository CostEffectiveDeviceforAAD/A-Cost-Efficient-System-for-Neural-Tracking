# Auditory BCI for real-time AAD

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
## Device diagram

1. Start streaming using python (LSL) > EEG acquisition start 
2. Start Experiment using psychopy
3. Send signal for sound onset to Arduino through Serial Port [ 1. line ]
4. play sound from WAV Trigger through arduino ( generate jitter ) [ 2. line ]
5. Flow Sound through earphon to subject and analog signal of sound to arduino [ 3. line ]
6. Detect sound onset timing in arduino by analog signal of sound [ 4. line ]
7. Send tigger about sound onset to Cyton + daisy [ 5. line ]
8. Send EEG data and Trigger to desktop, be received and recorded by python. [ 6. line ]
9. ~ 


![제목 없는 프레젠테이션](https://user-images.githubusercontent.com/85104167/120435956-4b9eda80-c3b9-11eb-9f6c-0167660a8a2e.jpg)
