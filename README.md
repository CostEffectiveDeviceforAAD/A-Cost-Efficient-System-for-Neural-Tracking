# A Cost-effective device for Online AAD

Proposed cost-effective device for online AAD Library

It performs three functions required for online AAD experiment: EEG acquisition, Stimuli Presentation and Stimuli Trigger.

This instruction contains information on how to conduct online AAD through the proposed device.

***
## EEG Acquisition Module
This device used [OpenBCI](https://openbci.com/?utm_source=google&utm_medium=cpc&utm_campaign=716348300&utm_content=openbci&gclid=Cj0KCQiA-eeMBhCpARIsAAZfxZBwfN8ei8seomxZ255WDN04UvwYix6hzXr-pJoc7drJViXE77-MirIaAnfWEALw_wcB) board (Cyton with Daisy) as EEG Acquisition Module. This can acquire 16 channels and send data to PC via bluetooth dongle.

Networking system for real time streaming selected the [Brainflow](https://github.com/brainflow-dev/brainflow) because of good accessibility and compatibility with OpenBCI. 

### Requirements
1. **OpenBCI Programming**

    For OpenBCI board running with Arduino IDE, see the OpenBCI Tutorial and Library.

>https://docs.openbci.com/docs/02Cyton/CytonProgram

>https://github.com/OpenBCI/OpenBCI_Cyton_Library


2. **OpenBCI Streaming**

    The board receive EEG data and trigger data through EEG channels and jumper wire, and send these data to PC via bluetooth dongle.
```
DataStreaming_OpenBCI.ino
```

***
## Stimuli Presentation
Stimuli are presented via [WAV Trigger](https://github.com/robertsonics/WAV-Trigger-Arduino-Serial-Library). 
Sound Stumuli for experiment are plugged in SD card. This used Arduino serial port for programming and is supplied power from arduino board.

See Arduino Serial Contol Tutorial for WAV Trigger
>http://robertsonics.com/2015/04/25/arduino-serial-control-tutorial/


***
## Stimuli Trigger





***
## Streaming and Visual presentatino via Python




## Enclosure for the proposed device
