# A Cost-effective device for Online AAD

We will propose a cost-effective device for online AAD library. It performs three functions required for online AAD experiment: EEG acquisition, Stimuli Presentation and Stimuli Trigger. To test this device, we conducted online AAD experiment using this device and custom-made code.

This instruction contains information on how to conduct online AAD experiment through the proposed device. 

refer to the detail process
> paper ... now preparing ...


***
## Essential Modules

### 1.  EEG Acquisition Module

This device used [OpenBCI](https://openbci.com/?utm_source=google&utm_medium=cpc&utm_campaign=716348300&utm_content=openbci&gclid=Cj0KCQiA-eeMBhCpARIsAAZfxZBwfN8ei8seomxZ255WDN04UvwYix6hzXr-pJoc7drJViXE77-MirIaAnfWEALw_wcB) (Cyton with Daisy) as EEG Acquisition Module. This can acquire 16 channels and send data to PC via bluetooth dongle. Networking system for real time streaming selected the [Brainflow](https://github.com/brainflow-dev/brainflow) because of good accessibility and compatibility with OpenBCI. 

  > For OpenBCI board running with Arduino IDE, see the OpenBCI Tutorial and Library.
  > 
  > https://docs.openbci.com/docs/02Cyton/CytonProgram
  > 
  > https://github.com/OpenBCI/OpenBCI_Cyton_Library


**Board Programming**

   The board receive EEG data and trigger data through EEG channels and jumper wire, and send these data to PC via bluetooth dongle.
   Refer to follow code:
    
`
AcquisitionOpenBCI.ino
`<br/>


### 2.  Stimuli Presentation
Stimuli are presented via [WAV Trigger](https://github.com/robertsonics/WAV-Trigger-Arduino-Serial-Library). Sound stimuli for experiment are plugged in SD card. This used Arduino serial port for programming and is supplied power from arduino board.

> See Arduino Serial Contol Tutorial for WAV Trigger
> 
> http://robertsonics.com/2015/04/25/arduino-serial-control-tutorial/



### 3.  Stimuli Trigger
We used Arduino UNO to synchronizate between trigger and EEG data. Arduino UNO communicate with laptop(or PC) for trial onset and with WAV Trigger for stimuli (i.e, speech) onset. Finally, this send the trigger for stimuli onset to OpenBCI board. All signals are deliveried to the IO pin through jumper wire.  

See custom-made code for Aduino UNO used in experiment.

`
ArduinoTrigger.ino
`

***
## Experimental processing

All experimental processing includes data streaming, prepocessing, decoding process with mTRF, visual presentation with Psychopy and communication with arduino. This experimental processing is operated throught custom-made code based on python. See the code for experimental processing.

If you want to start from the right of attention direction, 

`
OnlineAAD_EXP_R.py
`

If you want to start from the left of attention direction, 

`
OnlineAAD_EXP_L.py
`

To use the code above, you have to include all functional code in Python folder. Also, you need to check your arduino and bluetooth COM port number.

***
## Enclosure for the proposed device
Encloeure of proposed device is producted by 3D printer with custom-made design. The dimension of the Encloeure is within 24 cm x 15 cm x 6 cm including cover (width x length x height). A Cover is made of acrilic material separately.

'
Enclousere
'


## Detail
The proposed device consisted of


***

<img src="https://user-images.githubusercontent.com/85104167/142797442-7c8c5677-199c-4192-8cdf-e37cbf4d5fd9.jpg" width="600" height="400">
<img src="https://user-images.githubusercontent.com/85104167/142797446-1ed05680-9816-4fd7-a80c-fed93afa0ad8.jpg" width="600" height="400">
<img src="https://user-images.githubusercontent.com/85104167/142797452-4d86a22f-e608-44a9-a706-3fac1b7e39b9.jpg" width="600" height="400">

