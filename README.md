# A Cost-efficient System for Neural Tracking

The cost-efficicent system is devised for measure neural tracking in real time. It comprises an EEG acquisition module, a sound player module and a sound trigger module (refer to the below diagram).

This instruction provides details on implementing the devised system, and includes the information for performing the real-time Auditory Attention Decoding (AAD) task which is typical auditory task for measuring the neural tracking.


## Essential Modules
<img width="60%" src="https://user-images.githubusercontent.com/85104167/180901773-ea07b295-08eb-4860-a95a-27ae7a6ec707.jpg"/>

### EEG Acquisition Module

 [OpenBCI board](https://openbci.com/?utm_source=google&utm_medium=cpc&utm_campaign=716348300&utm_content=openbci&gclid=Cj0KCQiA-eeMBhCpARIsAAZfxZBwfN8ei8seomxZ255WDN04UvwYix6hzXr-pJoc7drJViXE77-MirIaAnfWEALw_wcB) (Biosensingboard Cyton with Daisy) is used as EEG Acquisition Module. The EEG acquisition module acquire sound onset trigger and EEG signals by 16 channels. Acquired data are sent to PC via bluetooth USB dongle and streamed through the [Brainflow](https://github.com/brainflow-dev/brainflow) which in open-source data acquisition software. In order to carry out the functions mentioned above, you need to proceed with uploading follow code to the EEG Acquisition Module:  `AcquisitionOpenBCI.ino`

   The OpenBCI Tutorial and Library can be found on the website provided by OpenBCI Inc., as indicated below.
  > 
  > https://docs.openbci.com/docs/02Cyton/CytonProgram
  > 
  > https://github.com/OpenBCI/OpenBCI_Cyton_Library

  

### Sound Player Module
For this system, [WAV Trigger](https://github.com/robertsonics/WAV-Trigger-Arduino-Serial-Library) is worked as Sound Player Module. 
When the sound player module received a command signal from the sound trigger module, it played .wav files saved on a micro-SD card.

See Arduino Serial Contol Tutorial to operate WAV Trigger
> 
> http://robertsonics.com/2015/04/25/arduino-serial-control-tutorial/


### Sound Trigger Module
For the sound trigger module, an Arduino UNO board was used. 
Firstly, it receives the trial start command from a laptop or PC through the COM port. Secondly, it sends the command for sound playback to the sound player module. Additionally, to synchronize the timing information of sound onset with EEG data, the sound trigger module generates a sound onset trigger that is sent to the EEG acquisition module. 
To operate it as described above, please upload the following code to the Sound Trigger Module: `ArduinoTrigger.ino` 

## Experimental processing

The real-time AAD task is performed using `OnlineAAD_EXP.py`. This script includes data streaming, prepocessing, decoding process with mTRF and communication with arduino board.

By utilizing `OnlineAAD_EXP.py` in conjunction with the custom codes provided in the Sub_Functions file, you can execute the real-time AAD task.
Also, you need to check your arduino, bluetooth COM port number and your directory of files.

+ Other requires
> https://github.com/SRSteinkamp/pymtrf

> https://github.com/brainflow-dev/brainflow

***
## Enclosure
Encloeure is producted by 3D printer with custom-made design. The dimension of the Encloeure is within 24 cm x 15 cm x 6 cm including cover (width x length x height). A Cover is made of acrilic material. 
You can use following 3D design file:

`DeviceEnclosure.step`
 or 
`DeviceEnclosure.stl`


## Detail
The outside of the system include DC power adapter, a switch of power for the EEG acquisition module, two LED and COM port adapter (described from left side, shown as first figure). Two LED are for checking the power on the EEG acquisition module and the trigger state on the sound trigger module. Second figure shown a parallel port adapter for EEG electrode cap, a volume controller using potentialmeter and an earphones stereo adapter. All modules are connected by jumper-wires (shown as final figure). The EEG acquisition module and sound trigger module are supplied from 5V DC power and COM port, and the sound player module is powered from the sound trigger module (5V).

***

<img src="https://user-images.githubusercontent.com/85104167/142797442-7c8c5677-199c-4192-8cdf-e37cbf4d5fd9.jpg" width="600" height="400">
<img src="https://user-images.githubusercontent.com/85104167/142797446-1ed05680-9816-4fd7-a80c-fed93afa0ad8.jpg" width="600" height="400">
<img src="https://user-images.githubusercontent.com/85104167/142797452-4d86a22f-e608-44a9-a706-3fac1b7e39b9.jpg" width="600" height="400">

