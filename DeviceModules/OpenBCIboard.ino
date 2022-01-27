///////////////////////////////////////
/////                             /////
/////         OpenBCI board       /////
/////     For acquisition data    /////
/////                             /////
///////////////////////////////////////

#include <DSPI.h>
#include <OBCI32_SD.h>
#include <EEPROM.h>
#include <OpenBCI_32bit_Library.h>
#include <OpenBCI_32Bit_Library_Definitions.h>

boolean addAccelToSD = false; // On writeDataToSDcard() call adds Accel data to SD card write
boolean addAuxToSD = false; // On writeDataToSDCard() call adds Aux data to SD card write
boolean SDfileOpen = false;
int sound_tri;
int sound_onset;
int defau;

void setup() {
  // Bring up the OpenBCI Board
  board.begin();
  //startFromScratch();
  board.useAccel(false);
}

void loop() {
  if (board.streaming) {
    if (board.channelDataAvailable) {
      // Read from the ADS(s), store data, set channelDataAvailable flag to false
      board.updateChannelData();

      // Read trigger  
      sound_tri = digitalRead(13);    // Trigger about onset of WAV Trigger
      sound_onset = digitalRead(12);  // Detect sound onset trigger
      defau = digitalRead(11);

      // Write trigger in cyton board
      board.auxData[0] = sound_tri;     // D13
      board.auxData[1] = sound_onset;   // D12
      board.auxData[2] = defau;         // D11
      
   // Send packet with channel data and auxData contents
      board.sendChannelData();   
    }
  }
  // Check the serial port for new data
  if (board.hasDataSerial0()) {
    // Read one char from the serial port
    char newChar = board.getCharSerial0();

    // Send to the board library
    board.processChar(newChar);
  }
  board.loop();
}
