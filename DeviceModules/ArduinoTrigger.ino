/////////////////////////////////////
////     Arduino - WAV Trigger   ////
////    For presentation sound   ////
////       providing trigger     ////
/////////////////////////////////////

#include <Metro.h>
#include <AltSoftSerial.h>    // Arduino build environment requires this
#include <wavTrigger.h>
wavTrigger wTrig;             // Our WAV Trigger object
Metro gLedMetro(500);         // LED blink interval timer
Metro gSeqMetro(6000);        // Sequencer state machine interval timer

void setup() {
  Serial.begin(9600);
  pinMode(13, OUTPUT);  // WAV Trigger onset
  pinMode(12, OUTPUT);  // Real sound onset
  pinMode(11, OUTPUT);  // 
  pinMode(10, OUTPUT);
  
  // If the Arduino is powering the WAV Trigger, we should wait for the WAV
  //  Trigger to finish reset before trying to send commands.
  delay(1000);
  wTrig.start();  
  delay(10);

  // Send a stop-all command and reset the sample-rate offset, in case we have
  // reset while the WAV Trigger was already playing.
  wTrig.stopAllTracks();   
  wTrig.samplerateOffset(0);  
}

int input = 0;
int track = 2;      // First audio track number
int thres = 11;    // Input Analog Sound threshold - should you check this
int prac  = 32;
void loop() {
  digitalWrite(13, LOW);
  digitalWrite(12, LOW);
  digitalWrite(10, LOW);

  int sound_L = analogRead(A0);
  int sound_R = analogRead(A1);
   
  // Trigger for start trial
  int input = Serial.read();
  if (input==49){      // Receving int = 1 
    
    // start sound
    wTrig.trackPlaySolo(track);
    digitalWrite(13, HIGH);      // Toward Cyton board /Trigger for Checking WAV Trigger, not sound onset
          
    while(true) {    
      // Detect analog signal of sound(beep sound)
      int sound_L = analogRead(A0);
      int sound_R = analogRead(A1);

      if ( sound_L > thres || sound_R > thres) {    
        digitalWrite(12, HIGH);  // Toward Cyton board /Trigger for onset sound
        digitalWrite(10,HIGH);   // Check LED (ouside)
        delay(62000);            // Play Wav file during 63s 
        track++;                 // for next track
        break;          
      }          
    }
  }  
  else if (input == 51){    // practice
    wTrig.trackPlaySolo(prac);  // the track for practice
        digitalWrite(10,HIGH);
        delay(15000);     
        prac++;
      } 
}


              
