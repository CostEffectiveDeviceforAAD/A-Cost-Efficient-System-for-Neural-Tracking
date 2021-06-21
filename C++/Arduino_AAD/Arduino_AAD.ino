// Arduino - WAV Trigger

#include <Metro.h>
#include <AltSoftSerial.h>    // Arduino build environment requires this
#include <wavTrigger.h>

wavTrigger wTrig;             // Our WAV Trigger object

Metro gLedMetro(500);         // LED blink interval timer
Metro gSeqMetro(6000);        // Sequencer state machine interval timer

void setup() {

  // Serial monitor
  Serial.begin(9600);
  pinMode(13, OUTPUT);
  pinMode(12, OUTPUT);
  pinMode(11, OUTPUT);
  
  
  // If the Arduino is powering the WAV Trigger, we should wait for the WAV
  //  Trigger to finish reset before trying to send commands.
  delay(1000);

  // WAV Trigger startup at 57600
  wTrig.start();  
  delay(10);

  // Send a stop-all command and reset the sample-rate offset, in case we have
  //  reset while the WAV Trigger was already playing.
  wTrig.stopAllTracks();   
  wTrig.samplerateOffset(0);  

  // Sound volume control
  //wTrig.masterGain(0);
  
}

int input = 0;
int track = 1;

void loop() {
 
  digitalWrite(13, LOW);
  digitalWrite(11, LOW);
  
  // Trigger for start trial
  int input = Serial.read();
  if (input==49){  //int = 1 

    // start sound
    wTrig.trackPlaySolo(track);
    // Trigger for just
    digitalWrite(11, HIGH);
          
    // Check LED (ouside)
    digitalWrite(13,HIGH);
    
    while(true) {    

      // Detect analog signal of sound(beep sound)
      int sound = analogRead(A0);
      
      if ( sound > 11 ) { // threshold =11

        /*
        ////////////// TEST //////////////////
        for(int i = 0; i < 32; i++){
          digitalWrite(11, LOW);
          digitalWrite(12, HIGH);   // trigger for start speech!
          delay(1000);
          digitalWrite(11, HIGH);
          digitalWrite(12, LOW);
          delay(1000); }
        ///////////////////////////////////////
        */
        digitalWrite(12,HIGH);
        delay(60000);
        digitalWrite(12,LOW);
        
        
        // for next track
        track++;
        // for next trial
        break;          
      }       
      
    }
  }
  
}

              
