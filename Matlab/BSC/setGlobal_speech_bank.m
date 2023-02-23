function setGlobal_speech_bank
    
    temp1 = load('AAK_left.mat');
    global left_speech_bank
    left_speech_bank = temp1.AAK_left;
    
    temp2 = load('AAK_right.mat');
    global right_speech_bank
    right_speech_bank = temp2.AAK_right;
    
end    