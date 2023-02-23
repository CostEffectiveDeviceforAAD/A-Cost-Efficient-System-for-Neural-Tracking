function setGlobal_IsTestOne
    
    % turn-on test toggle
    global IsTestOne
    IsTestOne = 1;
    
    
    % load decoder models
    temp=load('C:\Users\Seung-Cheol Baek\AppData\Roaming\Neuroscan\Curry 8\Matlab\models.mat');
    
    % BOTH
    global model
    model = temp.models.model;
    
    % LEFT
    global model_l
    model_l = temp.models.model_l;
    
    % RIGHT
    global model_r
    model_r = temp.models.model_r;
    
    
    % initialize decoding weights and bias
    global rL; global rL_biased;
    rL = []; rL_biased = [];
    
    global rR; global rR_biased;
    rR = []; rR_biased = [];
    
    % initialize accs
    global acc; global acc_biased;
    acc = []; acc_biased = [];
    
end    