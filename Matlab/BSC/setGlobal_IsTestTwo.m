function setGlobal_IsTestTwo
    
    % turn-on training toggle
    global IsTestTwo
    IsTestTwo = 1;
    
    % initialize decoding weights and bias
    global rL; global rL_biased;
    rL = []; rL_biased = [];
    
    global rR; global rR_biased;
    rR = []; rR_biased = [];
    
    % initialize accs
    global acc; global acc_biased;
    acc = []; acc_biased = [];
    
end