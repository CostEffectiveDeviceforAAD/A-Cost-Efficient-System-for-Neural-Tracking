function setGlobal_IsTrain
    
    % turn-on training toggle
    global IsTrain
    IsTrain = 1;
    
    % initialize decoding weights and bias
    global w_l % for left
    w_l = zeros(60,17);
    
    global w_r % for right
    w_r = zeros(60,17);
    
    global b_l % for left
    b_l = zeros(1);
    
    global b_r % for right
    b_r = zeros(1);
    
end