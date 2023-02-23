function offGlobal_IsTrain
    
    % turn-on training toggle
    global IsTrain
    IsTrain = 0;
    
    % call decoder params
    global w_l; global b_l;
    global w_r; global b_r;
    
    
    % BOTH - decoder params into structure
    model.w    = (w_l + w_r)/(46*14); % divide through Nframes by NtrainTrials
    model.b    = (b_l + b_r)/(46*14); % divide through Nframes by NtrainTrials
    model.t    = [-250.000, -234.375, -218.750, -203.125, -187.500, ... % tmin = 0, tmax = 250 (17 steps)
                  -171.875, -156.250, -140.625, -125.000, -109.375, ... % dir = -1 (decoder model)
                  -93.750,  -78.125,  -62.500,  -46.875,  -31.250, -15.625, 0];
    model.fs   = 64;
    model.dir  = -1;
    model.type = 'multi';
    
    % LEFT - decoder params into structure
    model_l.w    = w_l/(46*7); % divide through Nframes by NtrainTrials
    model_l.b    = b_l/(46*7); % divide through Nframes by NtrainTrials
    model_l.t    = [-250.000, -234.375, -218.750, -203.125, -187.500, ... % tmin = 0, tmax = 250 (17 steps)
                  -171.875, -156.250, -140.625, -125.000, -109.375, ... % dir = -1 (decoder model)
                  -93.750,  -78.125,  -62.500,  -46.875,  -31.250, -15.625, 0];
    model_l.fs   = 64;
    model_l.dir  = -1;
    model_l.type = 'multi';
    
    % RIGHT - decoder params into structure
    model_r.w    = w_r/(46*7); % divide through Nframes by NtrainTrials
    model_r.b    = b_r/(46*7); % divide through Nframes by NtrainTrials
    model_r.t    = [-250.000, -234.375, -218.750, -203.125, -187.500, ... % tmin = 0, tmax = 250 (17 steps)
                  -171.875, -156.250, -140.625, -125.000, -109.375, ... % dir = -1 (decoder model)
                  -93.750,  -78.125,  -62.500,  -46.875,  -31.250, -15.625, 0];
    model_r.fs   = 64;
    model_r.dir  = -1;
    model_r.type = 'multi';
    
    
    % Concat models into single structure
    models.model   = model;
    models.model_l = model_l;
    models.model_r = model_r;
    
    % save decoder models
    save('C:\Users\Seung-Cheol Baek\AppData\Roaming\Neuroscan\Curry 8\Matlab\models.mat', 'models')
    
end