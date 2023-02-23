function offGlobal_IsTestOne
    
    % turn-on training toggle
    global IsTestOne
    IsTestOne = 0;
    
    
    % call correlation coefficients
    global rL; global rL_biased;    
    global rR; global rR_biased;

    rs_L = [rL; rL_biased];
    rs_R = [rR; rR_biased];
    
    
    % call accs
    global acc; global acc_biased;
    accs = [acc; acc_biased];
    
    
    % all results into structure
    test1_Result.rs_L = rs_L;
    test1_Result.rs_R = rs_R;
    test1_Result.accs = accs;
    
    
    % save result structure
    save('test1_Result.mat', 'test1_Result')
    
end    