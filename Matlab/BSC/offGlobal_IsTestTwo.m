function offGlobal_IsTestTwo
    
    % turn-on training toggle
    global IsTestTwo
    IsTestTwo = 0;
    
    % call correlation coefficients
    global rL; global rL_biased; 
    global rR; global rR_biased;

    rs_L = [rL; rL_biased];
    rs_R = [rR; rR_biased];
    
     % call accs
    global acc; global acc_biased;
    accs = [acc; acc_biased];
    
    
    % all results into structure
    test2_Result.rs_L = rs_L;
    test2_Result.rs_R = rs_R;
    test2_Result.accs = accs;
    
    
    % save result structure
    save('test2_Result.mat', 'test2_Result')
  
end    