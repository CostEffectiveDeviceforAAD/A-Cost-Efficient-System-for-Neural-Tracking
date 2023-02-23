function setGlobal_onTrial(inevents, insampleratehz)

    global onTrial
    onTrial = -4; % for correcting start point

    % find a time point of the first audio marker within a trial 
    idx = find(inevents(2,:) == 1000001);
    relative_pnt = inevents(1,idx(1))/insampleratehz;
    
    global f_onTrial
    f_onTrial = 1 - relative_pnt; % flag for clipping data
    
    global trial_count % update trial_count
    trial_count = trial_count + 1;
    
end