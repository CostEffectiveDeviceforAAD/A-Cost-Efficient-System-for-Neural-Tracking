function setGlobal_dir(inevents, dir_marker)
    
    % load dir marker
    global dir
    
    % assign direction number
    if any(ismember(inevents(2,:), dir_marker(1))) % Left
        dir = 1;
    elseif any(ismember(inevents(2,:), dir_marker(2))) % Right
        dir = 2;
    end
    
end    