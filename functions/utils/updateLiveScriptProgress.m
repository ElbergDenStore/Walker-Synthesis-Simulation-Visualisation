function updateLiveScriptProgress(total_pts, reset_flag)
    persistent p last_percent reverseStr
    
    % Initialization / Reset
    if nargin > 1 && reset_flag
        p = 0;
        last_percent = -1; 
        reverseStr = '';
        return;
    end
    
    if isempty(p)
        p = 0;
        last_percent = -1;
        reverseStr = '';
    end
    
    p = p + 1;
    current_percent = floor((p / total_pts) * 100);
    
    % Only update the screen when the percentage actually changes (prevents terminal lag)
    if current_percent > last_percent || p == total_pts
        bar_length = 40; % How wide you want the progress bar to be
        num_equals = round((current_percent / 100) * bar_length);
        num_spaces = bar_length - num_equals;
        
        % Build the string: e.g., [========          ]
        bar_str = ['[', repmat('=', 1, num_equals), repmat(' ', 1, num_spaces), ']'];
        
        % Create the full message
        msg = sprintf('Processing: %s %d%%', bar_str, current_percent);
        
        % Print backspaces to clear the old line, then print the new line
        fprintf([reverseStr, msg]);
        
        % Save the number of backspaces needed for the next loop
        reverseStr = repmat(sprintf('\b'), 1, length(msg)-1);
        
        last_percent = current_percent;
    end
    
    % Cap it off cleanly when finished and drop to a new line
    if p >= total_pts
        % fprintf('\n');
        p = 0; 
        last_percent = -1;
        reverseStr = ''; 
    end
end