function quick_push()
    % QUICK_PUSH Prompts for a commit message, stages all changes, commits, and pushes.
    
    % 1. Prompt for the commit message
    commit_msg = input('Enter commit message (or press Enter to abort): ', 's');
    
    % Abort if the message is empty
    if isempty(strtrim(commit_msg))
        fprintf('Action aborted: No commit message provided.\n');
        return;
    end
    
    fprintf('\n--- Staging changes (git add .) ---\n');
    [status_add, out_add] = system('git add .');
    disp(out_add);
    
    fprintf('--- Committing changes ---\n');
    % Safely format the commit command with the user's message
    commit_cmd = sprintf('git commit -m "%s"', commit_msg);
    [status_commit, out_commit] = system(commit_cmd);
    disp(out_commit);
    
    % Only push if the commit was successful (status 0) or if there was nothing to commit
    if status_commit == 0
        fprintf('--- Pushing to remote (git push) ---\n');
        [status_push, out_push] = system('git push');
        disp(out_push);
        
        if status_push == 0
            fprintf('Success! Code pushed to the repository.\n');
        else
            fprintf('Warning: Push failed. Check your connection or Git credentials.\n');
        end
    else
        fprintf('Notice: Commit failed or working tree is already clean.\n');
    end
end