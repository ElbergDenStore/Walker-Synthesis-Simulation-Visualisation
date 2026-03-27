% profile_validator.m
disp('======================================================');
disp('Starting Profiler for Validator Script');
disp('======================================================');

profile on;

% Run the validation script
validate_simulators;

% Save profiler output
disp('Wrapping up Profile report...');
profile off;
out_dir = fullfile('simulation_output', 'validator_profile_results');
if exist(out_dir, 'dir')
    rmdir(out_dir, 's');
end
profsave(profile('info'), out_dir);

disp('======================================================');
fprintf('Profiling complete! Open `%s` in your browser.\n', fullfile(out_dir, 'index.html'));
disp('======================================================');
