import sys, re
file_path = "c:\\Users\\StoreElberg\\OneDrive - Aalborg Universitet\\10_Semester\\matlab_code\\coverage_simulator_function.m"
with open(file_path, 'r') as f: text = f.read()

# 1. Signature Cfg.use_parallel
text = re.sub(r'(% Outputs: metrics - Struct containing key performance indicators for optimization\n)', r'\1    Cfg.use_parallel = use_parallel;\n', text)

# 2. DataQueue for plot removed
text = text.replace("plot_dq = parallel.pool.DataQueue;\n        updateLiveScriptProgress(num_plots, true); \n        afterEach(plot_dq, @(~) updateLiveScriptProgress(num_plots, false));", "updateLiveScriptProgress(num_plots, true);")
text = re.sub(r'^[ \t]*send\(plot_dq, \[\]\);', r'        updateLiveScriptProgress(num_plots, false);', text, flags=re.MULTILINE)

with open(file_path, 'w') as f: f.write(text)
