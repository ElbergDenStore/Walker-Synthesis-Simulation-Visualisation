import sys, re
file_path = "c:\\Users\\StoreElberg\\OneDrive - Aalborg Universitet\\10_Semester\\matlab_code\\coverage_simulator_function.m"
with open(file_path, 'r') as f: text = f.read()

# 1. Signature Cfg.use_parallel
text = re.sub(r'(% Outputs: metrics - Struct containing key performance indicators for optimization\n)', r'\1    Cfg.use_parallel = use_parallel;\n', text)

# 2. link_calc_matrix -> link_calc
text = text.replace('UL_Result = link_calc_matrix(el_mat, range_mat, lat_vec, lon_vec, Cfg.UL, Cfg);', 'UL_Result = link_calc(el_mat, range_mat, lat_vec, lon_vec, Cfg.UL, Cfg);')

# 3. throughput_mean_per_UE initialization
text = text.replace("metrics.throughput_mean   = NaN;", "metrics.throughput_mean   = NaN;\n    metrics.throughput_mean_per_UE   = NaN(NumUEs, 1);")

# 4. throughput_mean_per_UE assignment
text = text.replace("metrics.throughput_mean  = mean(all_thpt(:), 'omitnan');", "metrics.throughput_mean  = mean(all_thpt(:), 'omitnan');\n            metrics.throughput_mean_per_UE = mean(all_thpt, 2, 'omitnan');")

# 5. DataQueue for plot removed
text = text.replace("plot_dq = parallel.pool.DataQueue;\n        updateLiveScriptProgress(num_plots, true); \n        afterEach(plot_dq, @(~) updateLiveScriptProgress(num_plots, false));", "updateLiveScriptProgress(num_plots, true);")
text = re.sub(r'^[ \t]*send\(plot_dq, \[\]\);', r'        updateLiveScriptProgress(num_plots, false);', text, flags=re.MULTILINE)

# 6. ecdf fix
text = text.replace("ecdf(all_thpt./thpt_scale)", "ecdf(all_thpt(:)./thpt_scale)")

# 7. lat_vector cellfun
text = text.replace("lat_vector = cellfun(@(x) x.Lat, UEs)';\n        lon_vector = cellfun(@(x) x.Lon, UEs)';", "lat_vector = [UEs.Lat]';\n        lon_vector = [UEs.Lon]';")

# 8. Map 8 meanThroughput undefined
text = text.replace("thpt_vals = meanThroughput(:) / thpt_scale;", "thpt_vals = metrics.throughput_mean_per_UE(:) / thpt_scale;")

with open(file_path, 'w') as f: f.write(text)
