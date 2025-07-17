function generate_MC_paper_figures(name_var)

files = {...
        {['ekf_' name_var '.mat'] , ['MC_grafik_EKF_' name_var '.pdf']}, ...
        {['ukf_' name_var '.mat'] , ['MC_grafik_UKF_' name_var '.pdf']}, ...
        {['pf_'  name_var '.mat'] , ['MC_grafik_PF_'  name_var '.pdf']}};

for i = 1:length(files)
    inputFile = files{i}{1}; outputFile = files{i}{2};
    fprintf('\n--- Memproses file: %s ---\n', inputFile);
    plot_figure(inputFile, outputFile, 1);
end
end