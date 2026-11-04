%% 报告补充试验：负荷、功率容量、规划电压裕度与成本、个人B/C
ess_setup(fileparts(which('ess_setup')));
report_out=fullfile(pwd,'supplement');addpath(report_out,'-begin');set(groot,'defaultFigureVisible','off');
for trial={'supplement_basics','planning_default','planning_v095','planning_v0965','planning_energy2','personal_B','personal_C'}
 run(fullfile(report_out,[trial{1} '.m']));
end
disp('REPORT_TRIALS_COMPLETED');
