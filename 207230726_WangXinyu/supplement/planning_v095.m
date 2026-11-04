%%
lesson = ess_live_setup();
export_figures = false;
%%
Stu_ID = 207230726; %[control:slider:stu_id_0]{"position":[10,19]}
XY=mod(Stu_ID,100);
X_digit=floor(XY/10); Y_digit=mod(XY,10);
alpha=1+0.002*XY;
if ismember(Y_digit,[0 1 2 3])
    case_name='case_nongcun';
elseif ismember(Y_digit,[4 5 6])
    case_name='case_chengshi_A';
else
    case_name='case_kaifaqu_A';
end
if ismember(X_digit,[0 1 2])
    study_buses=[4 9 14]; model_type='A';
elseif ismember(X_digit,[3 4 5])
    study_buses=[5 11 16]; model_type='B';
else
    study_buses=[6 13 17]; model_type='C';
end
profiles=load(fullfile(lesson.cases,'profiles_24h.mat'));
T=profiles.T; dt=1;
load_curve_pu=profiles.load_curve_pu; pv_curve_pu=profiles.pv_curve_pu;
price_MWh=profiles.price_MWh; price_kWh=profiles.price_kWh;
disp(table(Stu_ID,string(case_name),alpha,string(model_type),string(mat2str(study_buses)),...
    'VariableNames',{'StudentID','Case','LoadMultiplier','Objective','StudyNodes'}));
mpc=loadcase(case_name);
cand_buses=study_buses;N_cand=numel(cand_buses);
yalmip('clear');
%%
demo_soc = 0.5; %[control:slider:demo_soc_1]{"position":[12,15]}
demo_w = 1; %[control:slider:demo_w_2]{"position":[10,11]}

Delta_Lower=max(0.2*demo_w,demo_soc-0.8*(1-demo_w));
Delta_Upper=min(0.8*demo_w,demo_soc-0.2*(1-demo_w));
disp(table(demo_soc,demo_w,Delta_Lower,Delta_Upper));
%%
P_opt_MW  = [1, 2, 3, 4, 5];
E_opt_MWh = [1, 2, 3, 4, 5];
K = length(E_opt_MWh);
eta_ch  = 0.90;
eta_dis = 0.95;
SoC_min = 0.20;
SoC_max = 0.80;
SoC_0   = 0.50;
c_P    = 0.08 * 1e6;
c_E    = 0.55 * 1e6;
c_om   = 200 * 1e3;
Y_life = 20;
d_rate = 0.05;
CRF    = (d_rate * (1 + d_rate)^Y_life) / ((1 + d_rate)^Y_life - 1);
D_days = 365;
V_min_Lin = 0.960;
V_max_Lin = 1.050;
%%
V_min_Lin = 0.95; %[control:slider:v_min_lin_3]{"position":[13,17]}
energy_cost_multiplier = 1; %[control:slider:energy_cost_multiplier_4]{"position":[26,29]}
c_E=0.55*1e6*energy_cost_multiplier;
%%
Sbase = mpc.baseMVA;
N_bus = size(mpc.bus, 1);
N_br  = size(mpc.branch, 1);
fbus  = mpc.branch(:, 1);
tbus  = mpc.branch(:, 2);
r_br  = mpc.branch(:, 3);
x_br  = mpc.branch(:, 4);
Pd_pu = (mpc.bus(:, 3) * alpha / Sbase) * load_curve_pu';
Qd_pu = (mpc.bus(:, 4) * alpha / Sbase) * load_curve_pu';
%%
z        = binvar(N_cand, 1, 'full');
w        = binvar(N_cand, K, 'full');
SoC      = sdpvar(N_cand, T, 'full');
delta    = sdpvar(N_cand, K, T, 'full');
P_ch_MW  = sdpvar(N_cand, T, 'full');
P_dis_MW = sdpvar(N_cand, T, 'full');
u_ch     = binvar(N_cand, T, 'full');
U_lin     = sdpvar(N_bus, T, 'full');
P_br      = sdpvar(N_br,  T, 'full');
Q_br      = sdpvar(N_br,  T, 'full');
P_grid_pu = sdpvar(1,     T, 'full');
Q_grid_pu = sdpvar(1,     T, 'full');
%%
%%
%%
Constraints = [];
P_rated     = w * P_opt_MW';
E_rated     = w * E_opt_MWh';
P_rated_MW  = P_rated;
E_rated_MWh = E_rated;
P_max_bound = max(P_opt_MW);
for idx = 1:N_cand
    Constraints = [Constraints, sum(w(idx, :)) == z(idx)];
    Constraints = [Constraints, SoC_min <= SoC(idx, :) <= SoC_max];
    Constraints = [Constraints, SoC(idx, T) == SoC_0];
    Constraints = [Constraints, 0 <= P_ch_MW(idx, :)  <= P_rated_MW(idx)];
    Constraints = [Constraints, 0 <= P_dis_MW(idx, :) <= P_rated_MW(idx)];
    Constraints = [Constraints, P_ch_MW(idx, :)  <= u_ch(idx, :) * P_max_bound];
    Constraints = [Constraints, P_dis_MW(idx, :) <= (1 - u_ch(idx, :)) * P_max_bound];
    for k = 1:K
        d_vec = reshape(delta(idx, k, :), 1, T);
        Constraints = [Constraints, SoC_min * w(idx, k) <= d_vec <= SoC_max * w(idx, k)];
        Constraints = [Constraints, SoC(idx, :) - SoC_max * (1 - w(idx, k)) <= d_vec <= ...
            SoC(idx, :) - SoC_min * (1 - w(idx, k))];
    end
    for t = 1:T
        if t == 1
            E_prev = E_opt_MWh * (SoC_0 * w(idx, :))';
        else
            E_prev = E_opt_MWh * reshape(delta(idx, :, t - 1), K, 1);
        end
        E_curr = E_opt_MWh * reshape(delta(idx, :, t), K, 1);
        Constraints = [Constraints, E_curr - E_prev == ...
            eta_ch * P_ch_MW(idx, t) - P_dis_MW(idx, t) / eta_dis];
    end
end
%%
A_from  = sparse(fbus, 1:N_br, 1, N_bus, N_br);
A_to    = sparse(tbus, 1:N_br, 1, N_bus, N_br);
A_cand  = sparse(cand_buses, 1:N_cand, 1, N_bus, N_cand);
A_slack = sparse(1, 1, 1, N_bus, 1);
Constraints = [Constraints, A_to * P_br - A_from * P_br + A_slack * P_grid_pu + ...
    A_cand * ((P_dis_MW - P_ch_MW) / Sbase) == Pd_pu];
Constraints = [Constraints, A_to * Q_br - A_from * Q_br + A_slack * Q_grid_pu == Qd_pu];
Constraints = [Constraints, U_lin(1, :) == 1.0];
Constraints = [Constraints, (V_min_Lin^2) <= U_lin(2:end, :) <= (V_max_Lin^2)];
for b = 1:N_br
    Constraints = [Constraints, U_lin(tbus(b), :) == U_lin(fbus(b), :) - ...
        2 * (r_br(b) * P_br(b, :) + x_br(b) * Q_br(b, :))];
end
%%
C_inv_ann  = CRF * sum(c_P * P_rated_MW + c_E * E_rated_MWh);
C_om_ann   = sum(c_om * P_rated_MW);
C_grid_ann = D_days * sum((Sbase * P_grid_pu') .* price_MWh);
Objective  = C_inv_ann + C_om_ann + C_grid_ann;
[sol,solver_name,solve_seconds]=ess_live_solve(Constraints,Objective,'MILP');
assert(sol.problem==0,'规划求解失败：%s',yalmiperror(sol.problem));
%%
V_lin_pu      = sqrt(value(U_lin));
V_ac_base     = zeros(N_bus, T);
V_ac_ess      = zeros(N_bus, T);
Ploss_ac_base = zeros(T, 1);
Ploss_ac_ess  = zeros(T, 1);
P_net_ess_MW  = value(P_dis_MW - P_ch_MW);
mpopt = mpoption('verbose', 0, 'out.all', 0);
for t = 1:T
    mpc_b = mpc;
    mpc_b.bus(:, 3) = Pd_pu(:, t) * Sbase;
    mpc_b.bus(:, 4) = Qd_pu(:, t) * Sbase;
    res_b = runpf(mpc_b, mpopt);
    assert(res_b.success==1,'无储能方案第 %d 时段潮流未收敛。',t);
    V_ac_base(:, t)  = res_b.bus(:, 8);
    Ploss_ac_base(t) = sum(real(get_losses(res_b)));
    mpc_e = mpc_b;
    mpc_e.bus(cand_buses, 3) = mpc_e.bus(cand_buses, 3) - P_net_ess_MW(:, t);
    res_e = runpf(mpc_e, mpopt);
    assert(res_e.success==1,'配置后第 %d 时段潮流未收敛。',t);
    V_ac_ess(:, t)  = res_e.bus(:, 8);
    Ploss_ac_ess(t) = sum(real(get_losses(res_e)));
end
B_loss_ann           = D_days * sum((Ploss_ac_base - Ploss_ac_ess) .* price_MWh);
max_voltage_error_pu = max(max(abs(V_lin_pu - V_ac_ess)));
fprintf('\n================ 第二周：储能选址定容与交流潮流验证报告 ================\n');
for idx = 1:N_cand
    fprintf('候选节点 %2d: 选址状态 z=%d | 额定功率 = %d MW | 额定容量 = %d MWh\n', ...
        cand_buses(idx), round(value(z(idx))), round(value(P_rated_MW(idx))), round(value(E_rated_MWh(idx))));
end
fprintf('年化投资成本: %.2f 万元 | 年运维成本: %.2f 万元 | 年网损降低收益: %.2f 万元\n', ...
    value(C_inv_ann) / 1e4, value(C_om_ann) / 1e4, B_loss_ann / 1e4);
fprintf('LinDistFlow 最低电压: %.4f p.u. | 交流潮流最低电压: %.4f p.u. | 最大电压偏差: %.4f p.u.\n', ...
    min(V_lin_pu(:)), min(V_ac_ess(:)), max_voltage_error_pu);
%%
Configuration=table(cand_buses(:),round(value(z)),value(P_rated_MW),value(E_rated_MWh),...
 'VariableNames',{'Node','Installed','RatedPower_MW','RatedEnergy_MWh'});
disp(Configuration);
AC_MinVoltage=min(V_ac_ess(:));AC_MaxVoltage=max(V_ac_ess(:));
disp(table(AC_MinVoltage,AC_MaxVoltage,max_voltage_error_pu,B_loss_ann));
%%

save(fullfile(report_out,"Planning_v095.mat"),"Configuration","V_lin_pu","V_ac_base","V_ac_ess","Ploss_ac_base","Ploss_ac_ess","B_loss_ann","max_voltage_error_pu","P_net_ess_MW","price_MWh","c_E","V_min_Lin","sol");
Plan=struct("z",value(z),"w",value(w),"Prated",value(P_rated),"Erated",value(E_rated),"SoC",value(SoC),"Pch",value(P_ch_MW),"Pdis",value(P_dis_MW),"Objective",value(Objective),"Cinv",value(C_inv_ann),"Com",value(C_om_ann),"Cgrid",value(C_grid_ann),"ConstraintMin",min(check(Constraints)));save(fullfile(report_out,"Planning_v095.mat"),"Plan","-append");
