%% 第二周：配电网储能选址定容优化 (Big-M 线性化 + LinDistFlow) 与交流潮流验证
% 课程名称：《储能优化配置与调度技术》课程设计 | 作者：汪昕宇（学号：207230726，储能231）

%% 0. 核心数学模型
% *1) 全生命周期等年值总成本目标函数*：
% 
% $$\min F = C_{\mathrm{inv}}^{\mathrm{ann}} + C_{\mathrm{om}}^{\mathrm{ann}} + C_{\mathrm{grid}}^{\mathrm{ann}} = \mathrm{CRF} \sum_{i\in\Omega_{\mathrm{cand}}} \left( c_P P_i^{\mathrm{rated}} + c_E E_i^{\mathrm{rated}} \right) + \sum_{i\in\Omega_{\mathrm{cand}}} c_{\mathrm{om}} P_i^{\mathrm{rated}} + D_{\mathrm{days}} \sum_{t=1}^{T} \lambda_t P_{\mathrm{grid},t} \Delta t$$
% 
% 其中资本回收系数（Capital Recovery Factor, CRF）为：
% 
% $$\mathrm{CRF} = \frac{d(1+d)^Y}{(1+d)^Y - 1}$$
% 
% *2) 双线性非凸项* $\mathrm{SoC}_{i,t} \cdot E_i^{\mathrm{rated}}$ *的离散档位编码与紧致 Big-M（McCormick 包络）线性化*：
% 
% 设离散容量备选档位集合为 $E^{\mathrm{opt}} = [E_1^{\mathrm{opt}}, \dots, E_K^{\mathrm{opt}}]$，引入二进制档位选择变量 $w_{i,k} \in \{0,1\}$ 与连续辅助变量 $\delta_{i,k,t} = \mathrm{SoC}_{i,t} \cdot w_{i,k}$，满足：
% 
% $$\sum_{k=1}^{K} w_{i,k} = z_i, \quad E_i^{\mathrm{rated}} = \sum_{k=1}^{K} w_{i,k} E_k^{\mathrm{opt}}, \quad P_i^{\mathrm{rated}} = \sum_{k=1}^{K} w_{i,k} P_k^{\mathrm{opt}}$$
% 
% $$\mathrm{SoC}_{\min} \cdot w_{i,k} \le \delta_{i,k,t} \le \mathrm{SoC}_{\max} \cdot w_{i,k}$$
% 
% $$\mathrm{SoC}_{i,t} - \mathrm{SoC}_{\max}\left(1 - w_{i,k}\right) \le \delta_{i,k,t} \le \mathrm{SoC}_{i,t} - \mathrm{SoC}_{\min}\left(1 - w_{i,k}\right)$$
% 
% 线性化后的储能时序能量递推方程为：
% 
% $$\sum_{k=1}^{K} E_k^{\mathrm{opt}} \left( \delta_{i,k,t} - \delta_{i,k,t-1} \right) = \left( \eta_{\mathrm{ch}} P_{i,t}^{\mathrm{ch}} - \frac{1}{\eta_{\mathrm{dis}}} P_{i,t}^{\mathrm{dis}} \right) \Delta t$$
% 
% *3) 配电网 LinDistFlow 线性化潮流方程与工程电压安全裕度约束*：
% 
% $$U_{j,t}^{\mathrm{Lin}} = U_{i,t}^{\mathrm{Lin}} - 2\left( r_{ij} P_{ij,t} + x_{ij} Q_{ij,t} \right), \quad \left(V_{\min}^{\mathrm{Lin}}\right)^2 \le U_{i,t}^{\mathrm{Lin}} \le \left(V_{\max}^{\mathrm{Lin}}\right)^2$$

clear; clc; close all;
ess_setup(fileparts(which('ess_setup')));
yalmip('clear');

%% 2. 学号哈希映射与网络参数读取
% 根据个人 9 位学号自动映射配电网算例、负荷缩放系数 $\alpha = 1 + 0.002 \times \mathrm{XY}$ 与储能候选安装节点集合 $\Omega_{\mathrm{cand}}$。
Stu_ID  = 207230726;
XY      = mod(Stu_ID, 100);
X_digit = floor(XY / 10);
Y_digit = mod(XY, 10);
alpha   = 1 + 0.002 * XY;

if ismember(Y_digit, [0, 1, 2, 3])
    case_name = 'case_nongcun';
elseif ismember(Y_digit, [4, 5, 6])
    case_name = 'case_chengshi_A';
else
    case_name = 'case_kaifaqu_A';
end
mpc = loadcase(case_name);
load('profiles_24h.mat');

if ismember(X_digit, [0, 1, 2])
    cand_buses = [4, 9, 14];
elseif ismember(X_digit, [3, 4, 5])
    cand_buses = [5, 11, 16];
else
    cand_buses = [6, 13, 17];
end
N_cand = length(cand_buses);

%% 3. 规划经济参数与储能离散容量档位
% 为补偿 LinDistFlow 忽略支路电流平方项 $L_{ij,t}$ 导致的节点电压高估偏差，线性规划内部将电压下限保守收紧为 $V_{\min}^{\mathrm{Lin}} = 0.960\;\mathrm{p.u.}$，确保回代交流潮流后真实电压满足 $V_{i,t}^{\mathrm{AC}} \ge 0.950\;\mathrm{p.u.}$。
P_opt_MW  = [1, 2, 3, 4, 5]; % 离散功率档位 P_k^{opt} (MW)
E_opt_MWh = [1, 2, 3, 4, 5]; % 离散能量档位 E_k^{opt} (MWh)
K = length(E_opt_MWh);

eta_ch  = 0.90;
eta_dis = 0.95;
SoC_min = 0.20;
SoC_max = 0.80;
SoC_0   = 0.50;

c_P    = 0.08 * 1e6;   % 单位功率投资成本: 80,000 元/MW
c_E    = 0.55 * 1e6;   % 单位能量投资成本: 550,000 元/MWh
c_om   = 200 * 1e3;    % 单位功率年运维成本: 200,000 元/(MW*年)
Y_life = 20;           % 规划年限: 20 年
d_rate = 0.05;         % 贴现率: 5%
CRF    = (d_rate * (1 + d_rate)^Y_life) / ((1 + d_rate)^Y_life - 1);
D_days = 365;

V_min_Lin = 0.960;     % LinDistFlow 节点电压安全下限 (p.u.)
V_max_Lin = 1.050;     % LinDistFlow 节点电压安全上限 (p.u.)

%% 4. 配电网拓扑与负荷标幺值计算
Sbase = mpc.baseMVA;
N_bus = size(mpc.bus, 1);
N_br  = size(mpc.branch, 1);
fbus  = mpc.branch(:, 1);
tbus  = mpc.branch(:, 2);
r_br  = mpc.branch(:, 3);
x_br  = mpc.branch(:, 4);

Pd_pu = (mpc.bus(:, 3) * alpha / Sbase) * load_curve_pu';
Qd_pu = (mpc.bus(:, 4) * alpha / Sbase) * load_curve_pu';

%% 5. 声明 YALMIP 决策变量
% 保留 $\mathrm{SoC}_{i,t}$ 为运行决策变量，引入档位选择二进制变量 $w_{i,k} \in \{0,1\}$ 与辅助变量 $\delta_{i,k,t} = \mathrm{SoC}_{i,t} \cdot w_{i,k}$。
z        = binvar(N_cand, 1, 'full');    % 候选节点选址二进制变量 z_i
w        = binvar(N_cand, K, 'full');    % 容量档位选择二进制变量 w_{i,k}
SoC      = sdpvar(N_cand, T, 'full');    % 储能荷电状态 SoC_{i,t}
delta    = sdpvar(N_cand, K, T, 'full'); % Big-M 辅助变量 \delta_{i,k,t} = SoC_{i,t} * w_{i,k}
P_ch_MW  = sdpvar(N_cand, T, 'full');    % 充电有功功率 P_{i,t}^{ch} (MW)
P_dis_MW = sdpvar(N_cand, T, 'full');    % 放电有功功率 P_{i,t}^{dis} (MW)
u_ch     = binvar(N_cand, T, 'full');    % 充放电互斥二进制变量 u_{i,t}^{ch}

U_lin     = sdpvar(N_bus, T, 'full');    % LinDistFlow 节点电压平方 U_{i,t}^{Lin} (p.u.^2)
P_br      = sdpvar(N_br,  T, 'full');    % 支路有功潮流 P_{ij,t} (p.u.)
Q_br      = sdpvar(N_br,  T, 'full');    % 支路无功潮流 Q_{ij,t} (p.u.)
P_grid_pu = sdpvar(1,     T, 'full');    % 关口购电有功功率 P_{grid,t} (p.u.)
Q_grid_pu = sdpvar(1,     T, 'full');    % 关口注入无功功率 Q_{grid,t} (p.u.)

%% 6. 构建储能选址定容与 Big-M 线性化约束
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

%% 7. 构建配电网 LinDistFlow 线性潮流约束
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

%% 8. 定义年综合费用目标函数并调用 MATLAB 原生 MILP 求解器 (intlinprog)
% 任务四选址定容模型为混合整数线性规划（MILP），默认调用 MATLAB Optimization Toolbox 原生 |intlinprog| 求解。为防止路径中 MOSEK 的同名兼容函数 |intlinprog.m| 覆盖 MathWorks 原生函数，求解前将内置 |optim| 目录置顶。
C_inv_ann  = CRF * sum(c_P * P_rated_MW + c_E * E_rated_MWh);
C_om_ann   = sum(c_om * P_rated_MW);
C_grid_ann = D_days * sum((Sbase * P_grid_pu') .* price_MWh);
Objective  = C_inv_ann + C_om_ann + C_grid_ann;

native_optim_dir = fullfile(matlabroot, 'toolbox', 'optim', 'optim');
if exist(native_optim_dir, 'dir')
    addpath(native_optim_dir, '-begin');
end
solver_name = 'intlinprog';
ops = sdpsettings('solver', solver_name, 'verbose', 1);
ops.intlinprog.RelativeGapTolerance = 1e-6;
ops.intlinprog.AbsoluteGapTolerance = 1e-8;
sol = optimize(Constraints, Objective, ops);

if sol.problem ~= 0
    error('MILP 选址定容求解失败: %s', yalmiperror(sol.problem));
end

%% 9. 基于 MATPOWER runpf 的事后交流潮流闭环验证与真实降损收益核算
% 将 MILP 规划求得的各节点净放电出力 $P_{i,t}^{\mathrm{net}*} = P_{i,t}^{\mathrm{dis}*} - P_{i,t}^{\mathrm{ch}*}$ 代回 MATPOWER |runpf| 进行 24 小时精确牛顿-拉夫逊交流潮流回代校验，核算真实年化降损收益 $B_{\mathrm{loss}}^{\mathrm{AC}}$ 与最大线性化电压偏差 $\Delta V_{\max}$：
% 
% $$B_{\mathrm{loss}}^{\mathrm{AC}} = D_{\mathrm{days}} \sum_{t=1}^{T} \left( P_{\mathrm{loss},t}^{\mathrm{AC,base}} - P_{\mathrm{loss},t}^{\mathrm{AC,ess}} \right) \lambda_t \Delta t, \quad \Delta V_{\max} = \max_{i\in\Omega_{\mathrm{N}}, t\in\mathcal{T}} \left| \sqrt{U_{i,t}^{\mathrm{Lin}}} - V_{i,t}^{\mathrm{AC,ess}} \right|$$
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
    assert(res_b.success==1,'基态交流潮流第 %d 小时未收敛',t);
    V_ac_base(:, t)  = res_b.bus(:, 8);
    Ploss_ac_base(t) = sum(real(get_losses(res_b)));
    
    mpc_e = mpc_b;
    mpc_e.bus(cand_buses, 3) = mpc_e.bus(cand_buses, 3) - P_net_ess_MW(:, t);
    res_e = runpf(mpc_e, mpopt);
    assert(res_e.success==1,'规划交流潮流第 %d 小时未收敛',t);
    V_ac_ess(:, t)  = res_e.bus(:, 8);
    Ploss_ac_ess(t) = sum(real(get_losses(res_e)));
end

B_loss_ann           = D_days * sum((Ploss_ac_base - Ploss_ac_ess) .* price_MWh);
max_voltage_error_pu = max(max(abs(V_lin_pu - V_ac_ess)));

fprintf('\n================ 第二周：储能选址定容与交流潮流验证报告 ================\n');
fprintf('学号: %d | 算例: %s | 求解器: %s\n', Stu_ID, case_name, solver_name);
for idx = 1:N_cand
    fprintf('候选节点 %2d: 选址状态 z=%d | 额定功率 = %d MW | 额定容量 = %d MWh\n', ...
        cand_buses(idx), round(value(z(idx))), round(value(P_rated_MW(idx))), round(value(E_rated_MWh(idx))));
end
fprintf('年化投资成本: %.2f 万元 | 年运维成本: %.2f 万元 | 事后核算年降损收益: %.2f 万元\n', ...
    value(C_inv_ann) / 1e4, value(C_om_ann) / 1e4, B_loss_ann / 1e4);
fprintf('LinDistFlow 最低电压: %.4f p.u. | 交流潮流真实最低电压: %.4f p.u. | 最大电压偏差: %.4f p.u.\n', ...
    min(V_lin_pu(:)), min(V_ac_ess(:)), max_voltage_error_pu);
fprintf('========================================================================\n');

%% 10. 预埋标准化输出结构体与《中国电机工程学报》规范绘图
% 注意：|Result.Opt_Sites| 保存的是候选节点集合 |cand_buses| 中的序号（|1..3|），实际配电网节点编号为 |cand_buses(Result.Opt_Sites)|。
Result.StuID       = Stu_ID;                   % 学生学号
Result.Obj_Value   = value(Objective);         % 年综合费用最优值 (元)
Result.Opt_Sites   = find(round(value(z))==1); % 选址候选节点集合中的序号 (1..3，对应母线 cand_buses(Result.Opt_Sites))
Result.Opt_P_rated = value(P_rated);           % 各候选节点配置功率 (MW)
Result.Opt_E_rated = value(E_rated);           % 各候选节点配置能量 (MWh)
Result.Max_V_Err   = max_voltage_error_pu;     % LinDistFlow 与 AC 潮流最大电压误差 (p.u.)
save('result_summary.mat', 'Result');
fprintf('[OK] 标准化评测结果已自动保存至当前目录的 result_summary.mat\n');

%% 11. 保存完整规划结果并检验物理约束
Planning=struct('StuID',Stu_ID,'Case',case_name,'Alpha',alpha,'Buses',cand_buses,...
    'Solver',solver_name,'Objective',value(Objective),'Sites',round(value(z)),...
    'P_rated',value(P_rated),'E_rated',value(E_rated),'w',value(w),'SoC',value(SoC),...
    'delta',value(delta),'P_ch_MW',value(P_ch_MW),'P_dis_MW',value(P_dis_MW),...
    'U',value(U_lin),'P_br',value(P_br),'Q_br',value(Q_br),...
    'V_lin',V_lin_pu,'V_ac_base',V_ac_base,'V_ac_ess',V_ac_ess,...
    'Loss_base_MW',Ploss_ac_base,'Loss_ess_MW',Ploss_ac_ess,...
    'C_inv_ann',value(C_inv_ann),'C_om_ann',value(C_om_ann),'C_grid_ann',value(C_grid_ann),...
    'B_loss_ann',B_loss_ann,'Max_V_Err',max_voltage_error_pu,...
    'Price',price_MWh,'ConstraintMin',min(check(Constraints)));
assert(Planning.ConstraintMin>=-1e-6,'规划约束残差过大');
assert(max(min(Planning.P_ch_MW(:),Planning.P_dis_MW(:)))<1e-5,'存在同时充放电');
assert(min(V_ac_ess(:))>=0.95-1e-5 && max(V_ac_ess(:))<=1.05+1e-5,'规划交流电压越限');
save('Week2_details.mat','Planning');
ess_plot_planning(Planning);

for ext={'png','pdf','fig'},copyfile(['Fig_Week2_Planning_207230726.' ext{1}],['Fig_Week2_207230726.' ext{1}]);end
