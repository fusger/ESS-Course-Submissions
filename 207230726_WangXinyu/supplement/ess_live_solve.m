function [sol,solver_name,solve_seconds] = ess_live_solve(constraints,objective,model_type)
% LP/MILP/SOCP 使用 MATLAB 内置求解器；MISOCP 使用支持该模型的求解器。
model_type=upper(char(model_type));
if ismember(model_type,{'LP','MILP','SOCP'})
    original_path=path;
    restore_path=onCleanup(@()path(original_path)); %#ok<NASGU>
    % 优先使用 MathWorks 实现，避免同名的 MOSEK 兼容接口覆盖原生函数。
    native_dir=fullfile(matlabroot,'toolbox','optim','optim');
    addpath(native_dir,'-begin');
    switch model_type
        case 'LP',solver_name='linprog';
        case 'MILP',solver_name='intlinprog';
        case 'SOCP',solver_name='coneprog';
    end
    solver_file=which(solver_name);
    assert(startsWith(solver_file,matlabroot,'IgnoreCase',true),...
        '请安装并启用 Optimization Toolbox 的 %s。',solver_name);
    ops=sdpsettings('solver',solver_name,'verbose',1,'savesolveroutput',1);
    if strcmp(model_type,'MILP')
        ops.intlinprog.RelativeGapTolerance=1e-6;
        ops.intlinprog.AbsoluteGapTolerance=1e-8;
    end
    tic;sol=optimize(constraints,objective,ops);solve_seconds=toc;
    sol.native_solver_file=solver_file;
    assert(sol.problem==0,'%s 求解失败：%s',solver_name,sol.info);
    return
end
assert(strcmp(model_type,'MISOCP'),'未定义的模型类型：%s',model_type);
solver_candidates={};
if exist('mosekopt','file')~=0,solver_candidates{end+1}='mosek';end
if exist('gurobi','file')~=0,solver_candidates{end+1}='gurobi';end
assert(~isempty(solver_candidates),...
    '本任务含充放电二进制变量和二阶锥约束，需配置 MOSEK 或 Gurobi。');
messages=cell(0,1);tic;
objective_coefficients=getbase(objective);
if max(abs(objective_coefficients(:)))>1000,objective=objective/1000;end
for k=1:numel(solver_candidates)
    solver_name=solver_candidates{k};
    ops=sdpsettings('solver',solver_name,'verbose',1,'savesolveroutput',1,'warmstart',1);
    if strcmp(solver_name,'mosek')
        ops.mosek.MSK_DPAR_MIO_MAX_TIME=600;
        ops.mosek.MSK_IPAR_NUM_THREADS=1;
        ops.mosek.MSK_DPAR_MIO_TOL_REL_GAP=1e-6;
        ops.mosek.MSK_DPAR_MIO_TOL_ABS_GAP=1e-8;
        ops.mosek.MSK_DPAR_MIO_TOL_FEAS=1e-9;
        ops.mosek.MSK_DPAR_MIO_TOL_ABS_RELAX_INT=1e-9;
        ops.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-10;
        ops.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-10;
        ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-10;
    else
        ops.gurobi.TimeLimit=120;
        ops.gurobi.MIPGap=1e-6;
        ops.gurobi.FeasibilityTol=1e-9;
        ops.gurobi.IntFeasTol=1e-9;
        ops.gurobi.OptimalityTol=1e-8;
        ops.gurobi.BarQCPConvTol=1e-9;
    end
    sol=optimize(constraints,objective,ops);
    save(fullfile(fileparts(mfilename('fullpath')),'last_solver_diagnostic.mat'),'sol','ops');
    if sol.problem==0,solve_seconds=toc;return;end
    messages{end+1}=[solver_name ': ' sol.info]; %#ok<AGROW>
end
error('MISOCP 求解失败：%s',strjoin(messages,'; '));
end

