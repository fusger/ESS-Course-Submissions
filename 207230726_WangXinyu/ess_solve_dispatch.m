function [sol,solver_name,info]=ess_solve_dispatch(C,F,u,c,d) %#ok<INUSD>
% 完整 MISOCP：优先课程指定求解器，备用 SCIP，不删除互斥约束。
tic;names={};
if exist('mosekopt','file'),names{end+1}='mosek';end
if exist('gurobi','file'),names{end+1}='gurobi';end
if strcmp(getenv('ESS_REQUIRE_MOSEK'),'1')
 assert(exist('mosekopt','file')~=0,'未找到 MOSEK MATLAB 接口');names={'mosek'};
end
for k=1:numel(names)
 name=names{k};ops=sdpsettings('solver',name,'verbose',1,'savesolveroutput',1);
 if strcmp(name,'gurobi'),ops.gurobi.MIPGap=1e-6;ops.gurobi.FeasibilityTol=1e-9;end
 if strcmp(name,'mosek')
  ops.mosek.MSK_DPAR_MIO_TOL_REL_GAP=1e-6;
  ops.mosek.MSK_DPAR_MIO_TOL_ABS_GAP=1e-8;
  ops.mosek.MSK_DPAR_MIO_TOL_FEAS=1e-9;
  ops.mosek.MSK_DPAR_MIO_TOL_ABS_RELAX_INT=1e-9;
  ops.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-10;
  ops.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-10;
  ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-10;
 end
 % 沿用课程示范的官方YALMIP-MOSEK接口及整数误差容差。
 sol=optimize(C,F,ops);
 if sol.problem==0
  b=getbase(F);
  solver_name=name;info=struct('Seconds',toc,'Output',sol.solveroutput,...
   'ObjectiveScale',1,'ObjectiveOffset',full(b(1)));
  if strcmp(name,'mosek')
   [~,vr]=mosekopt('version echo(0)');
   info.Version=sprintf('%d.%d.%d',vr.version.major,vr.version.minor,vr.version.revision);
  end
  return
 end
 warning('%s: %s',name,sol.info);
end
assert(~strcmp(getenv('ESS_REQUIRE_MOSEK'),'1'),'本次要求使用 MOSEK，未成功求解，停止而不切换备用求解器。');
% 仅借用 coneprog 格式导出器；二进制列显式传给 SCIP。
ops=sdpsettings('solver','coneprog','relax',1,'verbose',0);
[~,recoverdata,diagnostic,p]=export(C,F,ops);
assert(isempty(diagnostic)||diagnostic.problem==0,'YALMIP 导出失败');
F_struc=p.F_struc;lb=p.lb;ub=p.ub;nf=p.K.f;nl=p.K.l;q=p.K.q;
[ok,bins]=ismember(getvariables(u),recoverdata.used_variables);
assert(all(ok)&&numel(bins)==72,'72 个互斥二进制变量导出失败');
c=p.c;objective_offset=p.f;
save('dispatch_solver_input.mat','F_struc','lb','ub','nf','nl','q','c','bins','objective_offset','-v7');
python=getenv('ESS_PYTHON');
if isempty(python)
 config=fullfile(fileparts(mfilename('fullpath')),'python_path.local.txt');
 if isfile(config)
  python=strtrim(fileread(config));
 else
  candidate=fullfile(getenv('USERPROFILE'),'.cache','codex-runtimes','codex-primary-runtime','dependencies','python','python.exe');
  if isfile(candidate),python=candidate;else,python='python';end
 end
end
script=fullfile(fileparts(mfilename('fullpath')),'ess_scip.py');
cmd=sprintf('"%s" "%s" "%s" "%s"',python,script,...
 fullfile(pwd,'dispatch_solver_input.mat'),fullfile(pwd,'dispatch_solver_output.mat'));
[status,msg]=system(cmd);fprintf('%s\n',msg);
assert(status==0,'SCIP 求解未完成，请按使用说明安装 Python 求解依赖。');
out=load('dispatch_solver_output.mat');
assign(recover(recoverdata.used_variables),recoverdata.x_equ+recoverdata.H*out.x(:));
assert(min(check(C))>=-1e-6,'SCIP 回代原始约束失败');
solver_name='SCIP (PySCIPOpt, full MISOCP)';info=struct('Seconds',toc,'Output',out);
sol=struct('problem',0,'info','SCIP full MISOCP solved and checked');
end
