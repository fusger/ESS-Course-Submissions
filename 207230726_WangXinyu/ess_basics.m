function ess_basics(Stu_ID)
% 第一周周一与周二任务：基态潮流和单节点套利。
S=load('profiles_24h.mat');mpc=loadcase('case_chengshi_A');alpha=1+0.002*mod(Stu_ID,100);
Buses=[4 9 14];T=S.T;N=size(mpc.bus,1);
V=zeros(N,T);Loss=zeros(1,T);Grid=zeros(1,T);
% 基态光伏 0.8 MW/节点，另以 1.2 MW/节点作调度的同场景无储能对照。
V12=zeros(N,T);Loss12=zeros(1,T);Grid12=zeros(1,T);
for t=1:T
 b=mpc;b.bus(:,3)=mpc.bus(:,3)*alpha*S.load_curve_pu(t);
 b.bus(:,4)=mpc.bus(:,4)*alpha*S.load_curve_pu(t);
 b.bus(Buses,3)=b.bus(Buses,3)-0.8*S.pv_curve_pu(t);
 r=runpf(b,mpoption('verbose',0,'out.all',0));assert(r.success==1,'基态潮流失败');
 V(:,t)=r.bus(:,8);Loss(t)=sum(real(get_losses(r)));Grid(t)=sum(r.gen(:,2));
 b.bus(Buses,3)=b.bus(Buses,3)-0.4*S.pv_curve_pu(t);
 r=runpf(b,mpoption('verbose',0,'out.all',0));assert(r.success==1,'调度对照潮流失败');
 V12(:,t)=r.bus(:,8);Loss12(t)=sum(real(get_losses(r)));Grid12(t)=sum(r.gen(:,2));
end
[v,ix]=min(V(:));[node,hour]=ind2sub(size(V),ix);
Baseline=struct('StuID',Stu_ID,'V',V,'Loss_MW',Loss,'Grid_MW',Grid,...
 'MinV',v,'WorstNode',node,'WorstHour',hour,'DailyLoss_MWh',sum(Loss),...
 'V12',V12,'Loss12_MW',Loss12,'Grid12_MW',Grid12,...
 'Alpha',alpha,'Price',S.price_MWh,'Load',S.load_curve_pu,'PV',S.pv_curve_pu);
% 同参数 LP 和 MILP，保留每次求解状态和实际求解结果。
Arbitrage=struct;
for binary=[false true]
 yalmip('clear');c=sdpvar(T,1);d=sdpvar(T,1);soc=sdpvar(T,1);
 C=[0<=c<=2,0<=d<=2,0.1<=soc<=1,soc(T)==0.5];
 if binary,u=binvar(T,1);C=[C,c<=2*u,d<=2*(1-u)];end
 for t=1:T
  if t==1,prev=0.5;else,prev=soc(t-1);end
  C=[C,soc(t)==prev+(0.90*c(t)-d(t)/0.95)/2];
 end
 F=sum(S.price_MWh.*(c-d));
 if binary,name='intlinprog';tag='MILP';else,name='linprog';tag='LP';end
 result=optimize(C,F,sdpsettings('solver',name,'verbose',0));
 assert(result.problem==0,'套利求解失败：%s',result.info);
 A=struct('Solver',name,'Revenue',-value(F),'P_ch',value(c),'P_dis',value(d),...
  'SoC',value(soc),'ConstraintMin',min(check(C)));
 assert(A.ConstraintMin>=-1e-6,'套利约束残差过大');Arbitrage.(tag)=A;
end
save('Basics_details.mat','Baseline','Arbitrage');
f=figure('Visible','off','Position',[100 100 860 650]);
tiledlayout(2,1,'TileSpacing','compact');
nexttile;plot(1:T,V(node,:),'b-o','LineWidth',1.8,'MarkerSize',4);ess_axes;
yline(0.95,':r');ylabel('电压 / p.u.');title(sprintf('0.8 MW 光伏基态薄弱节点 %d',node));
nexttile;plot(1:T,Loss,'LineWidth',1.8);ess_axes;xlabel('时段 / h');ylabel('有功损耗 / MW');
ess_export(f,sprintf('Fig_Baseline_%d',Stu_ID));close(f);
f=figure('Visible','off','Position',[100 100 860 760]);
tiledlayout(3,1,'TileSpacing','compact');
nexttile;stairs(1:T,S.price_MWh,'LineWidth',1.8);ess_axes;ylabel('电价 / 元每MWh');title('2 MW / 2 MWh 单节点储能套利');
nexttile;plot(1:T,[Arbitrage.MILP.P_ch Arbitrage.MILP.P_dis],'LineWidth',1.8);ess_axes;
ylabel('功率 / MW');legend('充电','放电','Location','best');
nexttile;plot(0:T,[0.5;Arbitrage.MILP.SoC],'LineWidth',1.8);ess_axes;xlim([0 24]);
ylabel('荷电状态 SoC');xlabel('时段 / h');
ess_export(f,sprintf('Fig_Arbitrage_%d',Stu_ID));close(f);
fprintf('BASICS_OK MinV=%.6f; daily loss=%.6f; arbitrage LP/MILP=%.6f/%.6f\n',...
 v,sum(Loss),Arbitrage.LP.Revenue,Arbitrage.MILP.Revenue);
end
