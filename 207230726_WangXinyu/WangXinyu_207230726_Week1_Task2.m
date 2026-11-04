%% 任务2：单节点分时电价套利
% 汪昕宇 207230726 储能231；根据教师示例实现个人参数计算。
clear;clc;close all;
ess_setup(fileparts(which('ess_setup')));
Stu_ID=207230726;
S=load('profiles_24h.mat');T=S.T;
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
f=figure('Visible','off','Position',[100 100 860 760]);
tiledlayout(3,1,'TileSpacing','compact');
nexttile;stairs(1:T,S.price_MWh,'LineWidth',1.8);ess_axes;ylabel('电价 / 元每MWh');title('2 MW / 2 MWh 单节点储能套利');
nexttile;plot(1:T,[Arbitrage.MILP.P_ch Arbitrage.MILP.P_dis],'LineWidth',1.8);ess_axes;
ylabel('功率 / MW');legend('充电','放电','Location','best');
nexttile;plot(0:T,[0.5;Arbitrage.MILP.SoC],'LineWidth',1.8);ess_axes;xlim([0 24]);
ylabel('荷电状态 SoC');xlabel('时段 / h');
ess_export(f,sprintf('Fig_Arbitrage_%d',Stu_ID));close(f);
save('Task2_details.mat','Arbitrage');disp(Arbitrage);
for ext={'png','pdf','fig'},copyfile(['Fig_Arbitrage_207230726.' ext{1}],['Fig_Week1_Task2_207230726.' ext{1}]);end
