%% 任务1：配电网24小时基态潮流
% 汪昕宇 207230726 储能231；根据教师示例实现个人参数计算。
clear;clc;close all;
ess_setup(fileparts(which('ess_setup')));
Stu_ID=207230726;
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
f=figure('Visible','off','Position',[100 100 860 650]);
tiledlayout(2,1,'TileSpacing','compact');
nexttile;plot(1:T,V(node,:),'b-o','LineWidth',1.8,'MarkerSize',4);ess_axes;
yline(0.95,':r');ylabel('电压 / p.u.');title(sprintf('0.8 MW 光伏基态薄弱节点 %d',node));
nexttile;plot(1:T,Loss,'LineWidth',1.8);ess_axes;xlabel('时段 / h');ylabel('有功损耗 / MW');
ess_export(f,sprintf('Fig_Baseline_%d',Stu_ID));close(f);
save('Task1_details.mat','Baseline');
for ext={'png','pdf','fig'},copyfile(['Fig_Baseline_207230726.' ext{1}],['Fig_Week1_Task1_207230726.' ext{1}]);end
