S=load('profiles_24h.mat');mpc=loadcase('case_chengshi_A');
AlphaTrials=cell(2,1);
for aidx=1:2
 a=[1.052 1.15];a=a(aidx);V=zeros(23,24);Loss=zeros(24,1);Grid=Loss;
 for t=1:24
  b=mpc;b.bus(:,3)=mpc.bus(:,3)*a*S.load_curve_pu(t);b.bus(:,4)=mpc.bus(:,4)*a*S.load_curve_pu(t);
  b.bus([4 9 14],3)=b.bus([4 9 14],3)-.8*S.pv_curve_pu(t);
  rr=runpf(b,mpoption('verbose',0,'out.all',0));assert(rr.success==1);V(:,t)=rr.bus(:,8);Loss(t)=sum(real(get_losses(rr)));Grid(t)=sum(rr.gen(:,2));
 end
 AlphaTrials{aidx}=struct('Alpha',a,'V',V,'Loss',Loss,'Grid',Grid);
end
ArbitrageTrials=cell(5,1);configs=[1 2;2 2;4 2;2 1;2 4];
for aidx=1:size(configs,1)
 Pr=configs(aidx,1);Er=configs(aidx,2);T=24;yalmip('clear');c=sdpvar(T,1);d=sdpvar(T,1);soc=sdpvar(T,1);u=binvar(T,1);
 C=[0<=c<=Pr,0<=d<=Pr,.1<=soc<=1,soc(T)==.5,c<=Pr*u,d<=Pr*(1-u)];
 for t=1:T
  if t==1,prev=.5;else,prev=soc(t-1);end
  C=[C,soc(t)==prev+(.9*c(t)-d(t)/.95)/Er];
 end
 F=sum(S.price_MWh.*(c-d));rr=optimize(C,F,sdpsettings('solver','intlinprog','verbose',0));assert(rr.problem==0);
 ArbitrageTrials{aidx}=struct('Prated',Pr,'Erated',Er,'Revenue',-value(F),'Pch',value(c),'Pdis',value(d),'SoC',value(soc),'Residual',min(check(C)));
end
save(fullfile(report_out,'Supplement_Basics.mat'),'AlphaTrials','ArbitrageTrials');
disp('BASICS_SUPPLEMENT_DONE');

