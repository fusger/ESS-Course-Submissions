function ess_plot_dispatch(D)
f=figure('Visible','off','Position',[100 100 860 760]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
labs=arrayfun(@(b)sprintf('节点 %d',b),D.Buses,'UniformOutput',false);
nexttile;plot(1:24,(D.P_dis_MW-D.P_ch_MW)','LineWidth',1.8);
ylabel('净放电功率 / MW');title('城市配电网模型 A 储能日前调度');ess_axes;legend(labs,'Location','eastoutside');
nexttile;plot(0:24,[0.5*ones(3,1) D.SoC]','LineWidth',1.8);
ylabel('荷电状态 SoC');ylim([0.15 0.85]);ess_axes;xlim([0 24]);
nexttile;plot(1:24,sum(D.PV_available_MW),'k--',1:24,sum(D.PV_MW),'b-',...
 1:24,D.P_grid_MW,'r-','LineWidth',1.8);ess_axes;
xlabel('时段 / h');ylabel('有功功率 / MW');legend('可用光伏','消纳光伏','关口购电','Location','eastoutside');
ess_export(f,sprintf('Fig_Week1_Dispatch_%d',D.StuID));close(f);
f=figure('Visible','off','Position',[100 100 860 650]);
tiledlayout(2,1,'TileSpacing','compact');
nexttile;plot(1:24,min(sqrt(D.U),[],1),'b-',1:24,min(D.AC_V,[],1),'r--','LineWidth',1.8);
ess_axes;ylabel('逐时最低电压 / p.u.');legend('MISOCP 模型','交流潮流','Location','best');
title(sprintf('调度电压验证 最大偏差 %.2e p.u.',D.AC_V_Error));
nexttile;imagesc(1:24,1:size(D.U,1),D.AC_V);axis xy;colorbar;
set(gca,'FontName','Microsoft YaHei');xlabel('时段 / h');ylabel('节点编号');
ess_export(f,sprintf('Fig_Week1_Voltage_%d',D.StuID));close(f);
end
