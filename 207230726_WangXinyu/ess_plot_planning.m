function ess_plot_planning(P)
[~,worst]=min(min(P.V_ac_base,[],2));
f=figure('Visible','off','Position',[100 100 860 650]);
tiledlayout(2,1,'TileSpacing','compact');
nexttile;plot(1:24,P.V_ac_base(worst,:),'k--',1:24,P.V_lin(worst,:),'b-^',...
 1:24,P.V_ac_ess(worst,:),'r-o','LineWidth',1.6,'MarkerSize',4);
ess_axes;yline(0.95,':r');ylabel('电压幅值 / p.u.');
title(sprintf('节点 %d 规划前后电压与线性化验证',worst));
legend('未配置储能 AC','配置储能 LinDistFlow','配置储能 AC','电压下限','Location','best');
nexttile;plot(1:24,P.Loss_base_MW,'k--',1:24,P.Loss_ess_MW,'b-','LineWidth',1.8);
ess_axes;ylabel('全网有功损耗 / MW');xlabel('时段 / h');legend('未配置储能','配置储能','Location','best');
ess_export(f,sprintf('Fig_Week2_Planning_%d',P.StuID));close(f);
f=figure('Visible','off','Position',[100 100 860 650]);
tiledlayout(2,1,'TileSpacing','compact');
nexttile;plot(1:24,(P.P_dis_MW-P.P_ch_MW)','LineWidth',1.8);ess_axes;
ylabel('净放电功率 / MW');legend(arrayfun(@(b)sprintf('节点 %d',b),P.Buses,'UniformOutput',false),'Location','best');
title('选址定容最优方案的逐时运行');
nexttile;plot(0:24,[0.5*ones(3,1) P.SoC]','LineWidth',1.8);ess_axes;xlim([0 24]);
xlabel('时段 / h');ylabel('荷电状态 SoC');ylim([0.15 0.85]);
ess_export(f,sprintf('Fig_Week2_Operation_%d',P.StuID));close(f);
end
