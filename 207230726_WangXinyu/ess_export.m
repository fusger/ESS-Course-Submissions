function ess_export(fig,name)
% 白底 600 dpi PNG、矢量 PDF 及可编辑 FIG。
set(fig,'Color','w');
exportgraphics(fig,[name '.png'],'Resolution',600,'BackgroundColor','white');
exportgraphics(fig,[name '.pdf'],'ContentType','vector','BackgroundColor','white');
savefig(fig,[name '.fig']);
end
