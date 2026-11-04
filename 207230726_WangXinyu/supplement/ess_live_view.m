function fig = ess_live_view(kind,varargin)
% 图形角色及色值统一读取本项目颜色表。
st=ess_live_style;
fig=figure('Color','w','Position',[100 100 st.defaults.width_px st.defaults.height_px]);
if isprop(fig,'Theme'),theme(fig,'light');end
switch kind
    case 'voltage'
        mpc=varargin{1};V=varargin{2};t=round(varargin{3});node=round(varargin{4});
        t=max(1,min(size(V,2),t));node=max(1,min(size(V,1),node));
        tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        ax=nexttile;
        g=graph(mpc.branch(:,1),mpc.branch(:,2));
        ph=plot(g,'Layout','layered','Sources',1,'NodeCData',V(:,t),...
            'MarkerSize',6,'NodeLabel',string(mpc.bus(:,1)),...
            'EdgeColor',[0.72 0.74 0.77],'LineWidth',0.8);
        highlight(ph,node,'MarkerSize',10);
        colormap(ax,st.sequential);clim([min([V(:);0.95]),max([V(:);1.00])]);
        cb=colorbar;cb.Label.String='节点电压 / p.u.';cb.Label.FontName=st.font_cjk;
        title(sprintf('第 %d 时段',t),'FontName',st.font_cjk);axis off;
        nexttile;
        h=plot(1:size(V,2),V(node,:));line_style(h,st,'single',1);hold on;
        h=plot(t,V(node,t));line_style(h,st,'current',1);
        xlabel('\fontname{Microsoft YaHei}时间 \fontname{Times New Roman}{\itt}/h');
        ylabel('\fontname{Microsoft YaHei}节点电压幅值 \fontname{Times New Roman}{\itV}/p.u.');
        title(sprintf('节点 %d',mpc.bus(node,1)),'FontName',st.font_cjk);
        xlim([1 size(V,2)]);voltage_limits(V(node,:));
    case {'arbitrage','dispatch_storage'}
        price=varargin{1};pch=varargin{2};pdis=varargin{3};soc=varargin{4};
        t=max(1,min(numel(price),round(varargin{5})));dt=varargin{6};erated=varargin{7};
        fig.Position(4)=730;
        tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
        nexttile;
        h=stairs(1:numel(price),price);line_style(h,st,'price',1);hold on;xline(t,':','Color',[0.65 0.65 0.65]);
        ylabel('\fontname{Microsoft YaHei}电价 \fontname{Times New Roman}{\itc}/(\fontname{Microsoft YaHei}元\fontname{Times New Roman}/kWh)');
        xlim([0.5 numel(price)+0.5]);xticklabels({});
        title('分时电价','FontName',st.font_cjk);
        nexttile;
        if strcmp(kind,'dispatch_storage')
            h1=plot(1:numel(price),pch(:));line_style(h1,st,'charge',1);hold on;
            h2=plot(1:numel(price),pdis(:));line_style(h2,st,'discharge',1);
            legend([h1 h2],{'充电功率','放电功率'},'Location','northoutside','Orientation','horizontal');
        else
            pnet=pdis(:)-pch(:);
            h=bar(1:numel(price),pnet,0.64,'FaceColor','flat','EdgeColor',[0.25 0.28 0.32],'LineWidth',0.45);
            h.CData=repmat(st.roles.discharge.color,numel(price),1);
            h.CData(pnet<0,:)=repmat(st.roles.charge.color,sum(pnet<0),1);
            h.Tag='ess_net_power';hold on;
        end
        yline(0,'-','Color',[0.55 0.55 0.55],'LineWidth',0.65);
        xline(t,':','Color',[0.65 0.65 0.65]);
        ylabel('\fontname{Microsoft YaHei}功率 \fontname{Times New Roman}{\itP}/MW');
        xlim([0.5 numel(price)+0.5]);xticklabels({});
        if strcmp(kind,'arbitrage'),title('净放电功率：放电为正，充电为负','FontName',st.font_cjk);end
        nexttile;
        h=plot(0:numel(price),[0.5;soc(:)]*100);line_style(h,st,'soc',1);hold on;
        h=plot(t,soc(t)*100);line_style(h,st,'current',1);
        xlabel('\fontname{Microsoft YaHei}时间 \fontname{Times New Roman}{\itt}/h');
        ylabel('\fontname{Microsoft YaHei}荷电状态 \fontname{Times New Roman}SoC/%');
        title(sprintf('第 %d 时段：充电 %.3f MW，放电 %.3f MW，储存能量 %.3f MWh',...
            t,pch(t),pdis(t),soc(t)*erated),'FontName',st.font_cjk);
        xlim([0 numel(price)*dt]);ylim([0 100]);
    case 'gap'
        gap=varargin{1};t=max(1,min(size(gap,2),round(varargin{2})));
        tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        ax=nexttile;imagesc(1:size(gap,2),1:size(gap,1),gap);colormap(ax,st.sequential);colorbar;
        xlabel('\fontname{Microsoft YaHei}时间 \fontname{Times New Roman}{\itt}/h');
        ylabel('支路编号','FontName',st.font_cjk);title('二阶锥约束松弛间隙','FontName',st.font_cjk);
        nexttile;bar(gap(:,t),0.65,'FaceColor',st.roles.single.color,'EdgeColor','none');
        xlabel('支路编号','FontName',st.font_cjk);ylabel('|UL-P^2-Q^2|','Interpreter','tex');
        title(sprintf('第 %d 时段',t),'FontName',st.font_cjk);
    case 'bigm'
        soc=varargin{1};w=round(varargin{2});lo=0.2;hi=0.8;
        sx=linspace(lo,hi,101);
        dl=max(lo*w,sx-hi*(1-w));du=min(hi*w,sx-lo*(1-w));
        h1=plot(sx,dl);line_style(h1,st,'lower',1);hold on;
        h2=plot(sx,du);line_style(h2,st,'upper',1);
        h3=plot(soc,w*soc);line_style(h3,st,'current',1);
        xlabel('\fontname{Microsoft YaHei}荷电状态 \fontname{Times New Roman}SoC');
        ylabel('\fontname{Microsoft YaHei}辅助变量 \fontname{Times New Roman}\delta');
        title(sprintf('w=%d，SoC=%.2f，辅助变量值为 %.2f',w,soc,w*soc),'FontName',st.font_cjk);
        legend([h1 h2 h3],{'四条约束给出的下界','四条约束给出的上界','当前取值'},...
            'Location','northoutside','Orientation','horizontal');ylim([-0.05 0.9]);
    case 'relaxation_voltage'
        vm=varargin{1};va=varargin{2};node=round(varargin{3});
        h1=plot(1:size(vm,2),vm(node,:));line_style(h1,st,'linear',1);hold on;
        h2=plot(1:size(va,2),va(node,:));line_style(h2,st,'ac',2);
        legend([h1 h2],{'二阶锥模型电压','交流潮流电压'},'Location','northoutside','Orientation','horizontal');
        xlabel('\fontname{Microsoft YaHei}时间 \fontname{Times New Roman}{\itt}/h');
        ylabel('\fontname{Microsoft YaHei}节点电压幅值 \fontname{Times New Roman}{\itV}/p.u.');
        title(sprintf('节点 %d',node),'FontName',st.font_cjk);
        xlim([1 size(vm,2)]);voltage_limits([vm(node,:) va(node,:)]);
    case 'planning'
        vb=varargin{1};vl=varargin{2};ve=varargin{3};node=max(1,min(size(vb,1),round(varargin{4})));
        h1=plot(1:size(vb,2),vb(node,:));line_style(h1,st,'baseline',1);hold on;
        h2=plot(1:size(vb,2),vl(node,:));line_style(h2,st,'linear',1);
        h3=plot(1:size(vb,2),ve(node,:));line_style(h3,st,'ac',2);
        legend([h1 h2 h3],{st.roles.baseline.label,st.roles.linear.label,st.roles.ac.label},...
            'Location','northoutside','NumColumns',3);
        xlabel('\fontname{Microsoft YaHei}时间 \fontname{Times New Roman}{\itt}/h');
        ylabel('\fontname{Microsoft YaHei}节点电压幅值 \fontname{Times New Roman}{\itV}/p.u.');
        title(sprintf('节点 %d',node),'FontName',st.font_cjk);
        xlim([1 size(vb,2)]);voltage_limits([vb(node,:) vl(node,:) ve(node,:)]);
    case 'matrix'
        mpc=varargin{1};b=max(1,min(size(mpc.branch,1),round(varargin{2})));
        nb=size(mpc.bus,1);nl=size(mpc.branch,1);
        af=sparse(mpc.branch(:,1),1:nl,1,nb,nl);at=sparse(mpc.branch(:,2),1:nl,1,nb,nl);
        tiledlayout(1,2,'TileSpacing','compact','Padding','compact');ax=nexttile;
        ph=plot(graph(mpc.branch(:,1),mpc.branch(:,2)),'Layout','layered','Sources',1,...
            'NodeColor',st.roles.single.color,'EdgeColor',[0.72 0.74 0.77]);
        highlight(ph,mpc.branch(b,1),mpc.branch(b,2),'EdgeColor',st.roles.current.color,'LineWidth',2.8);
        title(sprintf('支路 %d：节点 %d 到节点 %d',b,mpc.branch(b,1),mpc.branch(b,2)),'FontName',st.font_cjk);axis off;
        ax=nexttile;imagesc(full(at-af));colormap(ax,st.signed);clim([-1 1]);cb=colorbar;cb.Ticks=[-1 0 1];hold on;
        xline(b,'Color',st.roles.current.color,'LineWidth',1.2);
        xlabel('支路编号','FontName',st.font_cjk);ylabel('节点编号','FontName',st.font_cjk);
        title('末端关联矩阵减首端关联矩阵','FontName',st.font_cjk);
end
axes_handles=findall(fig,'Type','axes');
for k=1:numel(axes_handles)
    ax=axes_handles(k);
    set(ax,'Color','w','XColor','k','YColor','k','ZColor','k','FontName',st.font_latin,...
        'FontSize',st.defaults.font_pt,'LineWidth',st.defaults.axis_pt,'Box','on',...
        'XGrid','off','YGrid','on','GridColor',[0.65 0.68 0.72],'GridAlpha',0.13,'TickDir','out');
    ax.XLabel.FontSize=st.defaults.label_pt;ax.YLabel.FontSize=st.defaults.label_pt;
    ax.Title.FontSize=st.defaults.title_pt;ax.Title.FontWeight='normal';
end
lg=findall(fig,'Type','legend');
set(lg,'FontName',st.font_latin,'Interpreter','tex','FontSize',11.5,'Box','off','Color','w','TextColor',[0.15 0.15 0.15]);
for k=1:numel(lg)
    lg(k).String=mixed_font(lg(k).String,st);
end
text_handles=findall(fig,'Type','text');
for k=1:numel(text_handles)
    h=text_handles(k);
    if strcmp(h.Interpreter,'tex')
        h.String=mixed_font(h.String,st);h.FontName=st.font_latin;
    end
end
set(findall(fig,'Type','colorbar'),'FontName',st.font_latin);
% GraphPlot 的节点与支路标签分别设置字体。
set(findall(fig,'Type','graphplot'),'NodeFontName',st.font_latin,...
    'EdgeFontName',st.font_latin);
end

function out = mixed_font(in,st)
% 已含 TeX 字体或数学命令的标签保留原格式；普通混排分别设置字体。
if iscell(in)
    out=cellfun(@(s)mixed_font(s,st),in,'UniformOutput',false);return;
end
out=char(string(in));
if isempty(out) || contains(out,'\'),return;end
parts=regexp(out,'[^\x00-\x7F]+|[\x00-\x7F]+','match');out='';
for j=1:numel(parts)
    if any(double(parts{j})>127),font=st.font_cjk;else,font=st.font_latin;end
    out=[out '\fontname{' font '}' parts{j}]; %#ok<AGROW>
end
end

function line_style(h,st,role,offset)
r=st.roles.(role);
set(h,'Color',r.color,'LineStyle',r.line,'Marker',r.marker,'LineWidth',st.defaults.line_pt,...
    'MarkerSize',st.defaults.marker_pt,'MarkerFaceColor','w','Tag',['ess_' role]);
if ~strcmp(r.marker,'none') && numel(h.XData)>1
    count=numel(h.XData);stride=max(1,ceil(count/st.defaults.max_markers));
    h.MarkerIndices=min(offset,count):stride:count;
end
if strcmp(role,'current')
    set(h,'MarkerFaceColor',r.color,'MarkerSize',8,'LineWidth',1.1);
end
end

function voltage_limits(values)
values=values(isfinite(values));
lo=min(values(:));hi=max(values(:));
pad=max(0.001,0.08*(hi-lo));ylim([lo-pad hi+pad]);
end
