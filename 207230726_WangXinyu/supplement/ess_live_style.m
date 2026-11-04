function st = ess_live_style
% 课程图形规范的公共入口。颜色与图形角色见项目 04 号文档。
registry=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'figure_style_registry.json')));
st=registry;
keys=fieldnames(registry.roles);
for k=1:numel(keys)
    role=keys{k};st.roles.(role).color=hex_rgb(registry.roles.(role).hex);
end
st.sequential=hex_rgb(registry.colormaps.sequential.colours);
st.signed=hex_rgb(registry.colormaps.signed.colours);
end

function rgb = hex_rgb(hex)
if ischar(hex),hex=cellstr(hex);end
hex=cellstr(string(hex));rgb=zeros(numel(hex),3);
for k=1:numel(hex)
    h=char(hex{k});h=erase(h,'#');
    rgb(k,:)=[hex2dec(h(1:2)),hex2dec(h(3:4)),hex2dec(h(5:6))]/255;
end
end
