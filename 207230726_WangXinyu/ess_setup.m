function ess_setup(folder)
% 仅修改本次 MATLAB 会话路径，不调用 savepath，不改系统配置。
if nargin==0,folder=fileparts(mfilename('fullpath'));end
cd(folder); addpath(folder); addpath(fullfile(folder,'cases'));
bundle=fileparts(folder);
% MOSEK 独立解压安装，使用不覆盖 MATLAB 优化函数的 R2019bom 接口。
mc={getenv('ESS_MOSEK_ROOT'),fullfile(bundle,'deps','mosek','11.2'),...
    fullfile(getenv('USERPROFILE'),'Documents','Codex','Tools','MOSEK','11.2')};
for k=1:numel(mc)
    if ~isempty(mc{k}) && isfile(fullfile(mc{k},'toolbox','R2019bom','mosekopt.mexw64'))
        mbin=fullfile(mc{k},'tools','platform','win64x86','bin');
        setenv('PATH',[mbin pathsep getenv('PATH')]);
        addpath(fullfile(mc{k},'toolbox','R2019bom'),'-begin');break;
    end
end
if isempty(getenv('MOSEKLM_LICENSE_FILE'))
    lic=fullfile(getenv('USERPROFILE'),'Documents','Codex','Tools','MOSEK','licenses','mosek.lic');
    if isfile(lic),setenv('MOSEKLM_LICENSE_FILE',lic);end
end
pydeps=fullfile(bundle,'deps','python');
if exist(pydeps,'dir'),setenv('PYTHONPATH',[pydeps pathsep getenv('PYTHONPATH')]);end
if exist(fullfile(bundle,'deps','yalmip'),'dir'),addpath(genpath(fullfile(bundle,'deps','yalmip')));end
mp=fullfile(bundle,'deps','matpower');
if exist(mp,'dir')
    addpath(fullfile(mp,'lib'));addpath(fullfile(mp,'data'));
    addpath(fullfile(mp,'mips','lib'));addpath(fullfile(mp,'mp-opt-model','lib'));
    addpath(fullfile(mp,'mptest','lib'));addpath(fullfile(mp,'most','lib'));
end
addpath(fullfile(matlabroot,'toolbox','optim','optim'),'-begin');
assert(exist('runpf','file')~=0,'请配置 MATPOWER 8.1');
assert(exist('sdpvar','file')~=0,'请配置 YALMIP');
end
