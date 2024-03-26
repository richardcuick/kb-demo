clear all;clc;close all

% 计算目标温度、压力下的钻井液密度
% 输入参数：目标温度T、目标压力P、油水比OWR、盐度SaltContent、油相类型OilType
% 输出参数：钻井液密度ρl

T=100;      % 目标温度，℃
P=30*10^6;  % 目标压力，Pa

fc=0.03;    % 化学添加剂体积分数
rho_c=1030; % 化学添加剂密度，kg/m^3
fs=0.03;    % 固相体积分数
rho_s=2500; % 固相密度，kg/m^3
fl=1-fc-fs; % 液相体积分数=1-化学添加剂体积分数-固相体积分数
OWR=0/100;  % OWR=0/100表示纯水基，OWR=20/100表示油水比为20和80，OWR=100/0表示纯油基
fo=OWR*fl;  % 油相绝对体积分数
fw=(1-OWR)*fl;  % 水相绝对体积分数

T_0=20;     % 参考温度，℃
P_0=101325; % 参考压力，Pa

SaltContent=10;    % 盐度（重量，wt%），SaltContent=0表示纯水
rho_b_0=Density_Brine(T_0,P_0,SaltContent); % 求解参考温度和压力下清水及任意盐度下盐水的密度，kg/m^3

OilType=1;  % 取值为1,2,3；OilType=1表示柴油、OilType=2表示矿物油、OilType=3表示合成有机物
rho_o_0=Density_Oil(T_0,P_0,OilType);   % 求解参考温度和压力下油相密度，kg/m^3

rho_b=Density_Brine(T,P,SaltContent);   % 目标温度、压力下的盐水密度，kg/m^3
rho_o=Density_Oil(T,P,OilType);         % 目标温度、压力下的油相密度，kg/m^3

rho_l=(rho_o_0*fo+rho_b_0*fw+rho_s*fs+rho_c*fc)/(1+fo*(rho_o_0/rho_o)+fw*(rho_b_0/rho_b-1));    % 目标温度、压力下的钻井液密度，kg/m^3

