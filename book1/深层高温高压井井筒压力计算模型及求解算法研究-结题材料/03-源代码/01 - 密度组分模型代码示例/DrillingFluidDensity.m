clear all;clc;close all

% ����Ŀ���¶ȡ�ѹ���µ��꾮Һ�ܶ�
% ���������Ŀ���¶�T��Ŀ��ѹ��P����ˮ��OWR���ζ�SaltContent����������OilType
% ����������꾮Һ�ܶȦ�l

T=100;      % Ŀ���¶ȣ���
P=30*10^6;  % Ŀ��ѹ����Pa

fc=0.03;    % ��ѧ��Ӽ��������
rho_c=1030; % ��ѧ��Ӽ��ܶȣ�kg/m^3
fs=0.03;    % �����������
rho_s=2500; % �����ܶȣ�kg/m^3
fl=1-fc-fs; % Һ���������=1-��ѧ��Ӽ��������-�����������
OWR=0/100;  % OWR=0/100��ʾ��ˮ����OWR=20/100��ʾ��ˮ��Ϊ20��80��OWR=100/0��ʾ���ͻ�
fo=OWR*fl;  % ��������������
fw=(1-OWR)*fl;  % ˮ������������

T_0=20;     % �ο��¶ȣ���
P_0=101325; % �ο�ѹ����Pa

SaltContent=10;    % �ζȣ�������wt%����SaltContent=0��ʾ��ˮ
rho_b_0=Density_Brine(T_0,P_0,SaltContent); % ���ο��¶Ⱥ�ѹ������ˮ�������ζ�����ˮ���ܶȣ�kg/m^3

OilType=1;  % ȡֵΪ1,2,3��OilType=1��ʾ���͡�OilType=2��ʾ�����͡�OilType=3��ʾ�ϳ��л���
rho_o_0=Density_Oil(T_0,P_0,OilType);   % ���ο��¶Ⱥ�ѹ���������ܶȣ�kg/m^3

rho_b=Density_Brine(T,P,SaltContent);   % Ŀ���¶ȡ�ѹ���µ���ˮ�ܶȣ�kg/m^3
rho_o=Density_Oil(T,P,OilType);         % Ŀ���¶ȡ�ѹ���µ������ܶȣ�kg/m^3

rho_l=(rho_o_0*fo+rho_b_0*fw+rho_s*fs+rho_c*fc)/(1+fo*(rho_o_0/rho_o)+fw*(rho_b_0/rho_b-1));    % Ŀ���¶ȡ�ѹ���µ��꾮Һ�ܶȣ�kg/m^3

