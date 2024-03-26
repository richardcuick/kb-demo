%% ���ڹܹرã������
close all;clear all;clc;
%% �������ݼ�Ԥ�������ʵ�λ�ƣ�
data_rheology=csvread('RheologyData.csv',1,0); % data_rheology�洢�����¶Ⱥ�ѹ�������µ�ţ����������ճ�ȼ����ݣ���һ��Ϊת�٣��ڶ���Ϊ����������ڱ�ͷ��csvֻ�ܶ�ȡ�����ݣ�1����ӵ�2�п�ʼ��ȡ��0����ӵ�1�п�ʼ��ȡ
[mu_0]=Rheology_Newtonian(data_rheology); % �꾮Һճ�ȣ�Pa*s
mu_0=30*10^(-3); % �꾮Һճ�ȣ�Pa*s

welldepth=8000; % ���m
nx=100; % �ռ�������
Nx=nx+1; % �ռ�ڵ���
Depth(1)=0; % �ռ�ڵ㣬m
for i=1:Nx-1
    dx(i)=welldepth/(Nx-1); % �ռ䲽����m
    Depth(i+1)=Depth(i)+dx(i); % �ռ�ڵ㣬m
end

for i=1:1:Nx
    theta(i)=0; % ��б�ǣ���
    theta(i)=theta(i)/180*pi; % ����б�ǵ�λת��Ϊ����
    D_w(i)=0.14; % ����ֱ����m
    D_d_o(i)=0.1143; % �����⾶��m
    D_d_i(i)=0.10053; % �����ھ���m
    A_d(i)=1/4*pi*D_d_i(i)^2; % �����������m^2
    A_a(i)=1/4*pi*(D_w(i)^2-D_d_o(i)^2); % ���ս������m^2
end

T_0=20; % �꾮Һ�����¶ȣ���
P_0=1*10^5; % �꾮Һ����ѹ����Pa
rho_0=1500; % �����¶ȡ�ѹ���µ��꾮Һ�ܶȣ�kg/m^3
g=9.81; % �������ٶȣ�m/s^2

Nx_bot=101; % �����������ڽڵ㣨ȡֵ��Χ1��Nx��
L_bot=Depth(Nx_bot); % ����������ȣ�m
D_bot=D_d_o(Nx); % �����׶��⾶��m
L_p=9.144; % ÿ���������ȣ�m
t_p=100; % ÿ���������ʱ�䣬s
V_p=L_p/t_p; % ��������ٶȣ�m/s

Qv=V_p*(0.25*pi*(D_d_o(Nx)^2-D_d_i(Nx)^2)); % ������������������m^3/s

%% ��Ͳ�����¶�
T_g=2; % �����ݶȣ���/100m
T_s=20; % ע�������¶ȣ���

for x=1:1:Nx
    T_d(x)=T_s+Depth(x)*T_g/100; % ���������¶ȣ���
    T_a(x)=T_s+Depth(x)*T_g/100; % ���������¶ȣ���
end

%%
Qv_a_ass=0; % ������������ֵ��m^3/s

err_BPD=1; % ���վ���ѹ������ھ���ѹ��֮���������
COUNT_Q=1; % ��������������ֵ
a(1)=Qv_a_ass;
b(1)=err_BPD;
while abs(err_BPD)>1*10^(-4) && COUNT_Q<1*10^3
    COUNT_Q=COUNT_Q+1;  % ����ѹ����������
    
    %% ����ѹ������
    % 1�����ڣ�Nx������
    P_a(1)=0*10^6; % ����ѹ����Pa
    rho_a(1)=Density_TP(rho_0,T_0,P_0,T_a(1),P_a(1)); % �����꾮Һ�ܶȣ�kg/m^3
    V_a(1)=Qv_a_ass/A_a(1); % �����꾮Һ���٣�m/s
    mu_a(1)=Rheology_TP(mu_0,T_0,P_0,T_a(1),P_a(1)); % �����꾮Һճ�ȣ�Pa*s
    [PressureSurge_a(1),flow_pattern_a(1)]=PressureSurgeAnnulus(rho_a(1),V_a(1),mu_a(1),D_w(1),D_d_o(1)); % �������嵥λ���Ȳ���ѹ����Pa/m�����������ͣ�������������
    PressureSurge_a(1)=-PressureSurge_a(1); % �������嵥λ���Ȳ���ѹ����Pa/m
    
    for i=2:1:Nx_bot
        P_a_ass(i)=P_a(i-1)+rho_a(i-1)*g*cos(theta(i-1))*dx(i-1); % ����ѹ������ֵ��Pa
        err_P=1; % �������ֵ
        COUNT_P=0; % ����������ֵ
        while err_P>1*10^(-4) && COUNT_P<1*10^3
            COUNT_P=COUNT_P+1; % ��������
            
            rho_a(i)=Density_TP(rho_0,T_0,P_0,T_a(i),P_a_ass(i)); % �����꾮Һ�ܶȣ�kg/m^3
            V_a(i)=Qv_a_ass/A_a(i); % �����꾮Һ���٣�m/s
            mu_a(i)=Rheology_TP(mu_0,T_0,P_0,T_a(i),P_a_ass(i)); % �����꾮Һճ�ȣ�Pa*s
            [PressureSurge_a(i),flow_pattern_a(i)]=PressureSurgeAnnulus(rho_a(i),V_a(i),mu_a(i),D_w(i),D_d_o(i)); % �������嵥λ���Ȳ���ѹ����Pa/m�����������ͣ�������������
            PressureSurge_a(i)=-PressureSurge_a(i); % �������嵥λ���Ȳ���ѹ����Pa/m
            P_a(i)=P_a(i-1)-(rho_a(i)*V_a(i)^2-rho_a(i-1)*V_a(i-1)^2)+dx(i-1)*(rho_a(i-1)*g*cos(theta(i-1))+rho_a(i)*g*cos(theta(i)))/2+dx(i-1)*(PressureSurge_a(i-1)+PressureSurge_a(i))/2; % ����ѹ������ֵ��Pa
            err_P=abs((P_a(i)-P_a_ass(i))/P_a_ass(i)); % ����ѹ��������
            P_a_ass(i)=P_a(i); % �µĻ���ѹ������ֵ��Pa
        end
    end
    
    for i=Nx_bot+1:1:Nx
        P_a_ass(i)=P_a(i-1)+rho_a(i-1)*g*cos(theta(i-1))*dx(i-1); % ����ѹ������ֵ��Pa
        err_P=1; % �������ֵ
        COUNT_P=0; % ����������ֵ
        while err_P>1*10^(-4) && COUNT_P<1*10^3
            COUNT_P=COUNT_P+1; % ��������
            
            rho_a(i)=Density_TP(rho_0,T_0,P_0,T_a(i),P_a_ass(i)); % �����꾮Һ�ܶȣ�kg/m^3
            V_a(i)=0; % �����꾮Һ���٣�m/s
            mu_a(i)=Rheology_TP(mu_0,T_0,P_0,T_a(i),P_a_ass(i)); % �����꾮Һճ�ȣ�Pa*s
            [PressureSurge_a(i),flow_pattern_a(i)]=PressureSurgeAnnulus(rho_a(i),V_a(i),mu_a(i),D_w(i),D_d_o(i)); % �������嵥λ���Ȳ���ѹ����Pa/m�����������ͣ�������������
            P_a(i)=P_a(i-1)-(rho_a(i)*V_a(i)^2-rho_a(i-1)*V_a(i-1)^2)+dx(i-1)*(rho_a(i-1)*g*cos(theta(i-1))+rho_a(i)*g*cos(theta(i)))/2+dx(i-1)*(PressureSurge_a(i-1)+PressureSurge_a(i))/2; % ����ѹ������ֵ��Pa
            err_P=abs((P_a(i)-P_a_ass(i))/P_a_ass(i)); % ����ѹ��������
            P_a_ass(i)=P_a(i); % �µĻ���ѹ������ֵ��Pa
        end
    end
    
    PressureSurge_A(1)=0; % ���ղ���ѹ����Pa
    ECD_PressureSurge_a(1,1)=rho_a(1); % ����ECD��kg/m^3
    for i=2:1:Nx
        PressureSurge_A(i)=PressureSurge_A(i-1)+PressureSurge_a(i)*dx(i-1); % ���ղ���ѹ����Pa
        ECD_PressureSurge_a(i,1)=rho_a(i)+PressureSurge_A(i)/(g*Depth(i)); % ����ECD��kg/m^3
    end
    
    %%
    Qv_d=Qv-Qv_a_ass; % ����������m^3/s
    
    %% ����ѹ������
    % 1�����ڣ�Nx������
    P_d(1)=0*10^6; % ����ѹ����Pa
    rho_d(1)=Density_TP(rho_0,T_0,P_0,T_d(1),P_d(1)); % �����꾮Һ�ܶȣ�kg/m^3
    V_d(1)=Qv_d/A_d(1); % �����꾮Һ���٣�m/s
    mu_d(1)=Rheology_TP(mu_0,T_0,P_0,T_d(1),P_d(1)); % �����꾮Һճ�ȣ�Pa*s
    [PressureSurge_d(1),flow_pattern_d(1)]=PressureSurgeDrillpipe(rho_d(1),V_d(1),mu_d(1),D_d_i(1)); % �������嵥λ���Ȳ���ѹ����Pa/m�����������ͣ�������������
    PressureSurge_d(1)=-PressureSurge_d(1); % �������嵥λ���Ȳ���ѹ����Pa/m
    
    for i=2:1:Nx_bot
        P_d_ass(i)=P_d(i-1)+rho_d(i-1)*g*cos(theta(i-1))*dx(i-1); % ����ѹ������ֵ��Pa
        err_P=1; % �������ֵ
        COUNT_P=0; % ����������ֵ
        while err_P>1*10^(-4) && COUNT_P<1*10^3
            COUNT_P=COUNT_P+1; % ��������
            
            rho_d(i)=Density_TP(rho_0,T_0,P_0,T_d(i),P_d_ass(i)); % �����꾮Һ�ܶȣ�kg/m^3
            V_d(i)=Qv_d/A_d(i); % �����꾮Һ���٣�m/s
            mu_d(i)=Rheology_TP(mu_0,T_0,P_0,T_d(i),P_d_ass(i)); % �����꾮Һճ�ȣ�Pa*s
            [PressureSurge_d(i),flow_pattern_d(i)]=PressureSurgeDrillpipe(rho_d(i),V_d(i),mu_d(i),D_d_i(i)); % �������嵥λ���Ȳ���ѹ����Pa/m�����������ͣ�������������
            PressureSurge_d(i)=-PressureSurge_d(i); % �������嵥λ���Ȳ���ѹ����Pa/m
            P_d(i)=P_d(i-1)+((rho_d(i-1)*V_d(i-1)^2)-(rho_d(i)*V_d(i)^2))+((rho_d(i-1)*g*cos(theta(i-1))+PressureSurge_d(i-1))+(rho_d(i)*g*cos(theta(i))+PressureSurge_d(i)))*(dx(i-1)/2); % ����ѹ������ֵ��Pa
            err_P=abs((P_d(i)-P_d_ass(i))/P_d_ass(i)); % ����ѹ��������
            P_d_ass(i)=P_d(i); % �µĹ���ѹ������ֵ��Pa
        end
    end
    
    for i=Nx_bot+1:1:Nx
        P_d_ass(i)=P_d(i-1)+rho_d(i-1)*g*cos(theta(i-1))*dx(i-1); % ����ѹ������ֵ��Pa
        err_P=1; % �������ֵ
        COUNT_P=0; % ����������ֵ
        while err_P>1*10^(-4) && COUNT_P<1*10^3
            COUNT_P=COUNT_P+1; % ��������
            
            rho_d(i)=Density_TP(rho_0,T_0,P_0,T_d(i),P_d_ass(i)); % �����꾮Һ�ܶȣ�kg/m^3
            V_d(i)=0; % �����꾮Һ���٣�m/s
            mu_d(i)=Rheology_TP(mu_0,T_0,P_0,T_d(i),P_d_ass(i)); % �����꾮Һճ�ȣ�Pa*s
            [PressureSurge_d(i),flow_pattern_d(i)]=PressureSurgeDrillpipe(rho_d(i),V_d(i),mu_d(i),D_d_i(i)); % �������嵥λ���Ȳ���ѹ����Pa/m�����������ͣ�������������
            P_d(i)=P_d(i-1)+((rho_d(i-1)*V_d(i-1)^2)-(rho_d(i)*V_d(i)^2))+((rho_d(i-1)*g*cos(theta(i-1))+PressureSurge_d(i-1))+(rho_d(i)*g*cos(theta(i))+PressureSurge_d(i)))*(dx(i-1)/2); % ����ѹ������ֵ��Pa
            err_P=abs((P_d(i)-P_d_ass(i))/P_d_ass(i)); % ����ѹ��������
            P_d_ass(i)=P_d(i); % �µĹ���ѹ������ֵ��Pa
        end
    end
    
    PressureSurge_D(1)=0; % ���ڲ���ѹ����Pa
    ECD_PressureSurge_d(1,1)=rho_a(1); % ����ECD��kg/m^3
    for i=2:1:Nx
        PressureSurge_D(i)=PressureSurge_D(i-1)+PressureSurge_d(i)*dx(i-1); % ���ڲ���ѹ����Pa
        ECD_PressureSurge_d(i,1)=rho_a(i)+PressureSurge_D(i)/(g*Depth(i)); % ����ECD��kg/m^3
    end
    
    err_BPD=abs(P_a(Nx)-P_d(Nx))/P_a(Nx); % ���վ���ѹ������ھ���ѹ��֮���������
    
    if (P_a(Nx)-P_d(Nx))>0 % ���ݻ��վ���ѹ������ھ���ѹ��֮��Ĵ�С��ϵ���Ի�����������ֵ���е��ڣ������󻷿�������
        Qv_a_ass=Qv_a_ass+Qv_d/2; % �µĻ�����������ֵ��m^3/s
    else % �����С����������
        Qv_a_ass=Qv_a_ass-abs(a(COUNT_Q-1)-a(COUNT_Q-2))/2; % �µĻ�����������ֵ��m^3/s
    end
    a(COUNT_Q)=Qv_a_ass;
    b(COUNT_Q)=err_BPD;
end

%% ��������
fprintf('ţ��ģ�ͣ����ڹܹر������\n');
fprintf('ÿ���������ȣ�m����%8.5f\n',L_p); % ÿ���������ȣ�m
fprintf('ÿ���������ʱ�䣨s����%8.5f\n',t_p); % ÿ���������ʱ�䣬s
fprintf('��������ٶȣ�m/s����%8.5f\n',abs(V_p)); % ��������ٶȣ�m/s
fprintf('����ѹ����MPa����%8.5f\n',PressureSurge_A(Nx)/10^6); % �������ѹ����MPa
fprintf('����ѹ����MPa����%8.5f\n',P_a(Nx)/10^6); % �������ѹ����MPa
fprintf('����ECD��kg/m^3����%8.5f\n',ECD_PressureSurge_a(Nx,1)); % �������ECD��kg/m^3
fprintf('�ܵ�ECD��kg/m^3����%8.5f\n',ECD_PressureSurge_a(Nx_bot,1)); % ����ܵ�ECD��kg/m^3

%% ��ͼ
% figure(1); % ����ѹ��vs����
% plot(P_a/10^6,Depth); % ����ѹ����MPa
% xlabel('ѹ����MPa��','FontName','����');
% ylabel('���m��','FontName','����');
% set(gca,'YDir','reverse');
% set(gca,'fontsize',16);
% box on; % ��ʾ������ı߿�
% grid on; % ��ʾ���������������
% grid minor; % ��ʾ������Ĵ�������
% 
% figure(2); % ����ECDvs����
% plot(ECD_PressureSurge(2:Nx),Depth(2:Nx)); % ����ECD��kg/m^3
% xlabel('����ECD��kg/m^3��','FontName','����');
% ylabel('���m��','FontName','����');
% set(gca,'YDir','reverse');
% set(gca,'fontsize',16);
% box on; % ��ʾ������ı߿�
% grid on; % ��ʾ���������������
% grid minor; % ��ʾ������Ĵ�������