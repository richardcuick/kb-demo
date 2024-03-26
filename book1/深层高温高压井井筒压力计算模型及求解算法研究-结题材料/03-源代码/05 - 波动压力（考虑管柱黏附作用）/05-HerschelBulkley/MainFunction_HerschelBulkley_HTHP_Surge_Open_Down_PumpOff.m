%% 开口管关泵（下入）
close all;clear all;clc;
%% 基础数据及预处理（国际单位制）
data_rheology=csvread('RheologyData.csv',1,0); % data_rheology存储测试温度和压力条件下的赫巴流体六速粘度计数据，第一列为转速，第二列为读数；因存在表头且csv只能读取纯数据，1代表从第2行开始读取，0代表从第1列开始读取
[tau_y,K_0,n_0]=Rheology_HerschelBulkley(data_rheology); % 计算钻井液屈服应力（Pa）、稠度系数（Pa·s^n）和流性指数
tau_y=5; % 屈服应力，Pa
K_0=2; % 稠度系数，Pa·s^n
n_0=0.6; % 流性指数

welldepth=8000; % 井深，m
nx=100; % 空间网格数
Nx=nx+1; % 空间节点数
Depth(1)=0; % 空间节点，m
for i=1:Nx-1
    dx(i)=welldepth/(Nx-1); % 空间步长，m
    Depth(i+1)=Depth(i)+dx(i); % 空间节点，m
end

for i=1:1:Nx
    theta(i)=0; % 井斜角，°
    theta(i)=theta(i)/180*pi; % 将井斜角单位转换为弧度
    D_w(i)=0.14; % 井眼直径，m
    D_d_o(i)=0.1143; % 管柱外径，m
    D_d_i(i)=0.10053; % 管柱内径，m
    A_d(i)=1/4*pi*D_d_i(i)^2; % 管柱截面积，m^2
    A_a(i)=1/4*pi*(D_w(i)^2-D_d_o(i)^2); % 环空截面积，m^2
end

T_0=20; % 钻井液测试温度，℃
P_0=1*10^5; % 钻井液测试压力，Pa
rho_0=1500; % 测试温度、压力下的钻井液密度，kg/m^3
g=9.81; % 重力加速度，m/s^2

Nx_bot=101; % 管柱下深所在节点（取值范围1～Nx）
L_bot=Depth(Nx_bot); % 管柱下入深度，m
D_bot=D_d_o(Nx); % 管柱底端外径，m
L_p=9.144; % 每根管柱长度，m
t_p=100; % 每根管柱下入时间，s
V_p=L_p/t_p; % 管柱下入速度，m/s

Qv=V_p*(0.25*pi*(D_d_o(Nx)^2-D_d_i(Nx)^2))+0.5*V_p*(A_a(Nx)+A_d(Nx)); % 管柱下入引起的流量，m^3/s

%% 井筒流体温度
T_g=2; % 地温梯度，℃/100m
T_s=20; % 注入流体温度，℃

for x=1:1:Nx
    T_d(x)=T_s+Depth(x)*T_g/100; % 管内流体温度，℃
    T_a(x)=T_s+Depth(x)*T_g/100; % 环空流体温度，℃
end

%%
Qv_a_ass=0; % 环空流量假设值，m^3/s

err_BPD=1; % 环空井底压力与管内井底压力之间的相对误差
COUNT_Q=1; % 排量迭代次数初值
a(1)=Qv_a_ass;
b(1)=err_BPD;
while abs(err_BPD)>1*10^(-4) && COUNT_Q<1*10^3
    COUNT_Q=COUNT_Q+1;  % 出口压力迭代次数
    
    %% 环空压力计算
    % 1代表井口，Nx代表井底
    P_a(1)=0.1*10^6; % 井口压力，Pa
    rho_a(1)=Density_TP(rho_0,T_0,P_0,T_a(1),P_a(1)); % 井口钻井液密度，kg/m^3
    V_a(1)=Qv_a_ass/A_a(1); % 井口钻井液流速，m/s
    [K_a(1),n_a(1)]=Rheology_TP(K_0,n_0,T_0,P_0,T_a(1),P_a(1)); % 井口钻井液稠度系数（Pa·s^n）和流性指数
    [PressureSurge_a(1),flow_pattern_a(1)]=PressureSurgeAnnulus(rho_a(1),V_a(1),K_a(1),n_a(1),D_w(1),D_d_o(1)); % 井口流体单位长度波动压力（Pa/m）和流体流型（层流、过渡流、湍流）
    
    for i=2:1:Nx_bot
        P_a_ass(i)=P_a(i-1)+rho_a(i-1)*g*cos(theta(i-1))*dx(i-1); % 环空压力假设值，Pa
        err_P=1; % 相对误差初值
        COUNT_P=0; % 迭代次数初值
        while err_P>1*10^(-4) && COUNT_P<1*10^3
            COUNT_P=COUNT_P+1; % 迭代次数
            
            rho_a(i)=Density_TP(rho_0,T_0,P_0,T_a(i),P_a_ass(i)); % 环空钻井液密度，kg/m^3
            V_a(i)=Qv_a_ass/A_a(i); % 环空钻井液流速，m/s
            [K_a(i),n_a(i)]=Rheology_TP(K_0,n_0,T_0,P_0,T_a(i),P_a_ass(i)); % 环空钻井液稠度系数（Pa·s^n）和流性指数
            [PressureSurge_a(i),flow_pattern_a(i)]=PressureSurgeAnnulus(rho_a(i),V_a(i),K_a(i),n_a(i),D_w(i),D_d_o(i)); % 环空流体单位长度波动压力（Pa/m）和流体流型（层流、过渡流、湍流）
            P_a(i)=P_a(i-1)-(rho_a(i)*V_a(i)^2-rho_a(i-1)*V_a(i-1)^2)+dx(i-1)*(rho_a(i-1)*g*cos(theta(i-1))+rho_a(i)*g*cos(theta(i)))/2+dx(i-1)*(PressureSurge_a(i-1)+PressureSurge_a(i))/2; % 环空压力计算值，Pa
            err_P=abs((P_a(i)-P_a_ass(i))/P_a_ass(i)); % 环空压力相对误差
            P_a_ass(i)=P_a(i); % 新的环空压力假设值，Pa
        end
    end
    
    for i=Nx_bot+1:1:Nx
        P_a_ass(i)=P_a(i-1)+rho_a(i-1)*g*cos(theta(i-1))*dx(i-1); % 环空压力假设值，Pa
        err_P=1; % 相对误差初值
        COUNT_P=0; % 迭代次数初值
        while err_P>1*10^(-4) && COUNT_P<1*10^3
            COUNT_P=COUNT_P+1; % 迭代次数
            
            rho_a(i)=Density_TP(rho_0,T_0,P_0,T_a(i),P_a_ass(i)); % 环空钻井液密度，kg/m^3
            V_a(i)=0; % 环空钻井液流速，m/s
            [K_a(i),n_a(i)]=Rheology_TP(K_0,n_0,T_0,P_0,T_a(i),P_a_ass(i)); % 环空钻井液稠度系数（Pa·s^n）和流性指数
            [PressureSurge_a(i),flow_pattern_a(i)]=PressureSurgeAnnulus(rho_a(i),V_a(i),K_a(i),n_a(i),D_w(i),D_d_o(i)); % 环空流体单位长度波动压力（Pa/m）和流体流型（层流、过渡流、湍流）
            P_a(i)=P_a(i-1)-(rho_a(i)*V_a(i)^2-rho_a(i-1)*V_a(i-1)^2)+dx(i-1)*(rho_a(i-1)*g*cos(theta(i-1))+rho_a(i)*g*cos(theta(i)))/2+dx(i-1)*(PressureSurge_a(i-1)+PressureSurge_a(i))/2; % 环空压力计算值，Pa
            err_P=abs((P_a(i)-P_a_ass(i))/P_a_ass(i)); % 环空压力相对误差
            P_a_ass(i)=P_a(i); % 新的环空压力假设值，Pa
        end
    end
    
    PressureSurge_A(1)=0; % 波动压力，Pa
    ECD_PressureSurge(1,1)=rho_a(1); % 环空ECD，kg/m^3
    for i=2:1:Nx
        PressureSurge_A(i)=PressureSurge_A(i-1)+PressureSurge_a(i)*dx(i-1); % 波动压力，Pa
        ECD_PressureSurge(i,1)=rho_a(i)+PressureSurge_A(i)/(g*Depth(i)); % 环空ECD，kg/m^3
    end
    
    %%
    Qv_d=Qv-Qv_a_ass; % 管内流量，m^3/s
    
    %% 管内压力计算
    % 1代表井口，Nx代表井底
    P_d(1)=0*10^6; % 井口压力，Pa
    rho_d(1)=Density_TP(rho_0,T_0,P_0,T_d(1),P_d(1)); % 井口钻井液密度，kg/m^3
    V_d(1)=Qv_d/A_d(1); % 井口钻井液流速，m/s
    [K_d(1),n_d(1)]=Rheology_TP(K_0,n_0,T_0,P_0,T_d(1),P_d(1)); % 井口钻井液稠度系数（Pa·s^n）和流性指数
    [PressureSurge_d(1),flow_pattern_d(1)]=PressureSurgeDrillPipe(rho_d(1),V_d(1),K_d(1),n_d(1),D_d_i(1)); % 井口流体单位长度波动压力（Pa/m）和流体流型（层流、过渡流、湍流）
    
    for i=2:1:Nx_bot
        P_d_ass(i)=P_d(i-1)+rho_d(i-1)*g*cos(theta(i-1))*dx(i-1); % 管内压力假设值，Pa
        err_P=1; % 相对误差初值
        COUNT_P=0; % 迭代次数初值
        while err_P>1*10^(-4) && COUNT_P<1*10^3
            COUNT_P=COUNT_P+1; % 迭代次数
            
            rho_d(i)=Density_TP(rho_0,T_0,P_0,T_d(i),P_d_ass(i)); % 管内钻井液密度，kg/m^3
            V_d(i)=Qv_d/A_d(i); % 管内钻井液流速，m/s
            [K_d(i),n_d(i)]=Rheology_TP(K_0,n_0,T_0,P_0,T_d(i),P_d_ass(i)); % 环空钻井液稠度系数（Pa·s^n）和流性指数
            [PressureSurge_d(i),flow_pattern_d(i)]=PressureSurgeDrillPipe(rho_d(i),V_d(i),K_d(i),n_d(i),D_d_i(i)); % 环空流体单位长度波动压力（Pa/m）和流体流型（层流、过渡流、湍流）
            P_d(i)=P_d(i-1)+((rho_d(i-1)*V_d(i-1)^2)-(rho_d(i)*V_d(i)^2))+((rho_d(i-1)*g*cos(theta(i-1))+PressureSurge_d(i-1))+(rho_d(i)*g*cos(theta(i))+PressureSurge_d(i)))*(dx(i-1)/2); % 管内压力计算值，Pa
            err_P=abs((P_d(i)-P_d_ass(i))/P_d_ass(i)); % 管内压力相对误差
            P_d_ass(i)=P_d(i); % 新的管内压力假设值，Pa
        end
    end
    
    for i=Nx_bot+1:1:Nx
        P_d_ass(i)=P_d(i-1)+rho_d(i-1)*g*cos(theta(i-1))*dx(i-1); % 管内压力假设值，Pa
        err_P=1; % 相对误差初值
        COUNT_P=0; % 迭代次数初值
        while err_P>1*10^(-4) && COUNT_P<1*10^3
            COUNT_P=COUNT_P+1; % 迭代次数
            
            rho_d(i)=Density_TP(rho_0,T_0,P_0,T_d(i),P_d_ass(i)); % 管内钻井液密度，kg/m^3
            V_d(i)=0; % 管内钻井液流速，m/s
            [K_d(i),n_d(i)]=Rheology_TP(K_0,n_0,T_0,P_0,T_d(i),P_d_ass(i)); % 环空钻井液稠度系数（Pa·s^n）和流性指数
            [PressureSurge_d(i),flow_pattern_d(i)]=PressureSurgeDrillPipe(rho_d(i),V_d(i),K_d(i),n_d(i),D_d_i(i)); % 环空流体单位长度波动压力（Pa/m）和流体流型（层流、过渡流、湍流）
            P_d(i)=P_d(i-1)+((rho_d(i-1)*V_d(i-1)^2)-(rho_d(i)*V_d(i)^2))+((rho_d(i-1)*g*cos(theta(i-1))+PressureSurge_d(i-1))+(rho_d(i)*g*cos(theta(i))+PressureSurge_d(i)))*(dx(i-1)/2); % 管内压力计算值，Pa
            err_P=abs((P_d(i)-P_d_ass(i))/P_d_ass(i)); % 管内压力相对误差
            P_d_ass(i)=P_d(i); % 新的管内压力假设值，Pa
        end
    end
    
    PressureSurge_D(1)=0; % 管内波动压力，Pa
    ECD_PressureSurge_d(1,1)=rho_a(1); % 管内ECD，kg/m^3
    for i=2:1:Nx
        PressureSurge_D(i)=PressureSurge_D(i-1)+PressureSurge_d(i)*dx(i-1); % 管内波动压力，Pa
        ECD_PressureSurge_d(i,1)=rho_a(i)+PressureSurge_D(i)/(g*Depth(i)); % 管内ECD，kg/m^3
    end
    
    err_BPD=abs(P_a(Nx)-P_d(Nx))/P_a(Nx); % 环空井底压力与管内井底压力之间的相对误差
    
    if (P_a(Nx)-P_d(Nx))<0 % 根据环空井底压力与管内井底压力之间的大小关系，对环空流量假设值进行调节（需增大环空流量）
        Qv_a_ass=Qv_a_ass+Qv_d/2; % 新的环空流量假设值，m^3/s
    else % （需减小环空流量）
        Qv_a_ass=Qv_a_ass-abs(a(COUNT_Q-1)-a(COUNT_Q-2))/2; % 新的环空流量假设值，m^3/s
    end
    a(COUNT_Q)=Qv_a_ass;
    b(COUNT_Q)=err_BPD;
end

%% 命令窗口输出
fprintf('赫巴模型（开口管关泵下入）\n');
fprintf('每根管柱长度（m）：%8.5f\n',L_p); % 每根管柱长度，m
fprintf('每根管柱下入时间（s）：%8.5f\n',t_p); % 每根管柱下入时间，s
fprintf('管柱下入速度（m/s）：%8.5f\n',abs(V_p)); % 管柱下入速度，m/s
fprintf('波动压力（MPa）：%8.5f\n',PressureSurge_A(Nx)/10^6); % 输出波动压力，MPa
fprintf('井底压力（MPa）：%8.5f\n',P_a(Nx)/10^6); % 输出井底压力，MPa
fprintf('井底ECD（kg/m^3）：%8.5f\n',ECD_PressureSurge(Nx,1)); % 输出井底ECD，kg/m^3
fprintf('管底ECD（kg/m^3）：%8.5f\n',ECD_PressureSurge(Nx_bot,1)); % 输出管底ECD，kg/m^3

%% 绘图
% figure(1); % 环空压力vs井深
% plot(P_a/10^6,Depth); % 环空压力，MPa
% xlabel('井口压力（MPa）','FontName','黑体');
% ylabel('井深（m）','FontName','黑体');
% set(gca,'YDir','reverse');
% set(gca,'fontsize',16);
% box on; % 显示坐标轴的边框
% grid on; % 显示坐标轴的主网格线
% grid minor; % 显示坐标轴的次网格线
% 
% figure(2); % 环空ECDvs井深
% plot(ECD_PressureSurge(2:Nx),Depth(2:Nx)); % 环空ECD，kg/m^3
% xlabel('环空ECD（kg/m^3）','FontName','黑体');
% ylabel('井深（m）','FontName','黑体');
% set(gca,'YDir','reverse');
% set(gca,'fontsize',16);
% box on; % 显示坐标轴的边框
% grid on; % 显示坐标轴的主网格线
% grid minor; % 显示坐标轴的次网格线