%% 泡沫瞬态冲砂（国际单位制）
close all;clear all;clc;
%% 导入测斜数据，构造原始空间节点DEPTH
DATA=csvread("TrajectoryData.csv",1,0); % 导入测深（m）、井斜角（°）、井斜方位角（°）（导入）
Dm=DATA(:,1); % 测深，m
alpha=DATA(:,2); % 井斜角，°
phi=DATA(:,3); % 井斜方位角，°
[r,c]=size(DATA); % 判断有多少行、列数据
L_s_t=5800; % 砂床顶部深度（等于连续油管下入过程的最终下入深度），m（输入）
L_s_b=6200; % 砂床底部深度（等于冲洗钻进过程的最终下入深度），m（输入）
Lp=4000; % 回拖短起深度，m（输入）

nx_1=100; % 井口至回拖短起深度段（0～Lp）原始空间网格数（输入）
Nx_1=nx_1+1; % 井口至回拖短起深度段（0～Lp）空间节点数
dx_1=Lp/nx_1; % 井口至回拖短起深度段（0～Lp）空间步长，m

for i=1:1:Nx_1
    DEPTH(i)=dx_1*(i-1); % 原始空间节点，m
end

nx_2=50; % 回拖短起深度至砂床顶部深度段（Lp～L_s_t）原始空间网格数（输入）
Nx_2=nx_2+1; % 回拖短起深度至砂床顶部深度段（Lp～L_s_t）空间节点数
dx_2=(L_s_t-Lp)/nx_2; % 回拖短起深度至砂床顶部深度段（Lp～L_s_t）空间步长，m

for i=Nx_1+1:1:Nx_1+Nx_2-1
    DEPTH(i)=DEPTH(i-1)+dx_2; % 原始空间节点，m
end

nx_3=50; % 砂床顶部深度至砂床底部深度段（L_s_t～L_s_b）原始空间网格数（输入）
Nx_3=nx_3+1; % 砂床顶部深度至砂床底部深度段（L_s_t～L_s_b）空间节点数
dx_3=(L_s_b-L_s_t)/nx_3; % 砂床顶部深度至砂床底部深度段（L_s_t～L_s_b）空间步长，m

for i=Nx_1+Nx_2:1:Nx_1+Nx_2+Nx_3-2
    DEPTH(i)=DEPTH(i-1)+dx_3; % 原始空间节点，m
end

Nx=Nx_1+Nx_2+Nx_3-2; % 原始空间节点数

%% 特殊节点
Dsp_wc=[3980]; % 特殊点井深（变径处），m（完井数据产生）

t_pen=2; % 总冲洗钻进次数（输入）
L_pen=(L_s_b-L_s_t)/t_pen; % 每次冲洗钻进长度，m
for i=1:1:t_pen
    Dsp_pen(i)=L_s_t+L_pen*i; % 特殊点井深（变径处），m（冲洗钻进产生）
end

Dsp=[Dsp_wc,Dsp_pen]; % 特殊点井深（变径处），m（输入或由完井数据决定）
n_sp=length(Dsp); % 判断有多少个特殊点

%% 添加特殊点井深节点，构造实际使用空间节点Depth
Depth=DEPTH;
for i=1:1:n_sp
    for j=1:1:Nx-1
        if Dsp(i)>DEPTH(j) && Dsp(i)<DEPTH(j+1) % 理论上此时需要添加空间节点
            Depth(j+1)=Dsp(i); % 将特殊点井深存到Depth，m
            for k=j+1:1:Nx
                Depth(k+1)=DEPTH(k); % 将DEPTH(j+1)之后的空间节点赋给Depth
            end
            Nx=Nx+1; % 空间节点数加一
            DEPTH=Depth; % 更新DEPTH
        end
    end
end
nx=Nx-1; % 空间网格数

Nx_Lp=0; % 回拖短起深度以上的空间节点数
for x=1:1:Nx
    if Depth(x)<=Lp
        Nx_Lp=Nx_Lp+1; % 回拖短起深度以上的空间节点数
    end
end

Nx_L_s_t=0; % 砂床顶部深度以上的空间节点数
for x=1:1:Nx
    if Depth(x)<=L_s_t
        Nx_L_s_t=Nx_L_s_t+1; % 砂床顶部深度以上的空间节点数
    end
end

for i=1:1:length(Dsp_pen)
    Nx_Dsp_pen(i)=0; % 冲洗钻进深度以上的空间节点数
    for x=1:1:Nx
        if Depth(x)<=Dsp_pen(i)
            Nx_Dsp_pen(i)=Nx_Dsp_pen(i)+1; % 冲洗钻进深度以上的空间节点数
        end
    end
end

%% 插值计算空间节点Depth对应的井斜角theta
for x=1:1:Nx
    for i=1:1:r-1
        if Depth(x)>=Dm(i) && Depth(x)<Dm(i+1)
            theta(x)=alpha(i)+(alpha(i+1)-alpha(i))*(Depth(x)-Dm(i))/(Dm(i+1)-Dm(i)); % 线性插值得到各节点井斜角，°
        elseif Depth(x)==Dm(i+1)
            theta(x)=alpha(i+1); % 井斜角，°
        end
    end
end

data=zeros(length(Depth),2); % 计算使用井深及井斜角
data(:,1)=Depth; % 井深，m
data(:,2)=theta; % 井斜角，°



%% 冲洗钻进过程1（国际单位制）
%% 计算空间步长dx及相应时间步长dt
V1=(Dsp_pen(1)-L_s_t)/(1*3600); % 冲洗钻进速度，m/s（输入）
t_1=(Dsp_pen(1)-L_s_t)/V1; % 冲洗钻进总时长，s

Nt=Nx_Dsp_pen(1)-Nx_L_s_t+1; % 时间节点数
nt=Nt-1; % 时间网格数

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % 空间步长，m
end

Time(1)=0; % 冲洗钻进总时长初值，s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Dsp_pen(1)-Nt+t)/V1; % 每冲洗钻进一个空间步长所需时间，s
    Time(t+1)=Time(t)+dt(t); % 冲洗钻进至第(t+1)个空间节点所经历的总时长，s
end

%% 计算不同时刻连续管冲洗钻进深度L_coil及盘管段连续管长度L_reel
L=10000; % 连续油管总长，m（输入）
L_wg=8; % 井口到注入头顶部段连续管长度，m（输入）
L_goose=3; % 导向器段连续管长度，m（输入）
D_goose=2; % 导向器段半径，m（输入）
H_goose=10; % 导向器距地面高度，m（输入）
L_gr=20; % 导向器到滚筒段连续管长度，m（输入）
theta_gr=acosd(H_goose/L_gr); % 导向器到滚筒段连续管与铅垂线夹角，°
D_r_i=3; % 滚筒内径，m（输入）
D_r_o=5; % 滚筒外径，m（输入）
W_r=5; % 滚筒宽度，m（输入）
D_cable=0.005; % 电缆外径，m

L_coil(1)=L_s_t; % 初始时刻连续管底部深度，m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % 初始时刻盘管段连续管长度，m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)+dx(Nx_Dsp_pen(1)-Nt+t-1); % 冲洗钻进深度，m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % 盘管段连续管长度，m
end
L_cable=L_coil; % 电缆下入长度，m

%% 数据输入及预处理
D_ct_o_0=0.04445; % 连续油管外径，m（输入）
L1=2000; % 近出口第一段连续管长度，m（输入）
D_ct_i_1=0.03709; % 近出口第一段连续管内径，m（输入）
L2=2000; % 近出口第二段连续管长度，m（输入）
D_ct_i_2=0.03653; % 近出口第二段连续管内径，m（输入）
L3=2000; % 近出口第三段连续管长度，m（输入）
D_ct_i_3=0.03555; % 近出口第三段连续管内径，m（输入）
L4=2000; % 近出口第四段连续管长度，m（输入）
D_ct_i_4=0.03489; % 近出口第四段连续管内径，m（输入）
L5=2000; % 近出口第五段连续管长度，m（输入）
D_ct_i_5=0.03409; % 近出口第五段连续管内径，m（输入）

L_t_1=4000; % 上部油管（或套管或裸眼）长度，m（输入）
D_t_i_1=0.09718; %0.068;% 上部油管（或套管或裸眼）内径，m（输入）
L_t_2=2200; % 下部油管（或套管或裸眼）长度，m（输入）
D_t_i_2=0.09718; %0.13970; % 下部油管（或套管或裸眼）内径，m（输入）

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % 油管（或套管或裸眼）内径，m
        else
            D_t_i(t,x)=D_t_i_2; % 油管（或套管或裸眼）内径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(1)
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % 连续油管外径，m
            else
                D_ct_o(t,x)=0; % 连续油管外径，m
            end
        else
            D_ct_o(t,x)=0; % 连续油管外径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(1)
            if Depth(x)<=L_coil(t) % 决定哪几个点可以有内径
                if (L_coil(t)-Depth(x))<=L1 %&& (L_coil(t)-Depth(x))>=0
                    D_ct_i(t,x)=D_ct_i_1; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % 连续油管内径，m
                end
            else
                D_ct_i(t,x)=0; % 连续油管内径，m
            end
        else
            D_ct_i(t,x)=0; % 连续油管内径，m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % 决定哪几个点可以有数值
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % 连续油管内截面积，m^2
        else
            A_ct(t,x)=0; % 连续油管内截面积，m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % 环空水力直径，m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % 环空截面积，m^2
    end
end

%% 数据输入及预处理
h_t=25.4*10^(-6); % 油管（或套管或裸眼）绝对粗糙度，m（输入）
h_ct=25.4*10^(-6); % 连续油管绝对粗糙度，m（输入）
h_a=(h_t+h_ct)/2; % 环空平均绝对粗糙度，m
epsilon_e=1*10^(-3); % 迭代求解误差限（输入）
epsilon_t=1*10^3; % 最大迭代次数（输入）
g=9.81; % 重力加速度，m/s^2（默认）

T_0=20; % 基液测试温度，℃（输入）
P_0=0.1*10^6; % 基液测试压力，Pa（输入）
rho_l_0=1150; % T_0、P_0下基液密度，kg/m^3（输入）
mu_l_0=0.03; % T_0、P_0下基液粘度，Pa*s（输入）
Qv_l_0=0.1/60; % 基液体积流量，m^3/s（输入）
Qm_l_0=Qv_l_0*rho_l_0; % 基液质量流量，kg/s
rho_g_0=0.655; % 注入气密度，kg/m^3
mu_g_0=RheologyG(T_0,P_0); % 气体粘度，Pa*s
Qv_g_0=4/60; % 注入气体积流量，m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % 注入气质量流量，kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % 泡沫质量流量，kg/s

D_s=1*10^(-3); % 砂砾直径，m（输入）
rho_s=2000; % 砂砾密度，kg/m^3（输入）
H_s=L_s_b-L_s_t; % 底部砂床高度，m
PHI=0.6; % 砂床充盈度（输入）
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % 井底砂砾总质量，kg

D_nozzle=4/1000; % 喷嘴直径，m（输入）
N_nozzle=3; % 喷嘴个数（输入）
C=0.95; % 喷嘴流量系数，取0.95（输入）

C0=1.2; % 漂移流方程系数（默认）

M_s(1)=0; % 井底进砂量，kg/s
for t=2:1:Nt
    M_s(t)=PHI*rho_s*V1*0.25*pi*D_t_i(t,Nx_Dsp_pen(1)-Nt+t)^2; % 井底进砂量，kg/s
end

OutPressure=1*10^6; % 井口压力，Pa（输入）

%% 温度设置
T_i=20; % 泡沫注入温度，℃
T_g=0.02; % 地温梯度，℃/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % 连续管内泡沫温度（假设每个时刻都一样），℃
    end
end
T_a=T_ct; % 环空泡沫温度（假设每个时刻都一样），℃

%% 第1个时间节点（初始时刻）相关参数计算（环空）
% 第1个空间节点（井口）处相关参数计算
P_a(1,1)=OutPressure; % 环空压力，Pa
rho_g_a(1,1)=DensityG(T_a(1,1),P_a(1,1)); % 环空气体密度，kg/m^3
rho_l_a(1,1)=DensityL(rho_l_0,T_0,P_0,T_a(1,1),P_a(1,1)); % 环空基液密度，kg/m^3
mu_g_a(1,1)=RheologyG(T_a(1,1),P_a(1,1)); % 环空气体粘度，Pa*s
mu_l_a(1,1)=RheologyL(mu_l_0,T_0,P_0,T_a(1,1),P_a(1,1)); % 环空基液粘度，Pa*s
alpha_f_a(1,1)=1; % 泡沫含量
gamma_g_a(1,1)=(Qm_g_0/rho_g_a(1,1))/(Qm_g_0/rho_g_a(1,1)+Qm_l_0/rho_l_a(1,1)); % 泡沫质量
gamma_l_a(1,1)=1-gamma_g_a(1,1); % 液体滞留量
alpha_g_a(1,1)=alpha_f_a(1,1)*gamma_g_a(1,1); % 环空气体含量
alpha_l_a(1,1)=alpha_f_a(1,1)*gamma_l_a(1,1); % 环空基液含量
alpha_s(1,1)=0; % 固相体积含量
rho_f_a(1,1)=rho_g_a(1,1)*gamma_g_a(1,1)+rho_l_a(1,1)*gamma_l_a(1,1); % 环空泡沫密度，kg/m^3
mu_f_a(1,1)=mu_g_a(1,1)*gamma_g_a(1,1)+mu_l_a(1,1)*gamma_l_a(1,1); % 环空泡沫粘度，Pa*s
V_s(1,1)=0; % 固相速度，m/s
Va_s(1,1)=0; % 固相表观流速，m/s
mu_s(1,1)=mu_f_a(1,1); % 固相粘度，Pa*s
Vsr(1,1)=0; % 砂砾沉降末速，m/s
Va_f_a(1,1)=Qm_f_0/(rho_f_a(1,1)*A_a(1,1)); % 环空泡沫表观流速，m/s
Va_g_a(1,1)=Va_f_a(1,1); % 环空气体表观流速，m/s
Va_l_a(1,1)=Va_f_a(1,1); % 环空基液表观流速，m/s
V_f_a(1,1)=Va_f_a(1,1)/alpha_f_a(1,1); % 环空泡沫流速，m/s
V_g_a(1,1)=V_f_a(1,1); % 环空气体流速，m/s
V_l_a(1,1)=V_f_a(1,1); % 环空基液流速，m/s

V_m(1,1)=Va_s(1,1)+Va_f_a(1,1); % 环空混合物速度，m/s
rho_m(1,1)=alpha_s(1,1)*rho_s+alpha_f_a(1,1)*rho_f_a(1,1); % 环空混合物密度，kg/m^3
mu_m(1,1)=alpha_s(1,1)*mu_s(1,1)+alpha_f_a(1,1)*mu_f_a(1,1); % 环空混合物粘度，Pa*s
[Ff_a(1,1),f_a(1,1),Re_a(1,1),flow_pattern_a(1,1)]=Friction_annulus(rho_m(1,1),V_m(1,1),mu_m(1,1),D_h(1,1),h_a,rho_f_a(1,1),V_f_a(1,1)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态

% 第2～Nx_L_s_t个空间节点处相关参数计算
for x=2:1:Nx_L_s_t
    P_a_ass(1,x)=P_a(1,x-1)+rho_f_a(1,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
    
    err_AnnPressure=1; % 环空压力相对误差
    COUNT_AnnPressure=0; % 环空压力迭代次数初值
    while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
        COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
        
        rho_g_a(1,x)=DensityG(T_a(1,x),P_a_ass(1,x)); % 环空气体密度，kg/m^3
        rho_l_a(1,x)=DensityL(rho_l_0,T_0,P_0,T_a(1,x),P_a_ass(1,x)); % 环空基液密度，kg/m^3
        mu_g_a(1,x)=RheologyG(T_a(1,x),P_a_ass(1,x)); % 环空气体粘度，Pa*s
        mu_l_a(1,x)=RheologyL(mu_l_0,T_0,P_0,T_a(1,x),P_a_ass(1,x)); % 环空基液粘度，Pa*s
        alpha_f_a(1,x)=1; % 环空泡沫含量
        gamma_g_a(1,x)=(Qm_g_0/rho_g_a(1,x))/(Qm_g_0/rho_g_a(1,x)+Qm_l_0/rho_l_a(1,x)); % 泡沫质量
        gamma_l_a(1,x)=1-gamma_g_a(1,x); % 液体滞留量
        alpha_g_a(1,x)=alpha_f_a(1,x)*gamma_g_a(1,x); % 环空气体含量
        alpha_l_a(1,x)=alpha_f_a(1,x)*gamma_l_a(1,x); % 环空基液含量
        alpha_s(1,x)=0; % 固相体积含量
        rho_f_a(1,x)=rho_g_a(1,x)*gamma_g_a(1,x)+rho_l_a(1,x)*gamma_l_a(1,x); % 环空泡沫密度，kg/m^3
        mu_f_a(1,x)=mu_g_a(1,x)*gamma_g_a(1,x)+mu_l_a(1,x)*gamma_l_a(1,x); % 环空泡沫粘度，Pa*s
        V_s(1,x)=0; % 固相速度，m/s
        Va_s(1,x)=0; % 固相表观流速，m/s
        mu_s(1,x)=mu_f_a(1,x); % 固相粘度，Pa*s
        Vsr(1,x)=0; % 砂砾沉降末速，m/s
        Va_f_a(1,x)=Qm_f_0/(rho_f_a(1,x)*A_a(1,x)); % 环空泡沫表观流速，m/s
        Va_g_a(1,x)=Va_f_a(1,x); % 环空气体表观流速，m/s
        Va_l_a(1,x)=Va_f_a(1,x); % 环空基液表观流速，m/s
        V_f_a(1,x)=Va_f_a(1,x)/alpha_f_a(1,x); % 环空泡沫流速，m/s
        V_g_a(1,x)=V_f_a(1,x); % 环空气体流速，m/s
        V_l_a(1,x)=V_f_a(1,x); % 环空基液流速，m/s
        
        V_m(1,x)=Va_s(1,x)+Va_f_a(1,x); % 环空混合物速度，m/s
        rho_m(1,x)=alpha_s(1,x)*rho_s+alpha_f_a(1,x)*rho_f_a(1,x); % 环空混合物密度，kg/m^3
        mu_m(1,x)=alpha_s(1,x)*mu_s(1,x)+alpha_f_a(1,x)*mu_f_a(1,x); % 环空混合物粘度，Pa*s
        [Ff_a(1,x),f_a(1,x),Re_a(1,x),flow_pattern_a(1,x)]=Friction_annulus(rho_m(1,x),V_m(1,x),mu_m(1,x),D_h(1,x),h_a,rho_f_a(1,x),V_f_a(1,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
        P_a(1,x)=-rho_m(1,x)*V_m(1,x)^2+P_a(1,x-1)+rho_m(1,x-1)*V_m(1,x-1)^2+((rho_m(1,x)*g*cosd(theta(x))+Ff_a(1,x)+rho_m(1,x-1)*g*cosd(theta(x-1))+Ff_a(1,x-1))*dx(x-1))/2; % 环空压力，Pa
        
        err_AnnPressure=abs(P_a(1,x)-P_a_ass(1,x))/P_a_ass(1,x); % 计算环空压力假设值与计算值之间的相对误差
        P_a_ass(1,x)=P_a(1,x); % 新的环空压力假设值，Pa
    end
end

% 第Nx_L_s_t+1～Nx个空间节点处相关参数计算
for x=Nx_L_s_t+1:1:Nx
    P_a_ass(1,x)=P_a(1,x-1)+rho_f_a(1,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
    
    err_AnnPressure=1; % 环空压力相对误差
    COUNT_AnnPressure=0; % 环空压力迭代次数初值
    while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
        COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
        
        rho_g_a(1,x)=DensityG(T_a(1,x),P_a_ass(1,x)); % 环空气体密度，kg/m^3
        rho_l_a(1,x)=DensityL(rho_l_0,T_0,P_0,T_a(1,x),P_a_ass(1,x)); % 环空基液密度，kg/m^3
        mu_g_a(1,x)=RheologyG(T_a(1,x),P_a_ass(1,x)); % 环空气体粘度，Pa*s
        mu_l_a(1,x)=RheologyL(mu_l_0,T_0,P_0,T_a(1,x),P_a_ass(1,x)); % 环空基液粘度，Pa*s
        alpha_f_a(1,x)=1-PHI; % 环空泡沫含量
        gamma_g_a(1,x)=(Qm_g_0/rho_g_a(1,x))/(Qm_g_0/rho_g_a(1,x)+Qm_l_0/rho_l_a(1,x)); % 泡沫质量
        gamma_l_a(1,x)=1-gamma_g_a(1,x); % 液体滞留量
        alpha_g_a(1,x)=alpha_f_a(1,x)*gamma_g_a(1,x); % 环空气体含量
        alpha_l_a(1,x)=alpha_f_a(1,x)*gamma_l_a(1,x); % 环空基液含量
        alpha_s(1,x)=PHI; % 固相体积含量
        rho_f_a(1,x)=rho_g_a(1,x)*gamma_g_a(1,x)+rho_l_a(1,x)*gamma_l_a(1,x); % 环空泡沫密度，kg/m^3
        mu_f_a(1,x)=mu_g_a(1,x)*gamma_g_a(1,x)+mu_l_a(1,x)*gamma_l_a(1,x); % 环空泡沫粘度，Pa*s
        V_s(1,x)=0; % 固相速度，m/s
        Va_s(1,x)=0; % 固相表观流速，m/s
        mu_s(1,x)=mu_f_a(1,x); % 固相粘度，Pa*s
        Vsr(1,x)=0; % 砂砾沉降末速，m/s
        Va_f_a(1,x)=0;%Qm_f_0/(rho_f_a(1,x)*A_a(1,x)); % 环空泡沫表观流速，m/s
        Va_g_a(1,x)=Va_f_a(1,x); % 环空气体表观流速，m/s
        Va_l_a(1,x)=Va_f_a(1,x); % 环空基液表观流速，m/s
        V_f_a(1,x)=Va_f_a(1,x)/alpha_f_a(1,x); % 环空泡沫流速，m/s
        V_g_a(1,x)=V_f_a(1,x); % 环空气体流速，m/s
        V_l_a(1,x)=V_f_a(1,x); % 环空基液流速，m/s
        
        V_m(1,x)=Va_s(1,x)+Va_f_a(1,x); % 环空混合物速度，m/s
        rho_m(1,x)=alpha_s(1,x)*rho_s+alpha_f_a(1,x)*rho_f_a(1,x); % 环空混合物密度，kg/m^3
        mu_m(1,x)=alpha_s(1,x)*mu_s(1,x)+alpha_f_a(1,x)*mu_f_a(1,x); % 环空混合物粘度，Pa*s
        [Ff_a(1,x),f_a(1,x),Re_a(1,x),flow_pattern_a(1,x)]=Friction_annulus(rho_m(1,x),V_m(1,x),mu_m(1,x),D_h(1,x),h_a,rho_f_a(1,x),V_f_a(1,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
        P_a(1,x)=-rho_m(1,x)*V_m(1,x)^2+P_a(1,x-1)+rho_m(1,x-1)*V_m(1,x-1)^2+((rho_m(1,x)*g*cosd(theta(x))+Ff_a(1,x)+rho_m(1,x-1)*g*cosd(theta(x-1))+Ff_a(1,x-1))*dx(x-1))/2; % 环空压力，Pa
        
        err_AnnPressure=abs(P_a(1,x)-P_a_ass(1,x))/P_a_ass(1,x); % 计算环空压力假设值与计算值之间的相对误差
        P_a_ass(1,x)=P_a(1,x); % 新的环空压力假设值，Pa
    end
end

%% 第2～Nt个时间节点相关参数计算（环空）
for t=2:1:Nt
    P_a(t,Nx_Dsp_pen(1)-Nt+t)=P_a(t-1,Nx_Dsp_pen(1)-Nt+t);  % 环空管底压力假设值，Pa
    
    err_OutPressure=1; % 出口压力相对误差
    COUNT_OutPressure=0; % 出口压力迭代次数初值
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % 出口压力迭代次数
        
        % 管底（第Nx_Dsp_pen(1)-Nt+t个空间节点）处相关参数计算
        rho_g_a(t,Nx_Dsp_pen(1)-Nt+t)=DensityG(T_a(t,Nx_Dsp_pen(1)-Nt+t),P_a(t,Nx_Dsp_pen(1)-Nt+t)); % 环空气体密度，kg/m^3
        rho_l_a(t,Nx_Dsp_pen(1)-Nt+t)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)-Nt+t),P_a(t,Nx_Dsp_pen(1)-Nt+t)); % 环空基液密度，kg/m^3
        mu_g_a(t,Nx_Dsp_pen(1)-Nt+t)=RheologyG(T_a(t,Nx_Dsp_pen(1)-Nt+t),P_a(t,Nx_Dsp_pen(1)-Nt+t)); % 环空气体粘度，Pa*s
        mu_l_a(t,Nx_Dsp_pen(1)-Nt+t)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)-Nt+t),P_a(t,Nx_Dsp_pen(1)-Nt+t)); % 环空基液粘度，Pa*s
        gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t)=(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)-Nt+t))/(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)-Nt+t)+Qm_l_0/rho_l_a(t,Nx_Dsp_pen(1)-Nt+t)); % 泡沫质量
        gamma_l_a(t,Nx_Dsp_pen(1)-Nt+t)=1-gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t); % 液体滞留量
        alpha_g_a(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_a(t-1,Nx_Dsp_pen(1)-Nt+t)*gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空气体含量
        alpha_l_a(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_a(t-1,Nx_Dsp_pen(1)-Nt+t)*gamma_l_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空基液含量
        rho_f_a(t,Nx_Dsp_pen(1)-Nt+t)=rho_g_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t)+rho_l_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空泡沫密度，kg/m^3
        mu_f_a(t,Nx_Dsp_pen(1)-Nt+t)=mu_g_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t)+mu_l_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空泡沫粘度，Pa*s
        mu_s(t,Nx_Dsp_pen(1)-Nt+t)=mu_f_a(t,Nx_Dsp_pen(1)-Nt+t); % 固相粘度，Pa*s
        Va_f_a(t,Nx_Dsp_pen(1)-Nt+t)=Qm_f_0/(A_a(t,Nx_Dsp_pen(1)-Nt+t)*rho_f_a(t,Nx_Dsp_pen(1)-Nt+t)); % 环空泡沫表观流速，m/s
        Va_g_a(t,Nx_Dsp_pen(1)-Nt+t)=Va_f_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空气体表观流速，m/s
        Va_l_a(t,Nx_Dsp_pen(1)-Nt+t)=Va_f_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空基液表观流速，m/s
        Va_s(t,Nx_Dsp_pen(1)-Nt+t)=M_s(t)/(A_a(t,Nx_Dsp_pen(1)-Nt+t)*rho_s); % 岩屑表观速度，m/s
        Vsr(t,Nx_Dsp_pen(1)-Nt+t)=12*(mu_f_a(t,Nx_Dsp_pen(1)-Nt+t)/(rho_f_a(t,Nx_Dsp_pen(1)-Nt+t)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Dsp_pen(1)-Nt+t))/rho_f_a(t,Nx_Dsp_pen(1)-Nt+t))*((rho_f_a(t,Nx_Dsp_pen(1)-Nt+t)*D_s/mu_f_a(t,Nx_Dsp_pen(1)-Nt+t))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
        alpha_s(t,Nx_Dsp_pen(1)-Nt+t)=Va_s(t,Nx_Dsp_pen(1)-Nt+t)/(C0*(Va_s(t,Nx_Dsp_pen(1)-Nt+t)+Va_f_a(t,Nx_Dsp_pen(1)-Nt+t))-Vsr(t,Nx_Dsp_pen(1)-Nt+t));  % 固相体积分数
        V_s(t,Nx_Dsp_pen(1)-Nt+t)=Va_s(t,Nx_Dsp_pen(1)-Nt+t)/alpha_s(t,Nx_Dsp_pen(1)-Nt+t); % 岩屑速度，m/s
        
        alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t)=1-alpha_s(t,Nx_Dsp_pen(1)-Nt+t); % 环空泡沫含量
        alpha_g_a(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空气体含量
        alpha_l_a(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空基液含量
        V_f_a(t,Nx_Dsp_pen(1)-Nt+t)=Va_f_a(t,Nx_Dsp_pen(1)-Nt+t)/alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空泡沫流速，m/s
        V_g_a(t,Nx_Dsp_pen(1)-Nt+t)=V_f_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空气体流速，m/s
        V_l_a(t,Nx_Dsp_pen(1)-Nt+t)=V_f_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空基液流速，m/s
        
        V_m(t,Nx_Dsp_pen(1)-Nt+t)=Va_s(t,Nx_Dsp_pen(1)-Nt+t)+Va_f_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空混合物速度，m/s
        rho_m(t,Nx_Dsp_pen(1)-Nt+t)=alpha_s(t,Nx_Dsp_pen(1)-Nt+t)*rho_s+alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t)*rho_f_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空混合物密度，kg/m^3
        mu_m(t,Nx_Dsp_pen(1)-Nt+t)=alpha_s(t,Nx_Dsp_pen(1)-Nt+t)*mu_s(t,Nx_Dsp_pen(1)-Nt+t)+alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t)*mu_f_a(t,Nx_Dsp_pen(1)-Nt+t); % 环空混合物粘度，Pa*s
        [Ff_a(t,Nx_Dsp_pen(1)-Nt+t),f_a(t,Nx_Dsp_pen(1)-Nt+t),Re_a(t,Nx_Dsp_pen(1)-Nt+t),flow_pattern_a(t,Nx_Dsp_pen(1)-Nt+t)]=Friction_annulus(rho_m(t,Nx_Dsp_pen(1)-Nt+t),V_m(t,Nx_Dsp_pen(1)-Nt+t),mu_m(t,Nx_Dsp_pen(1)-Nt+t),D_h(t,Nx_Dsp_pen(1)-Nt+t),h_a,rho_f_a(t,Nx_Dsp_pen(1)-Nt+t),V_f_a(t,Nx_Dsp_pen(1)-Nt+t)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
        
        % 第Nx_Dsp_pen(1)-Nt+t-1～1个空间节点处相关参数计算
        for x=Nx_Dsp_pen(1)-Nt+t-1:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % 环空压力假设值，Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液粘度，Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % 环空基液含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            
            err_NodePressure=1; % 环空压力相对误差
            COUNT_NodePressure=0; % 环空压力迭代次数初值
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % 环空压力迭代次数
                
                % 割线法求解固相体积分数
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % 固相体积分数假设值1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % 固相体积分数假设值2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % 固相体积分数绝对误差
                COUNT_NodeEg=0; % 迭代次数
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % 固相体积分数为alpha_s_ass1时（假设值）
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % 固相速度，m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % 固相表观流速，m/s
                    alpha_f_ass1=1-alpha_s_ass1; % 泡沫体积分数
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % 泡沫速度，m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % 泡沫表观流速，m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % 环空混合物速度，m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % 固相体积分数计算值
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 固相体积分数为alpha_s_ass2时（假设值）
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % 固相质量守恒方程离散后公式中间值计算
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % 固相速度，m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % 固相表观流速，m/s
                    alpha_f_ass2=1-alpha_s_ass2; % 泡沫体积分数
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % 泡沫速度，m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % 泡沫表观流速，m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % 环空混合物速度，m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % 固相体积分数计算值
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 割线法计算固相体积分数
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % 新的固相体积分数假设值alpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % 固相体积分数绝对误差
                    alpha_s_ass1=alpha_s_ass2; % 新的固相体积分数假设值1
                    alpha_s_ass2=alpha_s_ass3; % 新的固相体积分数假设值2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % 将迭代求解得到的真实固相体积分数值赋给alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % 当固相体积分数小于一定值时，认为固相体积数为0，用于防止出现后续的计算误差
                    alpha_s(t,x)=0; % 固相体积分数
                    V_s(t,x)=0; % 固相速度，m/s
                    Va_s(t,x)=0; % 固相表观速度，m/s
                    mu_s(t,x)=0; % 固相粘度，Pa*s
                    Vsr(t,x)=0; % 岩屑沉降末速，m/s
                    
                    alpha_f_a(t,x)=1; % 泡沫体积分数                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                else
                    V_s(t,x)=V_s_ass2; % 固相速度，m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % 固相表观流速，m/s
                    Vsr(t,x)=Vsr_ass2; % 岩屑沉降末速，m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % 泡沫密度，kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % 泡沫体积分数
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    V_f_a(t,x)=V_f_ass2; % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % 环空混合物速度，m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % 环空压力计算值，Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % 环空压力相对误差
                P_a(t,x)=P_new; % 新的环空压力假设值，Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % 出口压力相对误差        
        if (P_a(t,1)-OutPressure)>0 % 根据出口压力误差的正负，对环空管底压力假设值进行调节
            P_a(t,Nx_Dsp_pen(1)-Nt+t)=P_a(t,Nx_Dsp_pen(1)-Nt+t)-(P_a(t,1)-OutPressure)/2; % 新的环空管底压力假设值，Pa
        else
            P_a(t,Nx_Dsp_pen(1)-Nt+t)=P_a(t,Nx_Dsp_pen(1)-Nt+t)-(P_a(t,1)-OutPressure)/2*0.3; % 新的环空管底压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(1)-Nt+t+1～Nx个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)-Nt+t+1:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1-PHI; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=PHI; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
end

%% 第1～Nt时间节点钻头压降计算
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Dsp_pen(1)-Nt+t); % 井底泡沫体积流量，m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Dsp_pen(1)-Nt+t)); % 计算冲砂工具压降（Pa）及射流喷嘴流速V_nozzle（m/s）
end

%% 第1～Nt时间节点相关参数计算（连续管内）
for t=1:1:Nt
    % 管底（第Nx_Dsp_pen(1)-Nt+t个空间节点）处相关参数计算
    P_ct(t,Nx_Dsp_pen(1)-Nt+t)=P_a(t,Nx_Dsp_pen(1)-Nt+t)+delta_P_SWT(t); % 管内压力，Pa
    rho_g_ct(t,Nx_Dsp_pen(1)-Nt+t)=DensityG(T_ct(t,Nx_Dsp_pen(1)-Nt+t),P_ct(t,Nx_Dsp_pen(1)-Nt+t)); % 管内气体密度，kg/m^3
    rho_l_ct(t,Nx_Dsp_pen(1)-Nt+t)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)-Nt+t),P_ct(t,Nx_Dsp_pen(1)-Nt+t)); % 管内基液密度，kg/m^3
    mu_g_ct(t,Nx_Dsp_pen(1)-Nt+t)=RheologyG(T_ct(t,Nx_Dsp_pen(1)-Nt+t),P_ct(t,Nx_Dsp_pen(1)-Nt+t)); % 管内气体粘度，Pa*s
    mu_l_ct(t,Nx_Dsp_pen(1)-Nt+t)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)-Nt+t),P_ct(t,Nx_Dsp_pen(1)-Nt+t)); % 管内基液粘度，Pa*s
    alpha_f_ct(t,Nx_Dsp_pen(1)-Nt+t)=1; % 管内泡沫含量
    gamma_g_ct(t,Nx_Dsp_pen(1)-Nt+t)=(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)-Nt+t))/(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)-Nt+t)+Qm_l_0/rho_l_ct(t,Nx_Dsp_pen(1)-Nt+t)); % 管内泡沫质量
    gamma_l_ct(t,Nx_Dsp_pen(1)-Nt+t)=1-gamma_g_ct(t,Nx_Dsp_pen(1)-Nt+t); % 管内液体滞留量
    alpha_g_ct(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_ct(t,Nx_Dsp_pen(1)-Nt+t); % 管内气体含量
    alpha_l_ct(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_ct(t,Nx_Dsp_pen(1)-Nt+t); % 管内基液含量
    rho_f_ct(t,Nx_Dsp_pen(1)-Nt+t)=rho_g_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_ct(t,Nx_Dsp_pen(1)-Nt+t)+rho_l_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_ct(t,Nx_Dsp_pen(1)-Nt+t); % 管内泡沫密度，kg/m^3
    mu_f_ct(t,Nx_Dsp_pen(1)-Nt+t)=mu_g_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_ct(t,Nx_Dsp_pen(1)-Nt+t)+mu_l_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_ct(t,Nx_Dsp_pen(1)-Nt+t); % 管内泡沫粘度，Pa*s
    V_f_ct(t,Nx_Dsp_pen(1)-Nt+t)=Qm_f_0/(rho_f_ct(t,Nx_Dsp_pen(1)-Nt+t)*A_ct(t,Nx_Dsp_pen(1)-Nt+t)); % 管内气体流速，m/s
    [Ff_ct(t,Nx_Dsp_pen(1)-Nt+t),f_ct(t,Nx_Dsp_pen(1)-Nt+t),Re_ct(t,Nx_Dsp_pen(1)-Nt+t),flow_pattern_ct(t,Nx_Dsp_pen(1)-Nt+t)]=Friction_coiledtubing(rho_f_ct(t,Nx_Dsp_pen(1)-Nt+t),V_f_ct(t,Nx_Dsp_pen(1)-Nt+t),mu_f_ct(t,Nx_Dsp_pen(1)-Nt+t),D_ct_i(t,Nx_Dsp_pen(1)-Nt+t),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    % 第Nx_Dsp_pen(1)-Nt+t-1～1个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)-Nt+t-1:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % 管内压力假设值，Pa
        
        err_DriPipePressure=1; % 管内压力相对误差
        COUNT_DriPipePressure=0; % 管内压力迭代次数初值
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % 管内压力迭代次数
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体密度，kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液密度，kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体粘度，Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液粘度，Pa*s
            alpha_f_ct(t,x)=1; % 管内泡沫含量
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % 管内泡沫质量
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % 管内气体含量
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % 管内基液含量
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫密度，kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫粘度，Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % 管内泡沫流速，m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % 管内压力，Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % 计算管内压力假设值与计算值之间的相对误差
            P_ct_ass(t,x)=P_ct(t,x); % 新的管内压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(1)-Nt+t+1～Nx个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)-Nt+t+1:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % 管内气体含量
        alpha_l_ct(t,x)=alpha_l_a(t,x); % 管内基液含量
        gamma_g_ct(t,x)=gamma_g_a(t,x); % 管内泡沫质量
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
        rho_f_ct(t,x)=rho_m(t,x); % 管内气体密度，kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % 管内气体粘度，Pa*s
        V_f_ct(t,t)=V_m(t,x); % 管内气体流速，sm/s
        Ff_ct(t,x)=Ff_a(t,x); % 管内流体单位长度摩擦压降（Pa/m）
        f_ct(t,x)=f_a(t,x); % 管内流体范宁摩擦因子
        Re_ct(t,x)=Re_a(t,x); % 管内流体雷诺数、流体流态
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % 管内流体流态
        P_ct(t,x)=P_a(t,x); % 管内压力，Pa
    end
end

%% 第1～Nt时间节点盘管段出口压力及泵压计算（连续管内）
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % 地面管内流体流速，m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % 地面管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % 盘管段流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % 导向器内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % 井口到注入头顶部段摩擦压降，Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % 导向器段摩擦压降，Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % 导向器到滚筒段摩擦压降，Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % 盘管段摩擦压降，Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % 盘管段出口压力，Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % 泵压，Pa
end

%% 计算井口累积出砂量
M_w_tem(1)=0; % 井口瞬时出砂量初值，kg
M_w_tot(1)=0; % 井口累积出砂量初值，kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口瞬时出砂量，kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口累积出砂量，kg
end

%% 计算末态井筒最大砂浓度
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(1)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% 计算环空泡沫返速平均值V_l_mre
V_f_mre=0; % 环空返速平均值，m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % 环空返速平均值，m/s
    end
end

%% 计算固相沉降末速平均值Vsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % 从Vsr矩阵中提取值非零的沉降末速，m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % 固相沉降末速平均值，m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % 固相沉降末速平均值，m/s
end

%% 判断是否有效冲砂（环空返速平均值大于2倍砂砾沉降末速平均值）、是否完成冲砂
fprintf("冲洗钻进过程1：\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“有效冲砂”
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % 命令行窗口说明当前条件下“完成冲砂”
    else
        fprintf("Sand Cleanout UnFinished!\n"); % 命令行窗口说明当前条件下“未完成冲砂”
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“无效冲砂”
end

%% 井筒ECD计算
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % 环空ECD，kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % 管内ECD，kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % 环空ECD，kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % 管内ECD，kg/m^3
end

%% 冲洗钻进过程末状态数据
alpha_f_a_1=alpha_f_a(Nt,:); % 冲洗钻进过程末状态泡沫含量
alpha_g_a_1=alpha_g_a(Nt,:); % 冲洗钻进过程末状态气相含量
alpha_l_a_1=alpha_l_a(Nt,:); % 冲洗钻进过程末状态液相含量
alpha_s_1=alpha_s(Nt,:); % 冲洗钻进过程末状态固相含量
f_a_1=f_a(Nt,:); % 冲洗钻进过程末状态环空摩擦因子
Ff_a_1=Ff_a(Nt,:); % 冲洗钻进过程末状态环空单位长度摩擦压降，Pa/m
flow_pattern_a_1=flow_pattern_a(Nt,:); % 冲洗钻进过程末状态环空流体流态
gamma_g_a_1=gamma_g_a(Nt,:); % 冲洗钻进过程末状态泡沫质量
gamma_l_a_1=gamma_l_a(Nt,:); % 冲洗钻进过程末状态液体滞留量
mu_f_a_1=mu_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫粘度，Pa*s
mu_g_a_1=mu_g_a(Nt,:); % 冲洗钻进过程末状态气相粘度，Pa*s
mu_l_a_1=mu_l_a(Nt,:); % 冲洗钻进过程末状态环空液体粘度，Pa*s
mu_m_1=mu_m(Nt,:); % 冲洗钻进过程末状态环空混合物粘度，Pa*s
mu_s_1=mu_s(Nt,:); % 冲洗钻进过程末状态固相粘度，Pa*s
P_a_1=P_a(Nt,:); % 冲洗钻进过程末状态环空压力，Pa
Re_a_1=Re_a(Nt,:); % 冲洗钻进过程末状态环空雷诺数
rho_f_a_1=rho_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫密度，kg/m^3
rho_g_a_1=rho_g_a(Nt,:); % 冲洗钻进过程末状态环空气相密度，kg/m^3
rho_l_a_1=rho_l_a(Nt,:); % 冲洗钻进过程末状态环空液相密度，kg/m^3
rho_m_1=rho_m(Nt,:); % 冲洗钻进过程末状态环空混合物密度，kg/m^3
V_f_a_1=V_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫流速，m/s
V_g_a_1=V_g_a(Nt,:); % 冲洗钻进过程末状态气相流速，m/s
V_l_a_1=V_l_a(Nt,:); % 冲洗钻进过程末状态环空液相流速，m/s
V_m_1=V_m(Nt,:); % 冲洗钻进过程末状态环空混合物流速，m/s
V_s_1=V_s(Nt,:); % 冲洗钻进过程末状态固相流速，m/s
Va_f_a_1=Va_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫表观流速，m/s
Va_g_a_1=Va_g_a(Nt,:); % 冲洗钻进过程末状态气相表观流速，m/s
Va_l_a_1=Va_l_a(Nt,:); % 冲洗钻进过程末状态环空液相表观流速，m/s
Va_s_1=Va_s(Nt,:); % 冲洗钻进过程末状态固相表观流速，m/s
Vsr_1=Vsr(Nt,:); % 冲洗钻进过程末状态固相滑移速度，m/s

%% 数据存储（冲洗钻进过程）
ANS_Nt_1=Nt; % 时间节点数
ANS_Time_1=Time(1:Nt); % 冲洗时间，s
ANS_alpha_g_a_1=alpha_g_a(1:Nt,:); % 环空气相体积分数
ANS_alpha_l_a_1=alpha_l_a(1:Nt,:); % 环空液相体积分数
ANS_alpha_s_1=alpha_s(1:Nt,:); % 固相含量
ANS_P_ct_1=P_ct(1:Nt,:); % 管内压力（Pa）
ANS_delta_P_SWT_1=delta_P_SWT(1:Nt); % 冲砂工具压降（Pa）
ANS_P_a_1=P_a(1:Nt,:); % 环空压力（Pa）
ANS_T_a_1=T_a(1:Nt,:); % 环空温度（℃）
ANS_T_ct_1=T_ct(1:Nt,:); % 管内温度（℃）
ANS_P_coil_1=P_coil(1:Nt); % 盘管段出口压力（Pa）
ANS_P_pump_1=P_pump(1:Nt); % 泵压（Pa）
ANS_M_w_tem_1=M_w_tem(1:Nt); % 井口瞬时出砂量（kg）
ANS_M_w_tot_1=M_w_tot(1:Nt); % 井口累积出砂量（kg）
ANS_Va_s_1=Va_s(1:Nt,:); % 岩屑表观速度（m/s）
ANS_Va_f_a_1=Va_f_a(1:Nt,:); % 环空泡沫表观流速（m/s）
ANS_V_s_1=V_s(1:Nt,:); % 岩屑沉降速度（m/s）
ANS_V_f_a_1=V_f_a(1:Nt,:); % 环空泡沫返速（m/s）
ANS_alpha_g_ct_1=alpha_g_ct(1:Nt,:); % 管内气相体积分数
ANS_alpha_l_ct_1=alpha_l_ct(1:Nt,:); % 管内液相体积分数
ANS_gamma_g_ct_1=gamma_g_ct(1:Nt,:); % 管内泡沫质量
ANS_gamma_g_a_1=gamma_g_a(1:Nt,:); % 环空泡沫质量
ANS_ECD_a_1=ECD_a(1:Nt,:); % 环空ECD（kg/m^3）
ANS_L_coil_1=L_coil(1:Nt); % 连续管下深（m）



%% 回拖短起过程1（国际单位制）
%% 变量数值清零
dt=zeros(); % 时间步长，s
dx=zeros(); % 空间步长，m
Time=zeros(); % 连续管下入总时长，s
L_coil=zeros(); % 连续管下深，m
L_reel=zeros(); % 盘管段长度，m
D_t_i=zeros(); % 油管内径，m
D_ct_o=zeros(); % 连续油管外径，m
D_ct_i=zeros(); % 连续油管内径，m
M_s=zeros(); % 井底进砂量，kg/s
P_a=zeros(); % 环空压力，Pa
rho_g_a=zeros(); % 气相密度，kg/m^3
rho_l_a=zeros(); % 环空液相密度，kg/m^3
rho_f_a=zeros(); % 环空泡沫密度，kg/m^3
mu_g_a=zeros(); % 气相粘度，Pa*s
mu_l_a=zeros(); % 环空液体粘度，Pa*s
mu_f_a=zeros(); % 环空泡沫粘度，Pa*s
mu_s=zeros(); % 固相粘度，Pa*s
alpha_g_a=zeros(); % 气相含量
alpha_l_a=zeros(); % 液相含量
alpha_f_a=zeros(); % 泡沫含量
alpha_s=zeros(); % 固相含量
Va_g_a=zeros(); % 气相表观流速，m/s
Va_l_a=zeros(); % 环空液相表观流速，m/s
Va_f_a=zeros(); % 环空泡沫表观流速，m/s
Va_s=zeros(); % 固相表观流速，m/s
V_g_a=zeros(); % 气相流速，m/s
V_l_a=zeros(); % 环空液相流速，m/s
V_f_a=zeros(); % 环空泡沫流速，m/s
V_s=zeros(); % 固相流速，m/s
Vsr=zeros(); % 固相滑移速度，m/s
gamma_g_a=zeros(); % 泡沫质量
gamma_l_a=zeros(); % 液体滞留量
V_m=zeros(); % 环空混合物流速，m/s
rho_m=zeros(); % 环空混合物密度，kg/m^3
mu_m=zeros(); % 环空混合物粘度，Pa*s
Ff_a=zeros(); % 环空单位长度摩擦压降，Pa/m
f_a=zeros(); % 环空摩擦因子
Re_a=zeros(); % 环空雷诺数
flow_pattern_a=zeros(); % 环空流体流态

%% 计算空间步长dx及相应时间步长dt
V2=(Dsp_pen(1)-Lp)/(3*3600); % 回拖短起速度，m/s（输入）
t_2=(Dsp_pen(1)-Lp)/V2; % 回拖短起总时长，s

Nt=Nx_Dsp_pen(1)-Nx_Lp+1; % 时间节点数
nt=Nt-1; % 时间网格数

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % 空间步长，m
end

Time(1)=0; % 连续管起出总时长初值，s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Dsp_pen(1)-t)/V2; % 每起出一个空间步长所需时间，s
    Time(t+1)=Time(t)+dt(t); % 连续管起出至第(t+1)个空间节点所经历的总时长，s
end

%% 计算不同时刻连续管起出深度L_coil及盘管段连续管长度L_reel
L=10000; % 连续油管总长，m（输入）
L_wg=8; % 井口到注入头顶部段连续管长度，m（输入）
L_goose=3; % 导向器段连续管长度，m（输入）
D_goose=2; % 导向器段半径，m（输入）
H_goose=10; % 导向器距地面高度，m（输入）
L_gr=20; % 导向器到滚筒段连续管长度，m（输入）
theta_gr=acosd(H_goose/L_gr); % 导向器到滚筒段连续管与铅垂线夹角，°
D_r_i=3; % 滚筒内径，m（输入）
D_r_o=5; % 滚筒外径，m（输入）
W_r=5; % 滚筒宽度，m（输入）
D_cable=0.005; % 电缆外径，m

L_coil(1)=Dsp_pen(1); % 初始时刻连续管底部深度，m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % 初始时刻盘管段连续管长度，m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)-dx(Nx_Dsp_pen(1)-t+1); % 连续管起出深度，m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % 盘管段连续管长度，m
end
L_cable=L_coil; % 电缆下入长度，m

%% 数据输入及预处理
D_ct_o_0=0.04445; % 连续油管外径，m（输入）
L1=2000; % 近出口第一段连续管长度，m（输入）
D_ct_i_1=0.03709; % 近出口第一段连续管内径，m（输入）
L2=2000; % 近出口第二段连续管长度，m（输入）
D_ct_i_2=0.03653; % 近出口第二段连续管内径，m（输入）
L3=2000; % 近出口第三段连续管长度，m（输入）
D_ct_i_3=0.03555; % 近出口第三段连续管内径，m（输入）
L4=2000; % 近出口第四段连续管长度，m（输入）
D_ct_i_4=0.03489; % 近出口第四段连续管内径，m（输入）
L5=2000; % 近出口第五段连续管长度，m（输入）
D_ct_i_5=0.03409; % 近出口第五段连续管内径，m（输入）

L_t_1=4000; % 上部油管（或套管或裸眼）长度，m（输入）
D_t_i_1=0.09718; %0.068;%0.09718; % 上部油管（或套管或裸眼）内径，m（输入）
L_t_2=2200; % 下部油管（或套管或裸眼）长度，m（输入）
D_t_i_2=0.09718; %0.13970; % 下部油管（或套管或裸眼）内径，m（输入）

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % 油管（或套管或裸眼）内径，m
        else
            D_t_i(t,x)=D_t_i_2; % 油管（或套管或裸眼）内径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(1)
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % 连续油管外径，m
            else
                D_ct_o(t,x)=0; % 连续油管外径，m
            end
        else
            D_ct_o(t,x)=0; % 连续油管外径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(1)
            if Depth(x)<=L_coil(t) % 决定哪几个点可以有内径
                if (L_coil(t)-Depth(x))<=L1 && (L_coil(t)-Depth(x))>=0
                    D_ct_i(t,x)=D_ct_i_1; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % 连续油管内径，m
                end
            else
                D_ct_i(t,x)=0; % 连续油管内径，m
            end
        else
            D_ct_i(t,x)=0; % 连续油管内径，m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % 决定哪几个点可以有数值
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % 连续油管内截面积，m^2
        else
            A_ct(t,x)=0; % 连续油管内截面积，m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % 环空水力直径，m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % 环空截面积，m^2
    end
end

%% 数据输入及预处理
h_t=25.4*10^(-6); % 油管（或套管或裸眼）绝对粗糙度，m（输入）
h_ct=25.4*10^(-6); % 连续油管绝对粗糙度，m（输入）
h_a=(h_t+h_ct)/2; % 环空平均绝对粗糙度，m
epsilon_e=1*10^(-3); % 迭代求解误差限（输入）
epsilon_t=1*10^3; % 最大迭代次数（输入）
g=9.81; % 重力加速度，m/s^2（默认）

T_0=20; % 基液测试温度，℃（输入）
P_0=0.1*10^6; % 基液测试压力，Pa（输入）
rho_l_0=1150; % T_0、P_0下基液密度，kg/m^3（输入）
mu_l_0=0.03; % T_0、P_0下基液粘度，Pa*s（输入）
Qv_l_0=0.1/60; % 基液体积流量，m^3/s（输入）
Qm_l_0=Qv_l_0*rho_l_0; % 基液质量流量，kg/s
rho_g_0=0.655; % 注入气密度，kg/m^3
mu_g_0=RheologyG(T_0,P_0); % 气体粘度，Pa*s
Qv_g_0=4/60; % 注入气体积流量，m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % 注入气质量流量，kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % 泡沫质量流量，kg/s

D_s=1*10^(-3); % 砂砾直径，m（输入）
rho_s=2000; % 砂砾密度，kg/m^3（输入）
H_s=L_s_b-Nx_Dsp_pen(1); % 底部砂床高度，m
PHI=0.6; % 砂床充盈度（输入）
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % 井底砂砾总质量，kg

D_nozzle=4/1000; % 喷嘴直径，m（输入）
N_nozzle=3; % 喷嘴个数（输入）
C=0.95; % 喷嘴流量系数，取0.95（输入）

C0=1.2; % 漂移流方程系数（默认）

M_s(1)=0; % 井底进砂量，kg/s
for t=2:1:Nt
    M_s(t)=0; % 井底进砂量，kg/s
end

OutPressure=1*10^6; % 井口压力，Pa（输入）

%% 温度设置
T_i=20; % 泡沫注入温度，℃
T_g=0.02; % 地温梯度，℃/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % 连续管内泡沫温度（假设每个时刻都一样），℃
    end
end
T_a=T_ct; % 环空泡沫温度（假设每个时刻都一样），℃

%% 第1个时间节点（初始时刻）相关参数插值计算（环空）
for x=1:1:Nx
    P_a(1,x)=P_a_1(x); % 环空压力，Pa
    rho_g_a(1,x)=rho_g_a_1(x); % 气相密度，kg/m^3
    rho_l_a(1,x)=rho_l_a_1(x); % 环空液相密度，kg/m^3
    rho_f_a(1,x)=rho_f_a_1(x); % 环空泡沫密度，kg/m^3
    mu_g_a(1,x)=mu_g_a_1(x); % 气相粘度，Pa*s
    mu_l_a(1,x)=mu_l_a_1(x); % 环空液体粘度，Pa*s
    mu_f_a(1,x)=mu_f_a_1(x); % 环空泡沫粘度，Pa*s
    mu_s(1,x)=mu_s_1(x); % 固相粘度，Pa*s
    alpha_g_a(1,x)=alpha_g_a_1(x); % 气相含量
    alpha_l_a(1,x)=alpha_l_a_1(x); % 液相含量
    alpha_f_a(1,x)=alpha_f_a_1(x); % 泡沫含量
    alpha_s(1,x)=alpha_s_1(1,x); % 固相含量
    Va_g_a(1,x)=Va_g_a_1(x); % 气相表观流速，m/s
    Va_l_a(1,x)=Va_l_a_1(x); % 环空液相表观流速，m/s
    Va_f_a(1,x)=Va_f_a_1(x); % 环空泡沫表观流速，m/s
    Va_s(1,x)=Va_s_1(x); % 固相表观流速，m/s
    V_g_a(1,x)=V_g_a_1(x); % 气相流速，m/s
    V_l_a(1,x)=V_l_a_1(x); % 环空液相流速，m/s
    V_f_a(1,x)=V_f_a_1(x); % 环空泡沫流速，m/s
    V_s(1,x)=V_s_1(x); % 固相流速，m/s
    Vsr(1,x)=Vsr_1(x); % 固相滑移速度，m/s
    gamma_g_a(1,x)=gamma_g_a_1(x); % 泡沫质量
    gamma_l_a(1,x)=gamma_l_a_1(1,x); % 液体滞留量
    V_m(1,x)=V_m_1(x); % 环空混合物流速，m/s
    rho_m(1,x)=rho_m_1(x); % 环空混合物密度，kg/m^3
    mu_m(1,x)=mu_m_1(x); % 环空混合物粘度，Pa*s
    Ff_a(1,x)=Ff_a_1(x); % 环空单位长度摩擦压降，Pa/m
    f_a(1,x)=f_a_1(x); % 环空摩擦因子
    Re_a(1,x)=Re_a_1(x); % 环空雷诺数
    flow_pattern_a(1,x)=flow_pattern_a_1(x); % 环空流体流态
end

%% 第2～Nt个时间节点相关参数计算（环空）
for t=2:1:Nt
    P_a(t,Nx_Dsp_pen(1)-t+1)=P_a(t-1,Nx_Dsp_pen(1)-t+1);  % 环空管底压力假设值，Pa
    
    err_OutPressure=1; % 出口压力相对误差
    COUNT_OutPressure=0; % 出口压力迭代次数初值
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % 出口压力迭代次数
        
        % 管底（第Nx_Dsp_pen(1)-t+1个空间节点）处相关参数计算
        rho_g_a(t,Nx_Dsp_pen(1)-t+1)=DensityG(T_a(t,Nx_Dsp_pen(1)-t+1),P_a(t,Nx_Dsp_pen(1)-t+1)); % 环空气体密度，kg/m^3
        rho_l_a(t,Nx_Dsp_pen(1)-t+1)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)-t+1),P_a(t,Nx_Dsp_pen(1)-t+1)); % 环空基液密度，kg/m^3
        mu_g_a(t,Nx_Dsp_pen(1)-t+1)=RheologyG(T_a(t,Nx_Dsp_pen(1)-t+1),P_a(t,Nx_Dsp_pen(1)-t+1)); % 环空气体粘度，Pa*s
        mu_l_a(t,Nx_Dsp_pen(1)-t+1)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)-t+1),P_a(t,Nx_Dsp_pen(1)-t+1)); % 环空基液粘度，Pa*s
        gamma_g_a(t,Nx_Dsp_pen(1)-t+1)=(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)-t+1))/(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)-t+1)+Qm_l_0/rho_l_a(t,Nx_Dsp_pen(1)-t+1)); % 泡沫质量
        gamma_l_a(t,Nx_Dsp_pen(1)-t+1)=1-gamma_g_a(t,Nx_Dsp_pen(1)-t+1); % 液体滞留量
        alpha_g_a(t,Nx_Dsp_pen(1)-t+1)=alpha_f_a(t-1,Nx_Dsp_pen(1)-t+1)*gamma_g_a(t,Nx_Dsp_pen(1)-t+1); % 环空气体含量
        alpha_l_a(t,Nx_Dsp_pen(1)-t+1)=alpha_f_a(t-1,Nx_Dsp_pen(1)-t+1)*gamma_l_a(t,Nx_Dsp_pen(1)-t+1); % 环空基液含量
        rho_f_a(t,Nx_Dsp_pen(1)-t+1)=rho_g_a(t,Nx_Dsp_pen(1)-t+1)*gamma_g_a(t,Nx_Dsp_pen(1)-t+1)+rho_l_a(t,Nx_Dsp_pen(1)-t+1)*gamma_l_a(t,Nx_Dsp_pen(1)-t+1); % 环空泡沫密度，kg/m^3
        mu_f_a(t,Nx_Dsp_pen(1)-t+1)=mu_g_a(t,Nx_Dsp_pen(1)-t+1)*gamma_g_a(t,Nx_Dsp_pen(1)-t+1)+mu_l_a(t,Nx_Dsp_pen(1)-t+1)*gamma_l_a(t,Nx_Dsp_pen(1)-t+1); % 环空泡沫粘度，Pa*s
        mu_s(t,Nx_Dsp_pen(1)-t+1)=mu_f_a(t,Nx_Dsp_pen(1)-t+1); % 固相粘度，Pa*s
        Va_f_a(t,Nx_Dsp_pen(1)-t+1)=Qm_f_0/(A_a(t,Nx_Dsp_pen(1)-t+1)*rho_f_a(t,Nx_Dsp_pen(1)-t+1)); % 环空泡沫表观流速，m/s
        Va_g_a(t,Nx_Dsp_pen(1)-t+1)=Va_f_a(t,Nx_Dsp_pen(1)-t+1); % 环空气体表观流速，m/s
        Va_l_a(t,Nx_Dsp_pen(1)-t+1)=Va_f_a(t,Nx_Dsp_pen(1)-t+1); % 环空基液表观流速，m/s
        Va_s(t,Nx_Dsp_pen(1)-t+1)=M_s(t)/(A_a(t,Nx_Dsp_pen(1)-t+1)*rho_s); % 岩屑表观速度，m/s
        Vsr(t,Nx_Dsp_pen(1)-t+1)=12*(mu_f_a(t,Nx_Dsp_pen(1)-t+1)/(rho_f_a(t,Nx_Dsp_pen(1)-t+1)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Dsp_pen(1)-t+1))/rho_f_a(t,Nx_Dsp_pen(1)-t+1))*((rho_f_a(t,Nx_Dsp_pen(1)-t+1)*D_s/mu_f_a(t,Nx_Dsp_pen(1)-t+1))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
        alpha_s(t,Nx_Dsp_pen(1)-t+1)=Va_s(t,Nx_Dsp_pen(1)-t+1)/(C0*(Va_s(t,Nx_Dsp_pen(1)-t+1)+Va_f_a(t,Nx_Dsp_pen(1)-t+1))-Vsr(t,Nx_Dsp_pen(1)-t+1));  % 固相体积分数
        V_s(t,Nx_Dsp_pen(1)-t+1)=0;%Va_s(t,Nx_Dsp_pen(1)-t+1)/alpha_s(t,Nx_Dsp_pen(1)-t+1); % 岩屑速度，m/s
        
        alpha_f_a(t,Nx_Dsp_pen(1)-t+1)=1-alpha_s(t,Nx_Dsp_pen(1)-t+1); % 环空泡沫含量
        alpha_g_a(t,Nx_Dsp_pen(1)-t+1)=alpha_f_a(t,Nx_Dsp_pen(1)-t+1)*gamma_g_a(t,Nx_Dsp_pen(1)-t+1); % 环空气体含量
        alpha_l_a(t,Nx_Dsp_pen(1)-t+1)=alpha_f_a(t,Nx_Dsp_pen(1)-t+1)*gamma_l_a(t,Nx_Dsp_pen(1)-t+1); % 环空基液含量
        V_f_a(t,Nx_Dsp_pen(1)-t+1)=Va_f_a(t,Nx_Dsp_pen(1)-t+1)/alpha_f_a(t,Nx_Dsp_pen(1)-t+1); % 环空泡沫流速，m/s
        V_g_a(t,Nx_Dsp_pen(1)-t+1)=V_f_a(t,Nx_Dsp_pen(1)-t+1); % 环空气体流速，m/s
        V_l_a(t,Nx_Dsp_pen(1)-t+1)=V_f_a(t,Nx_Dsp_pen(1)-t+1); % 环空基液流速，m/s
        
        V_m(t,Nx_Dsp_pen(1)-t+1)=Va_s(t,Nx_Dsp_pen(1)-t+1)+Va_f_a(t,Nx_Dsp_pen(1)-t+1); % 环空混合物速度，m/s
        rho_m(t,Nx_Dsp_pen(1)-t+1)=alpha_s(t,Nx_Dsp_pen(1)-t+1)*rho_s+alpha_f_a(t,Nx_Dsp_pen(1)-t+1)*rho_f_a(t,Nx_Dsp_pen(1)-t+1); % 环空混合物密度，kg/m^3
        mu_m(t,Nx_Dsp_pen(1)-t+1)=alpha_s(t,Nx_Dsp_pen(1)-t+1)*mu_s(t,Nx_Dsp_pen(1)-t+1)+alpha_f_a(t,Nx_Dsp_pen(1)-t+1)*mu_f_a(t,Nx_Dsp_pen(1)-t+1); % 环空混合物粘度，Pa*s
        [Ff_a(t,Nx_Dsp_pen(1)-t+1),f_a(t,Nx_Dsp_pen(1)-t+1),Re_a(t,Nx_Dsp_pen(1)-t+1),flow_pattern_a(t,Nx_Dsp_pen(1)-t+1)]=Friction_annulus(rho_m(t,Nx_Dsp_pen(1)-t+1),V_m(t,Nx_Dsp_pen(1)-t+1),mu_m(t,Nx_Dsp_pen(1)-t+1),D_h(t,Nx_Dsp_pen(1)-t+1),h_a,rho_f_a(t,Nx_Dsp_pen(1)-t+1),V_f_a(t,Nx_Dsp_pen(1)-t+1)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
        
        % 第Nx_Dsp_pen(1)-t～1个空间节点处相关参数计算
        for x=Nx_Dsp_pen(1)-t:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % 环空压力假设值，Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液粘度，Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % 环空基液含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            
            err_NodePressure=1; % 环空压力相对误差
            COUNT_NodePressure=0; % 环空压力迭代次数初值
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % 环空压力迭代次数
                
                % 割线法求解固相体积分数
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % 固相体积分数假设值1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % 固相体积分数假设值2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % 固相体积分数绝对误差
                COUNT_NodeEg=0; % 迭代次数
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % 固相体积分数为alpha_s_ass1时（假设值）
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % 固相速度，m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % 固相表观流速，m/s
                    alpha_f_ass1=1-alpha_s_ass1; % 泡沫体积分数
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % 泡沫速度，m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % 泡沫表观流速，m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % 环空混合物速度，m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % 固相体积分数计算值
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 固相体积分数为alpha_s_ass2时（假设值）
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % 固相质量守恒方程离散后公式中间值计算
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % 固相速度，m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % 固相表观流速，m/s
                    alpha_f_ass2=1-alpha_s_ass2; % 泡沫体积分数
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % 泡沫速度，m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % 泡沫表观流速，m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % 环空混合物速度，m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % 固相体积分数计算值
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 割线法计算固相体积分数
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % 新的固相体积分数假设值alpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % 固相体积分数绝对误差
                    alpha_s_ass1=alpha_s_ass2; % 新的固相体积分数假设值1
                    alpha_s_ass2=alpha_s_ass3; % 新的固相体积分数假设值2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % 将迭代求解得到的真实固相体积分数值赋给alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % 当固相体积分数小于一定值时，认为固相体积数为0，用于防止出现后续的计算误差
                    alpha_s(t,x)=0; % 固相体积分数
                    V_s(t,x)=0; % 固相速度，m/s
                    Va_s(t,x)=0; % 固相表观速度，m/s
                    mu_s(t,x)=0; % 固相粘度，Pa*s
                    Vsr(t,x)=0; % 岩屑沉降末速，m/s
                    
                    alpha_f_a(t,x)=1; % 泡沫体积分数                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                else
                    V_s(t,x)=V_s_ass2; % 固相速度，m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % 固相表观流速，m/s
                    Vsr(t,x)=Vsr_ass2; % 岩屑沉降末速，m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % 泡沫密度，kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % 泡沫体积分数
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    V_f_a(t,x)=V_f_ass2; % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % 环空混合物速度，m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % 环空压力计算值，Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % 环空压力相对误差
                P_a(t,x)=P_new; % 新的环空压力假设值，Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % 出口压力相对误差
        if (P_a(t,1)-OutPressure)>0 % 根据出口压力误差的正负，对环空管底压力假设值进行调节
            P_a(t,Nx_Dsp_pen(1)-t+1)=P_a(t,Nx_Dsp_pen(1)-t+1)-(P_a(t,1)-OutPressure)/2; % 新的环空管底压力假设值，Pa
        else
            P_a(t,Nx_Dsp_pen(1)-t+1)=P_a(t,Nx_Dsp_pen(1)-t+1)-(P_a(t,1)-OutPressure)/2*0.3; % 新的环空管底压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(1)-t+2～Nx_Dsp_pen(1)个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)-t+2:1:Nx_Dsp_pen(1)
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=0; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(1)+1～Nx个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)+1:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1-PHI; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=PHI; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
end

%% 第1～Nt时间节点钻头压降计算
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Dsp_pen(1)-t+1); % 井底泡沫体积流量，m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Dsp_pen(1)-t+1)); % 计算冲砂工具压降（Pa）及射流喷嘴流速V_nozzle（m/s）
end

%% 第1～Nt时间节点相关参数计算（连续管内）
for t=1:1:Nt
    % 管底（第Nx_Dsp_pen(1)-t+1个空间节点）处相关参数计算
    P_ct(t,Nx_Dsp_pen(1)-t+1)=P_a(t,Nx_Dsp_pen(1)-t+1)+delta_P_SWT(t); % 管内压力，Pa
    rho_g_ct(t,Nx_Dsp_pen(1)-t+1)=DensityG(T_ct(t,Nx_Dsp_pen(1)-t+1),P_ct(t,Nx_Dsp_pen(1)-t+1)); % 管内气体密度，kg/m^3
    rho_l_ct(t,Nx_Dsp_pen(1)-t+1)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)-t+1),P_ct(t,Nx_Dsp_pen(1)-t+1)); % 管内基液密度，kg/m^3
    mu_g_ct(t,Nx_Dsp_pen(1)-t+1)=RheologyG(T_ct(t,Nx_Dsp_pen(1)-t+1),P_ct(t,Nx_Dsp_pen(1)-t+1)); % 管内气体粘度，Pa*s
    mu_l_ct(t,Nx_Dsp_pen(1)-t+1)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)-t+1),P_ct(t,Nx_Dsp_pen(1)-t+1)); % 管内基液粘度，Pa*s
    alpha_f_ct(t,Nx_Dsp_pen(1)-t+1)=1; % 管内泡沫含量
    gamma_g_ct(t,Nx_Dsp_pen(1)-t+1)=(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)-t+1))/(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)-t+1)+Qm_l_0/rho_l_ct(t,Nx_Dsp_pen(1)-t+1)); % 管内泡沫质量
    gamma_l_ct(t,Nx_Dsp_pen(1)-t+1)=1-gamma_g_ct(t,Nx_Dsp_pen(1)-t+1); % 管内液体滞留量
    alpha_g_ct(t,Nx_Dsp_pen(1)-t+1)=alpha_f_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(1)-t+1); % 管内气体含量
    alpha_l_ct(t,Nx_Dsp_pen(1)-t+1)=alpha_f_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(1)-t+1); % 管内基液含量
    rho_f_ct(t,Nx_Dsp_pen(1)-t+1)=rho_g_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(1)-t+1)+rho_l_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(1)-t+1); % 管内泡沫密度，kg/m^3
    mu_f_ct(t,Nx_Dsp_pen(1)-t+1)=mu_g_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(1)-t+1)+mu_l_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(1)-t+1); % 管内泡沫粘度，Pa*s
    V_f_ct(t,Nx_Dsp_pen(1)-t+1)=Qm_f_0/(rho_f_ct(t,Nx_Dsp_pen(1)-t+1)*A_ct(t,Nx_Dsp_pen(1)-t+1)); % 管内气体流速，m/s
    [Ff_ct(t,Nx_Dsp_pen(1)-t+1),f_ct(t,Nx_Dsp_pen(1)-t+1),Re_ct(t,Nx_Dsp_pen(1)-t+1),flow_pattern_ct(t,Nx_Dsp_pen(1)-t+1)]=Friction_coiledtubing(rho_f_ct(t,Nx_Dsp_pen(1)-t+1),V_f_ct(t,Nx_Dsp_pen(1)-t+1),mu_f_ct(t,Nx_Dsp_pen(1)-t+1),D_ct_i(t,Nx_Dsp_pen(1)-t+1),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    % 第Nx_Dsp_pen(1)-t～1个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)-t:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % 管内压力假设值，Pa
        
        err_DriPipePressure=1; % 管内压力相对误差
        COUNT_DriPipePressure=0; % 管内压力迭代次数初值
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % 管内压力迭代次数
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体密度，kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液密度，kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体粘度，Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液粘度，Pa*s
            alpha_f_ct(t,x)=1; % 管内泡沫含量
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % 管内泡沫质量
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % 管内气体含量
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % 管内基液含量
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫密度，kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫粘度，Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % 管内泡沫流速，m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % 管内压力，Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % 计算管内压力假设值与计算值之间的相对误差
            P_ct_ass(t,x)=P_ct(t,x); % 新的管内压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(1)-t+2～Nx个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)-t+2:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % 管内气体含量
        alpha_l_ct(t,x)=alpha_l_a(t,x); % 管内基液含量
        gamma_g_ct(t,x)=gamma_g_a(t,x); % 管内泡沫质量
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
        rho_f_ct(t,x)=rho_m(t,x); % 管内气体密度，kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % 管内气体粘度，Pa*s
        V_f_ct(t,t)=V_m(t,x); % 管内气体流速，sm/s
        Ff_ct(t,x)=Ff_a(t,x); % 管内流体单位长度摩擦压降（Pa/m）
        f_ct(t,x)=f_a(t,x); % 管内流体范宁摩擦因子
        Re_ct(t,x)=Re_a(t,x); % 管内流体雷诺数、流体流态
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % 管内流体流态
        P_ct(t,x)=P_a(t,x); % 管内压力，Pa
    end
end

%% 第1～Nt时间节点盘管段出口压力及泵压计算（连续管内）
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % 地面管内流体流速，m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % 地面管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % 盘管段流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % 导向器内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % 井口到注入头顶部段摩擦压降，Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % 导向器段摩擦压降，Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % 导向器到滚筒段摩擦压降，Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % 盘管段摩擦压降，Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % 盘管段出口压力，Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % 泵压，Pa
end

%% 计算井口累积出砂量
M_w_tem(1)=0; % 井口瞬时出砂量初值，kg
M_w_tot(1)=0; % 井口累积出砂量初值，kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口瞬时出砂量，kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口累积出砂量，kg
end

%% 计算末态井筒最大砂浓度
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(1)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% 计算环空泡沫返速平均值V_l_mre
V_f_mre=0; % 环空返速平均值，m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % 环空返速平均值，m/s
    end
end

%% 计算固相沉降末速平均值Vsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % 从Vsr矩阵中提取值非零的沉降末速，m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % 固相沉降末速平均值，m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % 固相沉降末速平均值，m/s
end

%% 判断是否有效冲砂（环空返速平均值大于2倍砂砾沉降末速平均值）、是否完成冲砂
fprintf("回拖短起过程1：\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“有效冲砂”
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % 命令行窗口说明当前条件下“完成冲砂”
    else
        fprintf("Sand Cleanout UnFinished!\n"); % 命令行窗口说明当前条件下“未完成冲砂”
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“无效冲砂”
end

%% 井筒ECD计算
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % 环空ECD，kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % 管内ECD，kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % 环空ECD，kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % 管内ECD，kg/m^3
end

%% 回拖短起过程末状态数据
alpha_f_a_2=alpha_f_a(Nt,:); % 冲洗钻进过程末状态泡沫含量
alpha_g_a_2=alpha_g_a(Nt,:); % 冲洗钻进过程末状态气相含量
alpha_l_a_2=alpha_l_a(Nt,:); % 冲洗钻进过程末状态液相含量
alpha_s_2=alpha_s(Nt,:); % 冲洗钻进过程末状态固相含量
f_a_2=f_a(Nt,:); % 冲洗钻进过程末状态环空摩擦因子
Ff_a_2=Ff_a(Nt,:); % 冲洗钻进过程末状态环空单位长度摩擦压降，Pa/m
flow_pattern_a_2=flow_pattern_a(Nt,:); % 冲洗钻进过程末状态环空流体流态
gamma_g_a_2=gamma_g_a(Nt,:); % 冲洗钻进过程末状态泡沫质量
gamma_l_a_2=gamma_l_a(Nt,:); % 冲洗钻进过程末状态液体滞留量
mu_f_a_2=mu_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫粘度，Pa*s
mu_g_a_2=mu_g_a(Nt,:); % 冲洗钻进过程末状态气相粘度，Pa*s
mu_l_a_2=mu_l_a(Nt,:); % 冲洗钻进过程末状态环空液体粘度，Pa*s
mu_m_2=mu_m(Nt,:); % 冲洗钻进过程末状态环空混合物粘度，Pa*s
mu_s_2=mu_s(Nt,:); % 冲洗钻进过程末状态固相粘度，Pa*s
P_a_2=P_a(Nt,:); % 冲洗钻进过程末状态环空压力，Pa
Re_a_2=Re_a(Nt,:); % 冲洗钻进过程末状态环空雷诺数
rho_f_a_2=rho_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫密度，kg/m^3
rho_g_a_2=rho_g_a(Nt,:); % 冲洗钻进过程末状态环空气相密度，kg/m^3
rho_l_a_2=rho_l_a(Nt,:); % 冲洗钻进过程末状态环空液相密度，kg/m^3
rho_m_2=rho_m(Nt,:); % 冲洗钻进过程末状态环空混合物密度，kg/m^3
V_f_a_2=V_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫流速，m/s
V_g_a_2=V_g_a(Nt,:); % 冲洗钻进过程末状态气相流速，m/s
V_l_a_2=V_l_a(Nt,:); % 冲洗钻进过程末状态环空液相流速，m/s
V_m_2=V_m(Nt,:); % 冲洗钻进过程末状态环空混合物流速，m/s
V_s_2=V_s(Nt,:); % 冲洗钻进过程末状态固相流速，m/s
Va_f_a_2=Va_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫表观流速，m/s
Va_g_a_2=Va_g_a(Nt,:); % 冲洗钻进过程末状态气相表观流速，m/s
Va_l_a_2=Va_l_a(Nt,:); % 冲洗钻进过程末状态环空液相表观流速，m/s
Va_s_2=Va_s(Nt,:); % 冲洗钻进过程末状态固相表观流速，m/s
Vsr_2=Vsr(Nt,:); % 冲洗钻进过程末状态固相滑移速度，m/s

%% 数据存储（回拖短起过程）
ANS_Nt_2=Nt; % 时间节点数
ANS_Time_2=Time(1:Nt); % 连续管起出时间，s
ANS_alpha_g_a_2=alpha_g_a(1:Nt,:); % 环空气相体积分数
ANS_alpha_l_a_2=alpha_l_a(1:Nt,:); % 环空液相体积分数
ANS_alpha_s_2=alpha_s(1:Nt,:); % 固相含量
ANS_P_ct_2=P_ct(1:Nt,:); % 管内压力（Pa）
ANS_delta_P_SWT_2=delta_P_SWT(1:Nt); % 冲砂工具压降（Pa）
ANS_P_a_2=P_a(1:Nt,:); % 环空压力（Pa）
ANS_T_a_2=T_a(1:Nt,:); % 环空温度（℃）
ANS_T_ct_2=T_ct(1:Nt,:); % 管内温度（℃）
ANS_P_coil_2=P_coil(1:Nt); % 盘管段出口压力（Pa）
ANS_P_pump_2=P_pump(1:Nt); % 泵压（Pa）
ANS_M_w_tem_2=M_w_tem(1:Nt); % 井口瞬时出砂量（kg）
ANS_M_w_tot_2=M_w_tot(1:Nt); % 井口累积出砂量（kg）
ANS_Va_s_2=Va_s(1:Nt,:); % 岩屑表观速度（m/s）
ANS_Va_f_a_2=Va_f_a(1:Nt,:); % 环空泡沫表观流速（m/s）
ANS_V_s_2=V_s(1:Nt,:); % 岩屑沉降速度（m/s）
ANS_V_f_a_2=V_f_a(1:Nt,:); % 环空泡沫返速（m/s）
ANS_alpha_g_ct_2=alpha_g_ct(1:Nt,:); % 管内气相体积分数
ANS_alpha_l_ct_2=alpha_l_ct(1:Nt,:); % 管内液相体积分数
ANS_gamma_g_ct_2=gamma_g_ct(1:Nt,:); % 管内泡沫质量
ANS_gamma_g_a_2=gamma_g_a(1:Nt,:); % 环空泡沫质量
ANS_ECD_a_2=ECD_a(1:Nt,:); % 环空ECD（kg/m^3）
ANS_L_coil_2=L_coil(1:Nt); % 连续管下深（m）



%% 定点循环过程1（国际单位制）
%% 变量数值清零
dt=zeros(); % 时间步长，s
dx=zeros(); % 空间步长，m
Time=zeros(); % 连续管下入总时长，s
L_coil=zeros(); % 连续管下深，m
L_reel=zeros(); % 盘管段长度，m
D_t_i=zeros(); % 油管内径，m
D_ct_o=zeros(); % 连续油管外径，m
D_ct_i=zeros(); % 连续油管内径，m
M_s=zeros(); % 井底进砂量，kg/s
P_a=zeros(); % 环空压力，Pa
rho_g_a=zeros(); % 气相密度，kg/m^3
rho_l_a=zeros(); % 环空液相密度，kg/m^3
rho_f_a=zeros(); % 环空泡沫密度，kg/m^3
mu_g_a=zeros(); % 气相粘度，Pa*s
mu_l_a=zeros(); % 环空液体粘度，Pa*s
mu_f_a=zeros(); % 环空泡沫粘度，Pa*s
mu_s=zeros(); % 固相粘度，Pa*s
alpha_g_a=zeros(); % 气相含量
alpha_l_a=zeros(); % 液相含量
alpha_f_a=zeros(); % 泡沫含量
alpha_s=zeros(); % 固相含量
Va_g_a=zeros(); % 气相表观流速，m/s
Va_l_a=zeros(); % 环空液相表观流速，m/s
Va_f_a=zeros(); % 环空泡沫表观流速，m/s
Va_s=zeros(); % 固相表观流速，m/s
V_g_a=zeros(); % 气相流速，m/s
V_l_a=zeros(); % 环空液相流速，m/s
V_f_a=zeros(); % 环空泡沫流速，m/s
V_s=zeros(); % 固相流速，m/s
Vsr=zeros(); % 固相滑移速度，m/s
gamma_g_a=zeros(); % 泡沫质量
gamma_l_a=zeros(); % 液体滞留量
V_m=zeros(); % 环空混合物流速，m/s
rho_m=zeros(); % 环空混合物密度，kg/m^3
mu_m=zeros(); % 环空混合物粘度，Pa*s
Ff_a=zeros(); % 环空单位长度摩擦压降，Pa/m
f_a=zeros(); % 环空摩擦因子
Re_a=zeros(); % 环空雷诺数
flow_pattern_a=zeros(); % 环空流体流态

%% 计算空间步长dx及相应时间步长dt
t_3=0.5*3600; % 定点循环总时长，s

Nt=101; % 时间节点数（输入）
nt=Nt-1; % 时间网格数

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % 空间步长，m
end

for t=1:1:Nt-1
    dt(t)=t_3/nt; % 时间步长
end

Time(1)=0; % 定点循环总时长，s
for t=1:1:Nt-1
    Time(t+1)=Time(t)+dt(t); % 定点循环至第(t+1)个时刻所经历的总时长，s
end

%% 计算不同时刻连续管下入深度L_coil及盘管段连续管长度L_reel
L=10000; % 连续油管总长，m（输入）
L_wg=8; % 井口到注入头顶部段连续管长度，m（输入）
L_goose=3; % 导向器段连续管长度，m（输入）
D_goose=2; % 导向器段半径，m（输入）
H_goose=10; % 导向器距地面高度，m（输入）
L_gr=20; % 导向器到滚筒段连续管长度，m（输入）
theta_gr=acosd(H_goose/L_gr); % 导向器到滚筒段连续管与铅垂线夹角，°
D_r_i=3; % 滚筒内径，m（输入）
D_r_o=5; % 滚筒外径，m（输入）
W_r=5; % 滚筒宽度，m（输入）
D_cable=0.005; % 电缆外径，m

L_coil(1)=Lp; % 连续管下入深度，m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % 初始时刻盘管段连续管长度，m
for t=2:1:Nt
    L_coil(t)=Lp; % 连续管下入深度，m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % 盘管段连续管长度，m
end
L_cable=L_coil; % 电缆下入长度，m

%% 数据输入及预处理
D_ct_o_0=0.04445; % 连续油管外径，m（输入）
L1=2000; % 近出口第一段连续管长度，m（输入）
D_ct_i_1=0.03709; % 近出口第一段连续管内径，m（输入）
L2=2000; % 近出口第二段连续管长度，m（输入）
D_ct_i_2=0.03653; % 近出口第二段连续管内径，m（输入）
L3=2000; % 近出口第三段连续管长度，m（输入）
D_ct_i_3=0.03555; % 近出口第三段连续管内径，m（输入）
L4=2000; % 近出口第四段连续管长度，m（输入）
D_ct_i_4=0.03489; % 近出口第四段连续管内径，m（输入）
L5=2000; % 近出口第五段连续管长度，m（输入）
D_ct_i_5=0.03409; % 近出口第五段连续管内径，m（输入）

L_t_1=4000; % 上部油管（或套管或裸眼）长度，m（输入）
D_t_i_1=0.09718; %0.068;% 上部油管（或套管或裸眼）内径，m（输入）
L_t_2=2200; % 下部油管（或套管或裸眼）长度，m（输入）
D_t_i_2=0.09718; %0.13970; % 下部油管（或套管或裸眼）内径，m（输入）

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % 油管（或套管或裸眼）内径，m
        else
            D_t_i(t,x)=D_t_i_2; % 油管（或套管或裸眼）内径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            D_ct_o(t,x)=D_ct_o_0; % 连续油管外径，m
        else
            D_ct_o(t,x)=0; % 连续油管外径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            if (L_coil(t)-Depth(x))<=L1 && (L_coil(t)-Depth(x))>=0
                D_ct_i(t,x)=D_ct_i_1; % 连续油管内径，m
            elseif (L_coil(t)-Depth(x))<=L1+L2
                D_ct_i(t,x)=D_ct_i_2; % 连续油管内径，m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                D_ct_i(t,x)=D_ct_i_3; % 连续油管内径，m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                D_ct_i(t,x)=D_ct_i_4; % 连续油管内径，m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                D_ct_i(t,x)=D_ct_i_5; % 连续油管内径，m
            end
        else
            D_ct_i(t,x)=0; % 连续油管内径，m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % 决定哪几个点可以有数值
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % 连续油管内截面积，m^2
        else
            A_ct(t,x)=0; % 连续油管内截面积，m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % 环空水力直径，m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % 环空截面积，m^2
    end
end

%% 数据输入及预处理
h_t=25.4*10^(-6); % 油管（或套管或裸眼）绝对粗糙度，m（输入）
h_ct=25.4*10^(-6); % 连续油管绝对粗糙度，m（输入）
h_a=(h_t+h_ct)/2; % 环空平均绝对粗糙度，m
epsilon_e=1*10^(-3); % 迭代求解误差限（输入）
epsilon_t=1*10^3; % 最大迭代次数（输入）
g=9.81; % 重力加速度，m/s^2（默认）

T_0=20; % 基液测试温度，℃（输入）
P_0=0.1*10^6; % 基液测试压力，Pa（输入）
rho_l_0=1150; % T_0、P_0下基液密度，kg/m^3（输入）
mu_l_0=0.03; % T_0、P_0下基液粘度，Pa*s（输入）
Qv_l_0=0.1/60; % 基液体积流量，m^3/s（输入）
Qm_l_0=Qv_l_0*rho_l_0; % 基液质量流量，kg/s
rho_g_0=0.655; % 注入气密度，kg/m^3
mu_g_0=RheologyG(T_0,P_0); % 气体粘度，Pa*s
Qv_g_0=4/60; % 注入气体积流量，m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % 注入气质量流量，kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % 泡沫质量流量，kg/s

D_s=1*10^(-3); % 砂砾直径，m（输入）
rho_s=2000; % 砂砾密度，kg/m^3（输入）
H_s=L_s_b-Dsp_pen(1); % 底部砂床高度，m
PHI=0.6; % 砂床充盈度（输入）
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % 井底砂砾总质量，kg

D_nozzle=4/1000; % 喷嘴直径，m（输入）
N_nozzle=3; % 喷嘴个数（输入）
C=0.95; % 喷嘴流量系数，取0.95（输入）

C0=1.2; % 漂移流方程系数（默认）

M_s(1)=0; % 井底进砂量，kg/s
for t=2:1:Nt
    M_s(t)=0; % 井底进砂量，kg/s
end

OutPressure=1*10^6; % 井口压力，Pa（输入）

%% 温度设置
T_i=20; % 泡沫注入温度，℃
T_g=0.02; % 地温梯度，℃/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % 连续管内泡沫温度（假设每个时刻都一样），℃
    end
end
T_a=T_ct; % 环空泡沫温度（假设每个时刻都一样），℃

%% 第1个时间节点（初始时刻）相关参数插值计算（环空）
for x=1:1:Nx
    P_a(1,x)=P_a_2(x); % 环空压力，Pa
    rho_g_a(1,x)=rho_g_a_2(x); % 气相密度，kg/m^3
    rho_l_a(1,x)=rho_l_a_2(x); % 环空液相密度，kg/m^3
    rho_f_a(1,x)=rho_f_a_2(x); % 环空泡沫密度，kg/m^3
    mu_g_a(1,x)=mu_g_a_2(x); % 气相粘度，Pa*s
    mu_l_a(1,x)=mu_l_a_2(x); % 环空液体粘度，Pa*s
    mu_f_a(1,x)=mu_f_a_2(x); % 环空泡沫粘度，Pa*s
    mu_s(1,x)=mu_s_2(x); % 固相粘度，Pa*s
    alpha_g_a(1,x)=alpha_g_a_2(x); % 气相含量
    alpha_l_a(1,x)=alpha_l_a_2(x); % 液相含量
    alpha_f_a(1,x)=alpha_f_a_2(x); % 泡沫含量
    alpha_s(1,x)=alpha_s_2(1,x); % 固相含量
    Va_g_a(1,x)=Va_g_a_2(x); % 气相表观流速，m/s
    Va_l_a(1,x)=Va_l_a_2(x); % 环空液相表观流速，m/s
    Va_f_a(1,x)=Va_f_a_2(x); % 环空泡沫表观流速，m/s
    Va_s(1,x)=Va_s_2(x); % 固相表观流速，m/s
    V_g_a(1,x)=V_g_a_2(x); % 气相流速，m/s
    V_l_a(1,x)=V_l_a_2(x); % 环空液相流速，m/s
    V_f_a(1,x)=V_f_a_2(x); % 环空泡沫流速，m/s
    V_s(1,x)=V_s_2(x); % 固相流速，m/s
    Vsr(1,x)=Vsr_2(x); % 固相滑移速度，m/s
    gamma_g_a(1,x)=gamma_g_a_2(x); % 泡沫质量
    gamma_l_a(1,x)=gamma_l_a_2(1,x); % 液体滞留量
    V_m(1,x)=V_m_2(x); % 环空混合物流速，m/s
    rho_m(1,x)=rho_m_2(x); % 环空混合物密度，kg/m^3
    mu_m(1,x)=mu_m_2(x); % 环空混合物粘度，Pa*s
    Ff_a(1,x)=Ff_a_2(x); % 环空单位长度摩擦压降，Pa/m
    f_a(1,x)=f_a_2(x); % 环空摩擦因子
    Re_a(1,x)=Re_a_2(x); % 环空雷诺数
    flow_pattern_a(1,x)=flow_pattern_a_2(x); % 环空流体流态
end

%% 第2～Nt个时间节点相关参数计算（环空）
for t=2:1:Nt
    P_a(t,Nx_Lp)=P_a(t-1,Nx_Lp);  % 环空管底压力假设值，Pa
    
    err_OutPressure=1; % 出口压力相对误差
    COUNT_OutPressure=0; % 出口压力迭代次数初值
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % 出口压力迭代次数
        
        % 管底（第Nx_Lp个空间节点）处相关参数计算
        rho_g_a(t,Nx_Lp)=DensityG(T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % 环空气体密度，kg/m^3
        rho_l_a(t,Nx_Lp)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % 环空基液密度，kg/m^3
        mu_g_a(t,Nx_Lp)=RheologyG(T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % 环空气体粘度，Pa*s
        mu_l_a(t,Nx_Lp)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % 环空基液粘度，Pa*s
        gamma_g_a(t,Nx_Lp)=(Qm_g_0/rho_g_a(t,Nx_Lp))/(Qm_g_0/rho_g_a(t,Nx_Lp)+Qm_l_0/rho_l_a(t,Nx_Lp)); % 泡沫质量
        gamma_l_a(t,Nx_Lp)=1-gamma_g_a(t,Nx_Lp); % 液体滞留量
        alpha_g_a(t,Nx_Lp)=alpha_f_a(t-1,Nx_Lp)*gamma_g_a(t,Nx_Lp); % 环空气体含量
        alpha_l_a(t,Nx_Lp)=alpha_f_a(t-1,Nx_Lp)*gamma_l_a(t,Nx_Lp); % 环空基液含量
        rho_f_a(t,Nx_Lp)=rho_g_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp)+rho_l_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % 环空泡沫密度，kg/m^3
        mu_f_a(t,Nx_Lp)=mu_g_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp)+mu_l_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % 环空泡沫粘度，Pa*s
        mu_s(t,Nx_Lp)=mu_f_a(t,Nx_Lp); % 固相粘度，Pa*s
        Va_f_a(t,Nx_Lp)=Qm_f_0/(A_a(t,Nx_Lp)*rho_f_a(t,Nx_Lp)); % 环空泡沫表观流速，m/s
        Va_g_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp); % 环空气体表观流速，m/s
        Va_l_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp); % 环空基液表观流速，m/s
        Va_s(t,Nx_Lp)=M_s(t)/(A_a(t,Nx_Lp)*rho_s); % 岩屑表观速度，m/s
        Vsr(t,Nx_Lp)=12*(mu_f_a(t,Nx_Lp)/(rho_f_a(t,Nx_Lp)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Lp))/rho_f_a(t,Nx_Lp))*((rho_f_a(t,Nx_Lp)*D_s/mu_f_a(t,Nx_Lp))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
        alpha_s(t,Nx_Lp)=Va_s(t,Nx_Lp)/(C0*(Va_s(t,Nx_Lp)+Va_f_a(t,Nx_Lp))-Vsr(t,Nx_Lp));  % 固相体积分数
        V_s(t,Nx_Lp)=0;%Va_s(t,Nx_Lp)/alpha_s(t,Nx_Lp); % 岩屑速度，m/s
        
        alpha_f_a(t,Nx_Lp)=1-alpha_s(t,Nx_Lp); % 环空泡沫含量
        alpha_g_a(t,Nx_Lp)=alpha_f_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp); % 环空气体含量
        alpha_l_a(t,Nx_Lp)=alpha_f_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % 环空基液含量
        V_f_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp)/alpha_f_a(t,Nx_Lp); % 环空泡沫流速，m/s
        V_g_a(t,Nx_Lp)=V_f_a(t,Nx_Lp); % 环空气体流速，m/s
        V_l_a(t,Nx_Lp)=V_f_a(t,Nx_Lp); % 环空基液流速，m/s
        
        V_m(t,Nx_Lp)=Va_s(t,Nx_Lp)+Va_f_a(t,Nx_Lp); % 环空混合物速度，m/s
        rho_m(t,Nx_Lp)=alpha_s(t,Nx_Lp)*rho_s+alpha_f_a(t,Nx_Lp)*rho_f_a(t,Nx_Lp); % 环空混合物密度，kg/m^3
        mu_m(t,Nx_Lp)=alpha_s(t,Nx_Lp)*mu_s(t,Nx_Lp)+alpha_f_a(t,Nx_Lp)*mu_f_a(t,Nx_Lp); % 环空混合物粘度，Pa*s
        [Ff_a(t,Nx_Lp),f_a(t,Nx_Lp),Re_a(t,Nx_Lp),flow_pattern_a(t,Nx_Lp)]=Friction_annulus(rho_m(t,Nx_Lp),V_m(t,Nx_Lp),mu_m(t,Nx_Lp),D_h(t,Nx_Lp),h_a,rho_f_a(t,Nx_Lp),V_f_a(t,Nx_Lp)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
        
        % 第Nx_Lp-1～1个空间节点处相关参数计算
        for x=Nx_Lp-1:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % 环空压力假设值，Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液粘度，Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % 环空基液含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            
            err_NodePressure=1; % 环空压力相对误差
            COUNT_NodePressure=0; % 环空压力迭代次数初值
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % 环空压力迭代次数
                
                % 割线法求解固相体积分数
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % 固相体积分数假设值1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % 固相体积分数假设值2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % 固相体积分数绝对误差
                COUNT_NodeEg=0; % 迭代次数
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % 固相体积分数为alpha_s_ass1时（假设值）
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % 固相速度，m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % 固相表观流速，m/s
                    alpha_f_ass1=1-alpha_s_ass1; % 泡沫体积分数
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % 泡沫速度，m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % 泡沫表观流速，m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % 环空混合物速度，m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % 固相体积分数计算值
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 固相体积分数为alpha_s_ass2时（假设值）
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % 固相质量守恒方程离散后公式中间值计算
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % 固相速度，m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % 固相表观流速，m/s
                    alpha_f_ass2=1-alpha_s_ass2; % 泡沫体积分数
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % 泡沫速度，m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % 泡沫表观流速，m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % 环空混合物速度，m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % 固相体积分数计算值
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 割线法计算固相体积分数
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % 新的固相体积分数假设值alpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % 固相体积分数绝对误差
                    alpha_s_ass1=alpha_s_ass2; % 新的固相体积分数假设值1
                    alpha_s_ass2=alpha_s_ass3; % 新的固相体积分数假设值2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % 将迭代求解得到的真实固相体积分数值赋给alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % 当固相体积分数小于一定值时，认为固相体积数为0，用于防止出现后续的计算误差
                    alpha_s(t,x)=0; % 固相体积分数
                    V_s(t,x)=0; % 固相速度，m/s
                    Va_s(t,x)=0; % 固相表观速度，m/s
                    mu_s(t,x)=0; % 固相粘度，Pa*s
                    Vsr(t,x)=0; % 岩屑沉降末速，m/s
                    
                    alpha_f_a(t,x)=1; % 泡沫体积分数                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                else
                    V_s(t,x)=V_s_ass2; % 固相速度，m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % 固相表观流速，m/s
                    Vsr(t,x)=Vsr_ass2; % 岩屑沉降末速，m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % 泡沫密度，kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % 泡沫体积分数
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    V_f_a(t,x)=V_f_ass2; % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % 环空混合物速度，m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % 环空压力计算值，Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % 环空压力相对误差
                P_a(t,x)=P_new; % 新的环空压力假设值，Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % 出口压力相对误差        
        if (P_a(t,1)-OutPressure)>0 % 根据出口压力误差的正负，对环空管底压力假设值进行调节
            P_a(t,Nx_Lp)=P_a(t,Nx_Lp)-(P_a(t,1)-OutPressure)/2; % 新的环空管底压力假设值，Pa
        else
            P_a(t,Nx_Lp)=P_a(t,Nx_Lp)-(P_a(t,1)-OutPressure)/2*0.3; % 新的环空管底压力假设值，Pa
        end
    end
    
    % 第Nx_Lp+1～Nx_Dsp_pen(1)个空间节点处相关参数计算
    for x=Nx_Lp+1:1:Nx_Dsp_pen(1)
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=0; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(1)+1～Nx个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)+1:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1-PHI; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=PHI; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
end

%% 第1～Nt时间节点钻头压降计算
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Lp); % 井底泡沫体积流量，m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Lp)); % 计算冲砂工具压降（Pa）及射流喷嘴流速V_nozzle（m/s）
end

%% 第1～Nt时间节点相关参数计算（连续管内）
for t=1:1:Nt
    % 管底（第Nx_Lp个空间节点）处相关参数计算
    P_ct(t,Nx_Lp)=P_a(t,Nx_Lp)+delta_P_SWT(t); % 管内压力，Pa
    rho_g_ct(t,Nx_Lp)=DensityG(T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % 管内气体密度，kg/m^3
    rho_l_ct(t,Nx_Lp)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % 管内基液密度，kg/m^3
    mu_g_ct(t,Nx_Lp)=RheologyG(T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % 管内气体粘度，Pa*s
    mu_l_ct(t,Nx_Lp)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % 管内基液粘度，Pa*s
    alpha_f_ct(t,Nx_Lp)=1; % 管内泡沫含量
    gamma_g_ct(t,Nx_Lp)=(Qm_g_0/rho_g_ct(t,Nx_Lp))/(Qm_g_0/rho_g_ct(t,Nx_Lp)+Qm_l_0/rho_l_ct(t,Nx_Lp)); % 管内泡沫质量
    gamma_l_ct(t,Nx_Lp)=1-gamma_g_ct(t,Nx_Lp); % 管内液体滞留量
    alpha_g_ct(t,Nx_Lp)=alpha_f_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp); % 管内气体含量
    alpha_l_ct(t,Nx_Lp)=alpha_f_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % 管内基液含量
    rho_f_ct(t,Nx_Lp)=rho_g_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp)+rho_l_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % 管内泡沫密度，kg/m^3
    mu_f_ct(t,Nx_Lp)=mu_g_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp)+mu_l_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % 管内泡沫粘度，Pa*s
    V_f_ct(t,Nx_Lp)=Qm_f_0/(rho_f_ct(t,Nx_Lp)*A_ct(t,Nx_Lp)); % 管内气体流速，m/s
    [Ff_ct(t,Nx_Lp),f_ct(t,Nx_Lp),Re_ct(t,Nx_Lp),flow_pattern_ct(t,Nx_Lp)]=Friction_coiledtubing(rho_f_ct(t,Nx_Lp),V_f_ct(t,Nx_Lp),mu_f_ct(t,Nx_Lp),D_ct_i(t,Nx_Lp),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    % 第Nx_Lp-1～1个空间节点处相关参数计算
    for x=Nx_Lp-1:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % 管内压力假设值，Pa
        
        err_DriPipePressure=1; % 管内压力相对误差
        COUNT_DriPipePressure=0; % 管内压力迭代次数初值
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % 管内压力迭代次数
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体密度，kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液密度，kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体粘度，Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液粘度，Pa*s
            alpha_f_ct(t,x)=1; % 管内泡沫含量
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % 管内泡沫质量
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % 管内气体含量
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % 管内基液含量
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫密度，kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫粘度，Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % 管内泡沫流速，m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % 管内压力，Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % 计算管内压力假设值与计算值之间的相对误差
            P_ct_ass(t,x)=P_ct(t,x); % 新的管内压力假设值，Pa
        end
    end
    
    % 第Nx_Lp+1～Nx个空间节点处相关参数计算
    for x=Nx_Lp+1:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % 管内气体含量
        alpha_l_ct(t,x)=alpha_l_a(t,x); % 管内基液含量
        gamma_g_ct(t,x)=gamma_g_a(t,x); % 管内泡沫质量
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
        rho_f_ct(t,x)=rho_m(t,x); % 管内气体密度，kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % 管内气体粘度，Pa*s
        V_f_ct(t,t)=V_m(t,x); % 管内气体流速，sm/s
        Ff_ct(t,x)=Ff_a(t,x); % 管内流体单位长度摩擦压降（Pa/m）
        f_ct(t,x)=f_a(t,x); % 管内流体范宁摩擦因子
        Re_ct(t,x)=Re_a(t,x); % 管内流体雷诺数、流体流态
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % 管内流体流态
        P_ct(t,x)=P_a(t,x); % 管内压力，Pa
    end
end

%% 第1～Nt时间节点盘管段出口压力及泵压计算（连续管内）
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % 地面管内流体流速，m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % 地面管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % 盘管段流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % 导向器内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % 井口到注入头顶部段摩擦压降，Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % 导向器段摩擦压降，Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % 导向器到滚筒段摩擦压降，Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % 盘管段摩擦压降，Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % 盘管段出口压力，Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % 泵压，Pa
end

%% 计算井口累积出砂量
M_w_tem(1)=0; % 井口瞬时出砂量初值，kg
M_w_tot(1)=0; % 井口累积出砂量初值，kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口瞬时出砂量，kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口累积出砂量，kg
end

%% 计算末态井筒最大砂浓度
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(1)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% 计算环空泡沫返速平均值V_l_mre
V_f_mre=0; % 环空返速平均值，m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % 环空返速平均值，m/s
    end
end

%% 计算固相沉降末速平均值Vsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % 从Vsr矩阵中提取值非零的沉降末速，m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % 固相沉降末速平均值，m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % 固相沉降末速平均值，m/s
end

%% 判断是否有效冲砂（环空返速平均值大于2倍砂砾沉降末速平均值）、是否完成冲砂
fprintf("定点循环过程1：\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“有效冲砂”
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % 命令行窗口说明当前条件下“完成冲砂”
    else
        fprintf("Sand Cleanout UnFinished!\n"); % 命令行窗口说明当前条件下“未完成冲砂”
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“无效冲砂”
end

%% 井筒ECD计算
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % 环空ECD，kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % 管内ECD，kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % 环空ECD，kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % 管内ECD，kg/m^3
end

%% 定点循环过程末状态数据
alpha_f_a_3=alpha_f_a(Nt,:); % 定点循环过程末状态泡沫含量
alpha_g_a_3=alpha_g_a(Nt,:); % 定点循环过程末状态气相含量
alpha_l_a_3=alpha_l_a(Nt,:); % 定点循环过程末状态液相含量
alpha_s_3=alpha_s(Nt,:); % 定点循环过程末状态固相含量
f_a_3=f_a(Nt,:); % 定点循环过程末状态环空摩擦因子
Ff_a_3=Ff_a(Nt,:); % 定点循环过程末状态环空单位长度摩擦压降，Pa/m
flow_pattern_a_3=flow_pattern_a(Nt,:); % 定点循环过程末状态环空流体流态
gamma_g_a_3=gamma_g_a(Nt,:); % 定点循环过程末状态泡沫质量
gamma_l_a_3=gamma_l_a(Nt,:); % 定点循环过程末状态液体滞留量
mu_f_a_3=mu_f_a(Nt,:); % 定点循环过程末状态环空泡沫粘度，Pa*s
mu_g_a_3=mu_g_a(Nt,:); % 定点循环过程末状态气相粘度，Pa*s
mu_l_a_3=mu_l_a(Nt,:); % 定点循环过程末状态环空液体粘度，Pa*s
mu_m_3=mu_m(Nt,:); % 定点循环过程末状态环空混合物粘度，Pa*s
mu_s_3=mu_s(Nt,:); % 定点循环过程末状态固相粘度，Pa*s
P_a_3=P_a(Nt,:); % 定点循环过程末状态环空压力，Pa
Re_a_3=Re_a(Nt,:); % 定点循环过程末状态环空雷诺数
rho_f_a_3=rho_f_a(Nt,:); % 定点循环过程末状态环空泡沫密度，kg/m^3
rho_g_a_3=rho_g_a(Nt,:); % 定点循环过程末状态环空气相密度，kg/m^3
rho_l_a_3=rho_l_a(Nt,:); % 定点循环过程末状态环空液相密度，kg/m^3
rho_m_3=rho_m(Nt,:); % 定点循环过程末状态环空混合物密度，kg/m^3
V_f_a_3=V_f_a(Nt,:); % 定点循环过程末状态环空泡沫流速，m/s
V_g_a_3=V_g_a(Nt,:); % 定点循环过程末状态气相流速，m/s
V_l_a_3=V_l_a(Nt,:); % 定点循环过程末状态环空液相流速，m/s
V_m_3=V_m(Nt,:); % 定点循环过程末状态环空混合物流速，m/s
V_s_3=V_s(Nt,:); % 定点循环过程末状态固相流速，m/s
Va_f_a_3=Va_f_a(Nt,:); % 定点循环过程末状态环空泡沫表观流速，m/s
Va_g_a_3=Va_g_a(Nt,:); % 定点循环过程末状态气相表观流速，m/s
Va_l_a_3=Va_l_a(Nt,:); % 定点循环过程末状态环空液相表观流速，m/s
Va_s_3=Va_s(Nt,:); % 定点循环过程末状态固相表观流速，m/s
Vsr_3=Vsr(Nt,:); % 定点循环过程末状态固相滑移速度，m/s

%% 数据存储（定点循环）
ANS_Nt_3=Nt; % 时间节点数
ANS_Time_3=Time(1:Nt); % 定点循环时间，s
ANS_alpha_g_a_3=alpha_g_a(1:Nt,:); % 环空气相体积分数
ANS_alpha_l_a_3=alpha_l_a(1:Nt,:); % 环空液相体积分数
ANS_alpha_s_3=alpha_s(1:Nt,:); % 固相含量
ANS_P_ct_3=P_ct(1:Nt,:); % 管内压力（Pa）
ANS_delta_P_SWT_3=delta_P_SWT(1:Nt); % 冲砂工具压降（Pa）
ANS_P_a_3=P_a(1:Nt,:); % 环空压力（Pa）
ANS_T_a_3=T_a(1:Nt,:); % 环空温度（℃）
ANS_T_ct_3=T_ct(1:Nt,:); % 管内温度（℃）
ANS_P_coil_3=P_coil(1:Nt); % 盘管段出口压力（Pa）
ANS_P_pump_3=P_pump(1:Nt); % 泵压（Pa）
ANS_M_w_tem_3=M_w_tem(1:Nt); % 井口瞬时出砂量（kg）
ANS_M_w_tot_3=M_w_tot(1:Nt); % 井口累积出砂量（kg）
ANS_Va_s_3=Va_s(1:Nt,:); % 岩屑表观速度（m/s）
ANS_Va_f_a_3=Va_f_a(1:Nt,:); % 环空泡沫表观流速（m/s）
ANS_V_s_3=V_s(1:Nt,:); % 岩屑沉降速度（m/s）
ANS_V_f_a_3=V_f_a(1:Nt,:); % 环空泡沫返速（m/s）
ANS_alpha_g_ct_3=alpha_g_ct(1:Nt,:); % 管内气相体积分数
ANS_alpha_l_ct_3=alpha_l_ct(1:Nt,:); % 管内液相体积分数
ANS_gamma_g_ct_3=gamma_g_ct(1:Nt,:); % 管内泡沫质量
ANS_gamma_g_a_3=gamma_g_a(1:Nt,:); % 环空泡沫质量
ANS_ECD_a_3=ECD_a(1:Nt,:); % 环空ECD（kg/m^3）
ANS_L_coil_3=L_coil(1:Nt); % 连续管下深（m）



%% 连续管下入过程1（国际单位制）
%% 变量数值清零
dt=zeros(); % 时间步长，s
dx=zeros(); % 空间步长，m
Time=zeros(); % 连续管下入总时长，s
L_coil=zeros(); % 连续管下深，m
L_reel=zeros(); % 盘管段长度，m
D_t_i=zeros(); % 油管内径，m
D_ct_o=zeros(); % 连续油管外径，m
D_ct_i=zeros(); % 连续油管内径，m
M_s=zeros(); % 井底进砂量，kg/s
P_a=zeros(); % 环空压力，Pa
rho_g_a=zeros(); % 气相密度，kg/m^3
rho_l_a=zeros(); % 环空液相密度，kg/m^3
rho_f_a=zeros(); % 环空泡沫密度，kg/m^3
mu_g_a=zeros(); % 气相粘度，Pa*s
mu_l_a=zeros(); % 环空液体粘度，Pa*s
mu_f_a=zeros(); % 环空泡沫粘度，Pa*s
mu_s=zeros(); % 固相粘度，Pa*s
alpha_g_a=zeros(); % 气相含量
alpha_l_a=zeros(); % 液相含量
alpha_f_a=zeros(); % 泡沫含量
alpha_s=zeros(); % 固相含量
Va_g_a=zeros(); % 气相表观流速，m/s
Va_l_a=zeros(); % 环空液相表观流速，m/s
Va_f_a=zeros(); % 环空泡沫表观流速，m/s
Va_s=zeros(); % 固相表观流速，m/s
V_g_a=zeros(); % 气相流速，m/s
V_l_a=zeros(); % 环空液相流速，m/s
V_f_a=zeros(); % 环空泡沫流速，m/s
V_s=zeros(); % 固相流速，m/s
Vsr=zeros(); % 固相滑移速度，m/s
gamma_g_a=zeros(); % 泡沫质量
gamma_l_a=zeros(); % 液体滞留量
V_m=zeros(); % 环空混合物流速，m/s
rho_m=zeros(); % 环空混合物密度，kg/m^3
mu_m=zeros(); % 环空混合物粘度，Pa*s
Ff_a=zeros(); % 环空单位长度摩擦压降，Pa/m
f_a=zeros(); % 环空摩擦因子
Re_a=zeros(); % 环空雷诺数
flow_pattern_a=zeros(); % 环空流体流态

%% 计算空间步长dx及相应时间步长dt
V4=(Dsp_pen(1)-Lp)/(1*3600); % 连续油管下入速度，m/s（输入）
t_4=(Dsp_pen(1)-Lp)/V4; % 连续油管下入总时长，s

Nt=Nx_Dsp_pen(1)-Nx_Lp+1; % 时间节点数
nt=Nt-1; % 时间网格数

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % 空间步长，m
end

Time(1)=0; % 连续管下入总时长初值，s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Lp+t-1)/V4; % 每下入一个空间步长所需时间，s
    Time(t+1)=Time(t)+dt(t); % 连续管下入至第(t+1)个空间节点所经历的总时长，s
end

%% 计算不同时刻连续管下入深度L_coil及盘管段连续管长度L_reel
L=10000; % 连续油管总长，m（输入）
L_wg=8; % 井口到注入头顶部段连续管长度，m（输入）
L_goose=3; % 导向器段连续管长度，m（输入）
D_goose=2; % 导向器段半径，m（输入）
H_goose=10; % 导向器距地面高度，m（输入）
L_gr=20; % 导向器到滚筒段连续管长度，m（输入）
theta_gr=acosd(H_goose/L_gr); % 导向器到滚筒段连续管与铅垂线夹角，°
D_r_i=3; % 滚筒内径，m（输入）
D_r_o=5; % 滚筒外径，m（输入）
W_r=5; % 滚筒宽度，m（输入）
D_cable=0; % 电缆外径，m

L_coil(1)=Lp; % 初始时刻连续管下入深度，m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % 初始时刻盘管段连续管长度，m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)+dx(Nx_Lp+t-2); % 连续管下入深度，m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % 盘管段连续管长度，m
end
L_cable=L_coil; % 电缆下入长度，m

%% 数据输入及预处理
D_ct_o_0=0.04445; % 连续油管外径，m（输入）
L1=2000; % 近出口第一段连续管长度，m（输入）
D_ct_i_1=0.03709; % 近出口第一段连续管内径，m（输入）
L2=2000; % 近出口第二段连续管长度，m（输入）
D_ct_i_2=0.03653; % 近出口第二段连续管内径，m（输入）
L3=2000; % 近出口第三段连续管长度，m（输入）
D_ct_i_3=0.03555; % 近出口第三段连续管内径，m（输入）
L4=2000; % 近出口第四段连续管长度，m（输入）
D_ct_i_4=0.03489; % 近出口第四段连续管内径，m（输入）
L5=2000; % 近出口第五段连续管长度，m（输入）
D_ct_i_5=0.03409; % 近出口第五段连续管内径，m（输入）

L_t_1=4000; % 上部油管（或套管或裸眼）长度，m（输入）
D_t_i_1=0.09718; %0.068;%0.09718; % 上部油管（或套管或裸眼）内径，m（输入）
L_t_2=2200; % 下部油管（或套管或裸眼）长度，m（输入）
D_t_i_2=0.09718; %0.13970; % 下部油管（或套管或裸眼）内径，m（输入）

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % 油管（或套管或裸眼）内径，m
        else
            D_t_i(t,x)=D_t_i_2; % 油管（或套管或裸眼）内径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Dsp_pen(1)
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % 连续油管外径，m
            else
                D_ct_o(t,x)=0; % 连续油管外径，m
            end
        else
            D_ct_o(t,x)=0; % 连续油管外径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Dsp_pen(1)
            if Depth(x)<=L_coil(t) % 决定哪几个点可以有内径
                if (L_coil(t)-Depth(x))<=L1
                    D_ct_i(t,x)=D_ct_i_1; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % 连续油管内径，m
                end
            else
                D_ct_i(t,x)=0; % 连续油管内径，m
            end
        else
            D_ct_i(t,x)=0; % 连续油管内径，m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % 决定哪几个点可以有数值
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % 连续油管内截面积，m^2
        else
            A_ct(t,x)=0; % 连续油管内截面积，m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % 环空水力直径，m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % 环空截面积，m^2
    end
end

%% 数据输入及预处理
h_t=25.4*10^(-6); % 油管（或套管或裸眼）绝对粗糙度，m（输入）
h_ct=25.4*10^(-6); % 连续油管绝对粗糙度，m（输入）
h_a=(h_t+h_ct)/2; % 环空平均绝对粗糙度，m
epsilon_e=1*10^(-3); % 迭代求解误差限（输入）
epsilon_t=1*10^3; % 最大迭代次数（输入）
g=9.81; % 重力加速度，m/s^2（默认）

T_0=20; % 基液测试温度，℃（输入）
P_0=0.1*10^6; % 基液测试压力，Pa（输入）
rho_l_0=1150; % T_0、P_0下基液密度，kg/m^3（输入）
mu_l_0=0.03; % T_0、P_0下基液粘度，Pa*s（输入）
Qv_l_0=0.1/60; % 基液体积流量，m^3/s（输入）
Qm_l_0=Qv_l_0*rho_l_0; % 基液质量流量，kg/s
rho_g_0=0.655; % 注入气密度，kg/m^3
mu_g_0=RheologyG(T_0,P_0); % 气体粘度，Pa*s
Qv_g_0=4/60; % 注入气体积流量，m^3/s（输入）
Qm_g_0=Qv_g_0*rho_g_0;  % 注入气质量流量，kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % 泡沫质量流量，kg/s

D_s=1*10^(-3); % 砂砾直径，m（输入）
rho_s=2000; % 砂砾密度，kg/m^3（输入）
H_s=L_s_b-L_s_t; % 底部砂床高度，m
PHI=0.6; % 砂床充盈度（输入）
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % 井底砂砾总质量，kg

D_nozzle=4/1000; % 喷嘴直径，m（输入）
N_nozzle=3; % 喷嘴个数（输入）
C=0.95; % 喷嘴流量系数，取0.95（输入）

C0=1.2; % 漂移流方程系数（默认）

M_s(1)=0; % 井底进砂量，kg/s
for t=2:1:Nt
    M_s(t)=0; % 井底进砂量，kg/s
end

OutPressure=1*10^6; % 井口压力，Pa（输入）

%% 温度设置
T_i=20; % 泡沫注入温度，℃
T_g=0.02; % 地温梯度，℃/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % 连续管内泡沫温度（假设每个时刻都一样），℃
    end
end
T_a=T_ct; % 环空泡沫温度（假设每个时刻都一样），℃

%% 第1个时间节点（初始时刻）相关参数插值计算（环空）
for x=1:1:Nx
    P_a(1,x)=P_a_3(x); % 环空压力，Pa
    rho_g_a(1,x)=rho_g_a_3(x); % 气相密度，kg/m^3
    rho_l_a(1,x)=rho_l_a_3(x); % 环空液相密度，kg/m^3
    rho_f_a(1,x)=rho_f_a_3(x); % 环空泡沫密度，kg/m^3
    mu_g_a(1,x)=mu_g_a_3(x); % 气相粘度，Pa*s
    mu_l_a(1,x)=mu_l_a_3(x); % 环空液体粘度，Pa*s
    mu_f_a(1,x)=mu_f_a_3(x); % 环空泡沫粘度，Pa*s
    mu_s(1,x)=mu_s_3(x); % 固相粘度，Pa*s
    alpha_g_a(1,x)=alpha_g_a_3(x); % 气相含量
    alpha_l_a(1,x)=alpha_l_a_3(x); % 液相含量
    alpha_f_a(1,x)=alpha_f_a_3(x); % 泡沫含量
    alpha_s(1,x)=alpha_s_3(1,x); % 固相含量
    Va_g_a(1,x)=Va_g_a_3(x); % 气相表观流速，m/s
    Va_l_a(1,x)=Va_l_a_3(x); % 环空液相表观流速，m/s
    Va_f_a(1,x)=Va_f_a_3(x); % 环空泡沫表观流速，m/s
    Va_s(1,x)=Va_s_3(x); % 固相表观流速，m/s
    V_g_a(1,x)=V_g_a_3(x); % 气相流速，m/s
    V_l_a(1,x)=V_l_a_3(x); % 环空液相流速，m/s
    V_f_a(1,x)=V_f_a_3(x); % 环空泡沫流速，m/s
    V_s(1,x)=V_s_3(x); % 固相流速，m/s
    Vsr(1,x)=Vsr_3(x); % 固相滑移速度，m/s
    gamma_g_a(1,x)=gamma_g_a_3(x); % 泡沫质量
    gamma_l_a(1,x)=gamma_l_a_3(1,x); % 液体滞留量
    V_m(1,x)=V_m_3(x); % 环空混合物流速，m/s
    rho_m(1,x)=rho_m_3(x); % 环空混合物密度，kg/m^3
    mu_m(1,x)=mu_m_3(x); % 环空混合物粘度，Pa*s
    Ff_a(1,x)=Ff_a_3(x); % 环空单位长度摩擦压降，Pa/m
    f_a(1,x)=f_a_3(x); % 环空摩擦因子
    Re_a(1,x)=Re_a_3(x); % 环空雷诺数
    flow_pattern_a(1,x)=flow_pattern_a_3(x); % 环空流体流态
end

%% 第2～Nt个时间节点相关参数计算（环空）
for t=2:1:Nt
    P_a(t,Nx_Lp+t-1)=P_a(t-1,Nx_Lp+t-1);  % 环空管底压力假设值，Pa
    
    err_OutPressure=1; % 出口压力相对误差
    COUNT_OutPressure=0; % 出口压力迭代次数初值
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % 出口压力迭代次数
        
        % 管底（第Nx_Lp+t-1个空间节点）处相关参数计算
        rho_g_a(t,Nx_Lp+t-1)=DensityG(T_a(t,Nx_Lp+t-1),P_a(t,Nx_Lp+t-1)); % 环空气体密度，kg/m^3
        rho_l_a(t,Nx_Lp+t-1)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Lp+t-1),P_a(t,Nx_Lp+t-1)); % 环空基液密度，kg/m^3
        mu_g_a(t,Nx_Lp+t-1)=RheologyG(T_a(t,Nx_Lp+t-1),P_a(t,Nx_Lp+t-1)); % 环空气体粘度，Pa*s
        mu_l_a(t,Nx_Lp+t-1)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Lp+t-1),P_a(t,Nx_Lp+t-1)); % 环空基液粘度，Pa*s
        gamma_g_a(t,Nx_Lp+t-1)=(Qm_g_0/rho_g_a(t,Nx_Lp+t-1))/(Qm_g_0/rho_g_a(t,Nx_Lp+t-1)+Qm_l_0/rho_l_a(t,Nx_Lp+t-1)); % 泡沫质量
        gamma_l_a(t,Nx_Lp+t-1)=1-gamma_g_a(t,Nx_Lp+t-1); % 液体滞留量
        alpha_g_a(t,Nx_Lp+t-1)=alpha_f_a(t-1,Nx_Lp+t-1)*gamma_g_a(t,Nx_Lp+t-1); % 环空气体含量
        alpha_l_a(t,Nx_Lp+t-1)=alpha_f_a(t-1,Nx_Lp+t-1)*gamma_l_a(t,Nx_Lp+t-1); % 环空基液含量
        rho_f_a(t,Nx_Lp+t-1)=rho_g_a(t,Nx_Lp+t-1)*gamma_g_a(t,Nx_Lp+t-1)+rho_l_a(t,Nx_Lp+t-1)*gamma_l_a(t,Nx_Lp+t-1); % 环空泡沫密度，kg/m^3
        mu_f_a(t,Nx_Lp+t-1)=mu_g_a(t,Nx_Lp+t-1)*gamma_g_a(t,Nx_Lp+t-1)+mu_l_a(t,Nx_Lp+t-1)*gamma_l_a(t,Nx_Lp+t-1); % 环空泡沫粘度，Pa*s
        mu_s(t,Nx_Lp+t-1)=mu_f_a(t,Nx_Lp+t-1); % 固相粘度，Pa*s
        Va_f_a(t,Nx_Lp+t-1)=Qm_f_0/(A_a(t,Nx_Lp+t-1)*rho_f_a(t,Nx_Lp+t-1)); % 环空泡沫表观流速，m/s
        Va_g_a(t,Nx_Lp+t-1)=Va_f_a(t,Nx_Lp+t-1); % 环空气体表观流速，m/s
        Va_l_a(t,Nx_Lp+t-1)=Va_f_a(t,Nx_Lp+t-1); % 环空基液表观流速，m/s
        Va_s(t,Nx_Lp+t-1)=M_s(t)/(A_a(t,Nx_Lp+t-1)*rho_s); % 岩屑表观速度，m/s
        Vsr(t,Nx_Lp+t-1)=12*(mu_f_a(t,Nx_Lp+t-1)/(rho_f_a(t,Nx_Lp+t-1)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Lp+t-1))/rho_f_a(t,Nx_Lp+t-1))*((rho_f_a(t,Nx_Lp+t-1)*D_s/mu_f_a(t,Nx_Lp+t-1))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
        alpha_s(t,Nx_Lp+t-1)=Va_s(t,Nx_Lp+t-1)/(C0*(Va_s(t,Nx_Lp+t-1)+Va_f_a(t,Nx_Lp+t-1))-Vsr(t,Nx_Lp+t-1));  % 固相体积分数
        V_s(t,Nx_Lp+t-1)=0;%Va_s(t,Nx_Lp+t-1)/alpha_s(t,Nx_Lp+t-1); % 岩屑速度，m/s
        
        alpha_f_a(t,Nx_Lp+t-1)=1-alpha_s(t,Nx_Lp+t-1); % 环空泡沫含量
        alpha_g_a(t,Nx_Lp+t-1)=alpha_f_a(t,Nx_Lp+t-1)*gamma_g_a(t,Nx_Lp+t-1); % 环空气体含量
        alpha_l_a(t,Nx_Lp+t-1)=alpha_f_a(t,Nx_Lp+t-1)*gamma_l_a(t,Nx_Lp+t-1); % 环空基液含量
        V_f_a(t,Nx_Lp+t-1)=Va_f_a(t,Nx_Lp+t-1)/alpha_f_a(t,Nx_Lp+t-1); % 环空泡沫流速，m/s
        V_g_a(t,Nx_Lp+t-1)=V_f_a(t,Nx_Lp+t-1); % 环空气体流速，m/s
        V_l_a(t,Nx_Lp+t-1)=V_f_a(t,Nx_Lp+t-1); % 环空基液流速，m/s
        
        V_m(t,Nx_Lp+t-1)=Va_s(t,Nx_Lp+t-1)+Va_f_a(t,Nx_Lp+t-1); % 环空混合物速度，m/s
        rho_m(t,Nx_Lp+t-1)=alpha_s(t,Nx_Lp+t-1)*rho_s+alpha_f_a(t,Nx_Lp+t-1)*rho_f_a(t,Nx_Lp+t-1); % 环空混合物密度，kg/m^3
        mu_m(t,Nx_Lp+t-1)=alpha_s(t,Nx_Lp+t-1)*mu_s(t,Nx_Lp+t-1)+alpha_f_a(t,Nx_Lp+t-1)*mu_f_a(t,Nx_Lp+t-1); % 环空混合物粘度，Pa*s
        [Ff_a(t,Nx_Lp+t-1),f_a(t,Nx_Lp+t-1),Re_a(t,Nx_Lp+t-1),flow_pattern_a(t,Nx_Lp+t-1)]=Friction_annulus(rho_m(t,Nx_Lp+t-1),V_m(t,Nx_Lp+t-1),mu_m(t,Nx_Lp+t-1),D_h(t,Nx_Lp+t-1),h_a,rho_f_a(t,Nx_Lp+t-1),V_f_a(t,Nx_Lp+t-1)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
        
        % 第Nx_Lp+t-2～1个空间节点处相关参数计算
        for x=Nx_Lp+t-2:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % 环空压力假设值，Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液粘度，Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % 环空基液含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            
            err_NodePressure=1; % 环空压力相对误差
            COUNT_NodePressure=0; % 环空压力迭代次数初值
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % 环空压力迭代次数
                
                % 割线法求解固相体积分数
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % 固相体积分数假设值1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % 固相体积分数假设值2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % 固相体积分数绝对误差
                COUNT_NodeEg=0; % 迭代次数
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % 固相体积分数为alpha_s_ass1时（假设值）
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % 固相速度，m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % 固相表观流速，m/s
                    alpha_f_ass1=1-alpha_s_ass1; % 泡沫体积分数
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % 泡沫速度，m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % 泡沫表观流速，m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % 环空混合物速度，m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % 固相体积分数计算值
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 固相体积分数为alpha_s_ass2时（假设值）
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % 固相质量守恒方程离散后公式中间值计算
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % 固相速度，m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % 固相表观流速，m/s
                    alpha_f_ass2=1-alpha_s_ass2; % 泡沫体积分数
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % 泡沫速度，m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % 泡沫表观流速，m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % 环空混合物速度，m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % 固相体积分数计算值
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 割线法计算固相体积分数
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % 新的固相体积分数假设值alpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % 固相体积分数绝对误差
                    alpha_s_ass1=alpha_s_ass2; % 新的固相体积分数假设值1
                    alpha_s_ass2=alpha_s_ass3; % 新的固相体积分数假设值2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % 将迭代求解得到的真实固相体积分数值赋给alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % 当固相体积分数小于一定值时，认为固相体积数为0，用于防止出现后续的计算误差
                    alpha_s(t,x)=0; % 固相体积分数
                    V_s(t,x)=0; % 固相速度，m/s
                    Va_s(t,x)=0; % 固相表观速度，m/s
                    mu_s(t,x)=0; % 固相粘度，Pa*s
                    Vsr(t,x)=0; % 岩屑沉降末速，m/s
                    
                    alpha_f_a(t,x)=1; % 泡沫体积分数                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                else
                    V_s(t,x)=V_s_ass2; % 固相速度，m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % 固相表观流速，m/s
                    Vsr(t,x)=Vsr_ass2; % 岩屑沉降末速，m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % 泡沫密度，kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % 泡沫体积分数
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    V_f_a(t,x)=V_f_ass2; % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % 环空混合物速度，m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % 环空压力计算值，Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % 环空压力相对误差
                P_a(t,x)=P_new; % 新的环空压力假设值，Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % 出口压力相对误差        
        if (P_a(t,1)-OutPressure)>0 % 根据出口压力误差的正负，对环空管底压力假设值进行调节
            P_a(t,Nx_Lp+t-1)=P_a(t,Nx_Lp+t-1)-(P_a(t,1)-OutPressure)/2; % 新的环空管底压力假设值，Pa
        else
            P_a(t,Nx_Lp+t-1)=P_a(t,Nx_Lp+t-1)-(P_a(t,1)-OutPressure)/2*0.3; % 新的环空管底压力假设值，Pa
        end
    end
    
    % 第Nx_Lp+t～Nx_Dsp_pen(1)个空间节点处相关参数计算
    for x=Nx_Lp+t:1:Nx_Dsp_pen(1)
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=0; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(1)+1～Nx个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)+1:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1-PHI; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=PHI; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
end

%% 第1～Nt时间节点钻头压降计算
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Lp+t-1); % 井底泡沫体积流量，m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Lp+t-1)); % 计算冲砂工具压降（Pa）及射流喷嘴流速V_nozzle（m/s）
end

%% 第1～Nt时间节点相关参数计算（连续管内）
for t=1:1:Nt
    % 管底（第Nx_Lp+t-1个空间节点）处相关参数计算
    P_ct(t,Nx_Lp+t-1)=P_a(t,Nx_Lp+t-1)+delta_P_SWT(t); % 管内压力，Pa
    rho_g_ct(t,Nx_Lp+t-1)=DensityG(T_ct(t,Nx_Lp+t-1),P_ct(t,Nx_Lp+t-1)); % 管内气体密度，kg/m^3
    rho_l_ct(t,Nx_Lp+t-1)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Lp+t-1),P_ct(t,Nx_Lp+t-1)); % 管内基液密度，kg/m^3
    mu_g_ct(t,Nx_Lp+t-1)=RheologyG(T_ct(t,Nx_Lp+t-1),P_ct(t,Nx_Lp+t-1)); % 管内气体粘度，Pa*s
    mu_l_ct(t,Nx_Lp+t-1)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Lp+t-1),P_ct(t,Nx_Lp+t-1)); % 管内基液粘度，Pa*s
    alpha_f_ct(t,Nx_Lp+t-1)=1; % 管内泡沫含量
    gamma_g_ct(t,Nx_Lp+t-1)=(Qm_g_0/rho_g_ct(t,Nx_Lp+t-1))/(Qm_g_0/rho_g_ct(t,Nx_Lp+t-1)+Qm_l_0/rho_l_ct(t,Nx_Lp+t-1)); % 管内泡沫质量
    gamma_l_ct(t,Nx_Lp+t-1)=1-gamma_g_ct(t,Nx_Lp+t-1); % 管内液体滞留量
    alpha_g_ct(t,Nx_Lp+t-1)=alpha_f_ct(t,Nx_Lp+t-1)*gamma_g_ct(t,Nx_Lp+t-1); % 管内气体含量
    alpha_l_ct(t,Nx_Lp+t-1)=alpha_f_ct(t,Nx_Lp+t-1)*gamma_l_ct(t,Nx_Lp+t-1); % 管内基液含量
    rho_f_ct(t,Nx_Lp+t-1)=rho_g_ct(t,Nx_Lp+t-1)*gamma_g_ct(t,Nx_Lp+t-1)+rho_l_ct(t,Nx_Lp+t-1)*gamma_l_ct(t,Nx_Lp+t-1); % 管内泡沫密度，kg/m^3
    mu_f_ct(t,Nx_Lp+t-1)=mu_g_ct(t,Nx_Lp+t-1)*gamma_g_ct(t,Nx_Lp+t-1)+mu_l_ct(t,Nx_Lp+t-1)*gamma_l_ct(t,Nx_Lp+t-1); % 管内泡沫粘度，Pa*s
    V_f_ct(t,Nx_Lp+t-1)=Qm_f_0/(rho_f_ct(t,Nx_Lp+t-1)*A_ct(t,Nx_Lp+t-1)); % 管内气体流速，m/s
    [Ff_ct(t,Nx_Lp+t-1),f_ct(t,Nx_Lp+t-1),Re_ct(t,Nx_Lp+t-1),flow_pattern_ct(t,Nx_Lp+t-1)]=Friction_coiledtubing(rho_f_ct(t,Nx_Lp+t-1),V_f_ct(t,Nx_Lp+t-1),mu_f_ct(t,Nx_Lp+t-1),D_ct_i(t,Nx_Lp+t-1),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    % 第Nx_Lp+t-2～1个空间节点处相关参数计算
    for x=Nx_Lp+t-2:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % 管内压力假设值，Pa
        
        err_DriPipePressure=1; % 管内压力相对误差
        COUNT_DriPipePressure=0; % 管内压力迭代次数初值
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % 管内压力迭代次数
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体密度，kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液密度，kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体粘度，Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液粘度，Pa*s
            alpha_f_ct(t,x)=1; % 管内泡沫含量
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % 管内泡沫质量
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % 管内气体含量
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % 管内基液含量
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫密度，kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫粘度，Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % 管内泡沫流速，m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % 管内压力，Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % 计算管内压力假设值与计算值之间的相对误差
            P_ct_ass(t,x)=P_ct(t,x); % 新的管内压力假设值，Pa
        end
    end
    
    % 第Nx_Lp+t～Nx个空间节点处相关参数计算
    for x=Nx_Lp+t:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % 管内气体含量
        alpha_l_ct(t,x)=alpha_l_a(t,x); % 管内基液含量
        gamma_g_ct(t,x)=gamma_g_a(t,x); % 管内泡沫质量
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
        rho_f_ct(t,x)=rho_m(t,x); % 管内气体密度，kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % 管内气体粘度，Pa*s
        V_f_ct(t,t)=V_m(t,x); % 管内气体流速，sm/s
        Ff_ct(t,x)=Ff_a(t,x); % 管内流体单位长度摩擦压降（Pa/m）
        f_ct(t,x)=f_a(t,x); % 管内流体范宁摩擦因子
        Re_ct(t,x)=Re_a(t,x); % 管内流体雷诺数、流体流态
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % 管内流体流态
        P_ct(t,x)=P_a(t,x); % 管内压力，Pa
    end
end

%% 第1～Nt时间节点盘管段出口压力及泵压计算（连续管内）
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % 地面管内流体流速，m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % 地面管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % 盘管段流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % 导向器内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % 井口到注入头顶部段摩擦压降，Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % 导向器段摩擦压降，Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % 导向器到滚筒段摩擦压降，Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % 盘管段摩擦压降，Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % 盘管段出口压力，Pa
    P_pumb(t)=P_coil(t)+delta_P_coil(t); % 泵压，Pa
end

%% 计算井口累积出砂量
M_w_tem(1)=0; % 井口瞬时出砂量初值，kg
M_w_tot(1)=0; % 井口累积出砂量初值，kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口瞬时出砂量，kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口累积出砂量，kg
end

%% 计算末态井筒最大砂浓度
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(2)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% 计算环空泡沫返速平均值V_l_mre
V_f_mre=0; % 环空返速平均值，m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % 环空返速平均值，m/s
    end
end

%% 计算固相沉降末速平均值Vsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % 从Vsr矩阵中提取值非零的沉降末速，m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % 固相沉降末速平均值，m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % 固相沉降末速平均值，m/s
end

%% 判断是否有效冲砂（环空返速平均值大于2倍砂砾沉降末速平均值）、是否完成冲砂
fprintf("连续管下入过程1：\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“有效冲砂”
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % 命令行窗口说明当前条件下“完成冲砂”
    else
        fprintf("Sand Cleanout UnFinished!\n"); % 命令行窗口说明当前条件下“未完成冲砂”
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“无效冲砂”
end

%% 井筒ECD计算
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % 环空ECD，kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % 管内ECD，kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % 环空ECD，kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % 管内ECD，kg/m^3
end

%% 连续管下入过程末状态数据
alpha_f_a_4=alpha_f_a(Nt,:); % 连续管下入过程末状态泡沫含量
alpha_g_a_4=alpha_g_a(Nt,:); % 连续管下入过程末状态气相含量
alpha_l_a_4=alpha_l_a(Nt,:); % 连续管下入过程末状态液相含量
alpha_s_4=alpha_s(Nt,:); % 连续管下入过程末状态固相含量
f_a_4=f_a(Nt,:); % 连续管下入过程末状态环空摩擦因子
Ff_a_4=Ff_a(Nt,:); % 连续管下入过程末状态环空单位长度摩擦压降，Pa/m
flow_pattern_a_4=flow_pattern_a(Nt,:); % 连续管下入过程末状态环空流体流态
gamma_g_a_4=gamma_g_a(Nt,:); % 连续管下入过程末状态泡沫质量
gamma_l_a_4=gamma_l_a(Nt,:); % 连续管下入过程末状态液体滞留量
mu_f_a_4=mu_f_a(Nt,:); % 连续管下入过程末状态环空泡沫粘度，Pa*s
mu_g_a_4=mu_g_a(Nt,:); % 连续管下入过程末状态气相粘度，Pa*s
mu_l_a_4=mu_l_a(Nt,:); % 连续管下入过程末状态环空液体粘度，Pa*s
mu_m_4=mu_m(Nt,:); % 连续管下入过程末状态环空混合物粘度，Pa*s
mu_s_4=mu_s(Nt,:); % 连续管下入过程末状态固相粘度，Pa*s
P_a_4=P_a(Nt,:); % 连续管下入过程末状态环空压力，Pa
Re_a_4=Re_a(Nt,:); % 连续管下入过程末状态环空雷诺数
rho_f_a_4=rho_f_a(Nt,:); % 连续管下入过程末状态环空泡沫密度，kg/m^3
rho_g_a_4=rho_g_a(Nt,:); % 连续管下入过程末状态环空气相密度，kg/m^3
rho_l_a_4=rho_l_a(Nt,:); % 连续管下入过程末状态环空液相密度，kg/m^3
rho_m_4=rho_m(Nt,:); % 连续管下入过程末状态环空混合物密度，kg/m^3
V_f_a_4=V_f_a(Nt,:); % 连续管下入过程末状态环空泡沫流速，m/s
V_g_a_4=V_g_a(Nt,:); % 连续管下入过程末状态气相流速，m/s
V_l_a_4=V_l_a(Nt,:); % 连续管下入过程末状态环空液相流速，m/s
V_m_4=V_m(Nt,:); % 连续管下入过程末状态环空混合物流速，m/s
V_s_4=V_s(Nt,:); % 连续管下入过程末状态固相流速，m/s
Va_f_a_4=Va_f_a(Nt,:); % 连续管下入过程末状态环空泡沫表观流速，m/s
Va_g_a_4=Va_g_a(Nt,:); % 连续管下入过程末状态气相表观流速，m/s
Va_l_a_4=Va_l_a(Nt,:); % 连续管下入过程末状态环空液相表观流速，m/s
Va_s_4=Va_s(Nt,:); % 连续管下入过程末状态固相表观流速，m/s
Vsr_4=Vsr(Nt,:); % 连续管下入过程末状态固相滑移速度，m/s

%% 数据存储（连续管下入）
ANS_Nt_4=Nt; % 时间节点数
ANS_Time_4=Time(1:Nt); % 定点循环时间，s
ANS_alpha_g_a_4=alpha_g_a(1:Nt,:); % 环空气相体积分数
ANS_alpha_l_a_4=alpha_l_a(1:Nt,:); % 环空液相体积分数
ANS_alpha_s_4=alpha_s(1:Nt,:); % 固相含量
ANS_P_ct_4=P_ct(1:Nt,:); % 管内压力（Pa）
ANS_delta_P_SWT_4=delta_P_SWT(1:Nt); % 冲砂工具压降（Pa）
ANS_P_a_4=P_a(1:Nt,:); % 环空压力（Pa）
ANS_T_a_4=T_a(1:Nt,:); % 环空温度（℃）
ANS_T_ct_4=T_ct(1:Nt,:); % 管内温度（℃）
ANS_P_coil_4=P_coil(1:Nt); % 盘管段出口压力（Pa）
ANS_P_pump_4=P_pump(1:Nt); % 泵压（Pa）
ANS_M_w_tem_4=M_w_tem(1:Nt); % 井口瞬时出砂量（kg）
ANS_M_w_tot_4=M_w_tot(1:Nt); % 井口累积出砂量（kg）
ANS_Va_s_4=Va_s(1:Nt,:); % 岩屑表观速度（m/s）
ANS_Va_f_a_4=Va_f_a(1:Nt,:); % 环空泡沫表观流速（m/s）
ANS_V_s_4=V_s(1:Nt,:); % 岩屑沉降速度（m/s）
ANS_V_f_a_4=V_f_a(1:Nt,:); % 环空泡沫返速（m/s）
ANS_alpha_g_ct_4=alpha_g_ct(1:Nt,:); % 管内气相体积分数
ANS_alpha_l_ct_4=alpha_l_ct(1:Nt,:); % 管内液相体积分数
ANS_gamma_g_ct_4=gamma_g_ct(1:Nt,:); % 管内泡沫质量
ANS_gamma_g_a_4=gamma_g_a(1:Nt,:); % 环空泡沫质量
ANS_ECD_a_4=ECD_a(1:Nt,:); % 环空ECD（kg/m^3）
ANS_L_coil_4=L_coil(1:Nt); % 连续管下深（m）



%% 冲洗钻进过程2（国际单位制）
%% 变量数值清零
dt=zeros(); % 时间步长，s
dx=zeros(); % 空间步长，m
Time=zeros(); % 连续管下入总时长，s
L_coil=zeros(); % 连续管下深，m
L_reel=zeros(); % 盘管段长度，m
D_t_i=zeros(); % 油管内径，m
D_ct_o=zeros(); % 连续油管外径，m
D_ct_i=zeros(); % 连续油管内径，m
M_s=zeros(); % 井底进砂量，kg/s
P_a=zeros(); % 环空压力，Pa
rho_g_a=zeros(); % 气相密度，kg/m^3
rho_l_a=zeros(); % 环空液相密度，kg/m^3
rho_f_a=zeros(); % 环空泡沫密度，kg/m^3
mu_g_a=zeros(); % 气相粘度，Pa*s
mu_l_a=zeros(); % 环空液体粘度，Pa*s
mu_f_a=zeros(); % 环空泡沫粘度，Pa*s
mu_s=zeros(); % 固相粘度，Pa*s
alpha_g_a=zeros(); % 气相含量
alpha_l_a=zeros(); % 液相含量
alpha_f_a=zeros(); % 泡沫含量
alpha_s=zeros(); % 固相含量
Va_g_a=zeros(); % 气相表观流速，m/s
Va_l_a=zeros(); % 环空液相表观流速，m/s
Va_f_a=zeros(); % 环空泡沫表观流速，m/s
Va_s=zeros(); % 固相表观流速，m/s
V_g_a=zeros(); % 气相流速，m/s
V_l_a=zeros(); % 环空液相流速，m/s
V_f_a=zeros(); % 环空泡沫流速，m/s
V_s=zeros(); % 固相流速，m/s
Vsr=zeros(); % 固相滑移速度，m/s
gamma_g_a=zeros(); % 泡沫质量
gamma_l_a=zeros(); % 液体滞留量
V_m=zeros(); % 环空混合物流速，m/s
rho_m=zeros(); % 环空混合物密度，kg/m^3
mu_m=zeros(); % 环空混合物粘度，Pa*s
Ff_a=zeros(); % 环空单位长度摩擦压降，Pa/m
f_a=zeros(); % 环空摩擦因子
Re_a=zeros(); % 环空雷诺数
flow_pattern_a=zeros(); % 环空流体流态

%% 计算空间步长dx及相应时间步长dt
V5=(Dsp_pen(2)-Dsp_pen(1))/(1*3600); % 冲洗钻进速度，m/s（输入）
t_5=(Dsp_pen(2)-Dsp_pen(1))/V5; % 冲洗钻进总时长，s

Nt=Nx_Dsp_pen(2)-Nx_Dsp_pen(1)+1; % 时间节点数
nt=Nt-1; % 时间网格数

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % 空间步长，m
end

Time(1)=0; % 冲洗钻进总时长初值，s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Dsp_pen(1)+t-1)/V5; % 每冲洗钻进一个空间步长所需时间，s
    Time(t+1)=Time(t)+dt(t); % 冲洗钻进至第(t+1)个空间节点所经历的总时长，s
end

%% 计算不同时刻连续管冲洗钻进深度L_coil及盘管段连续管长度L_reel
L=10000; % 连续油管总长，m（输入）
L_wg=8; % 井口到注入头顶部段连续管长度，m（输入）
L_goose=3; % 导向器段连续管长度，m（输入）
D_goose=2; % 导向器段半径，m（输入）
H_goose=10; % 导向器距地面高度，m（输入）
L_gr=20; % 导向器到滚筒段连续管长度，m（输入）
theta_gr=acosd(H_goose/L_gr); % 导向器到滚筒段连续管与铅垂线夹角，°
D_r_i=3; % 滚筒内径，m（输入）
D_r_o=5; % 滚筒外径，m（输入）
W_r=5; % 滚筒宽度，m（输入）
D_cable=0.005; % 电缆外径，m

L_coil(1)=Dsp_pen(1); % 初始时刻连续管底部深度，m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % 初始时刻盘管段连续管长度，m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)+dx(Nx_Dsp_pen(1)+t-2); % 冲洗钻进深度，m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % 盘管段连续管长度，m
end
L_cable=L_coil; % 电缆下入长度，m

%% 数据输入及预处理
D_ct_o_0=0.04445; % 连续油管外径，m（输入）
L1=2000; % 近出口第一段连续管长度，m（输入）
D_ct_i_1=0.03709; % 近出口第一段连续管内径，m（输入）
L2=2000; % 近出口第二段连续管长度，m（输入）
D_ct_i_2=0.03653; % 近出口第二段连续管内径，m（输入）
L3=2000; % 近出口第三段连续管长度，m（输入）
D_ct_i_3=0.03555; % 近出口第三段连续管内径，m（输入）
L4=2000; % 近出口第四段连续管长度，m（输入）
D_ct_i_4=0.03489; % 近出口第四段连续管内径，m（输入）
L5=2000; % 近出口第五段连续管长度，m（输入）
D_ct_i_5=0.03409; % 近出口第五段连续管内径，m（输入）

L_t_1=4000; % 上部油管（或套管或裸眼）长度，m（输入）
D_t_i_1=0.09718; %0.068;% 上部油管（或套管或裸眼）内径，m（输入）
L_t_2=2200; % 下部油管（或套管或裸眼）长度，m（输入）
D_t_i_2=0.09718; %0.13970; % 下部油管（或套管或裸眼）内径，m（输入）

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % 油管（或套管或裸眼）内径，m
        else
            D_t_i(t,x)=D_t_i_2; % 油管（或套管或裸眼）内径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(2)
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % 连续油管外径，m
            else
                D_ct_o(t,x)=0; % 连续油管外径，m
            end
        else
            D_ct_o(t,x)=0; % 连续油管外径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(2)
            if Depth(x)<=L_coil(t) % 决定哪几个点可以有内径
                if (L_coil(t)-Depth(x))<=L1 %&& (L_coil(t)-Depth(x))>=0
                    D_ct_i(t,x)=D_ct_i_1; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % 连续油管内径，m
                end
            else
                D_ct_i(t,x)=0; % 连续油管内径，m
            end
        else
            D_ct_i(t,x)=0; % 连续油管内径，m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % 决定哪几个点可以有数值
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % 连续油管内截面积，m^2
        else
            A_ct(t,x)=0; % 连续油管内截面积，m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % 环空水力直径，m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % 环空截面积，m^2
    end
end

%% 数据输入及预处理
h_t=25.4*10^(-6); % 油管（或套管或裸眼）绝对粗糙度，m（输入）
h_ct=25.4*10^(-6); % 连续油管绝对粗糙度，m（输入）
h_a=(h_t+h_ct)/2; % 环空平均绝对粗糙度，m
epsilon_e=1*10^(-3); % 迭代求解误差限（输入）
epsilon_t=1*10^3; % 最大迭代次数（输入）
g=9.81; % 重力加速度，m/s^2（默认）

T_0=20; % 基液测试温度，℃（输入）
P_0=0.1*10^6; % 基液测试压力，Pa（输入）
rho_l_0=1150; % T_0、P_0下基液密度，kg/m^3（输入）
mu_l_0=0.03; % T_0、P_0下基液粘度，Pa*s（输入）
Qv_l_0=0.1/60; % 基液体积流量，m^3/s（输入）
Qm_l_0=Qv_l_0*rho_l_0; % 基液质量流量，kg/s
rho_g_0=0.655; % 注入气密度，kg/m^3
mu_g_0=RheologyG(T_0,P_0); % 气体粘度，Pa*s
Qv_g_0=4/60; % 注入气体积流量，m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % 注入气质量流量，kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % 泡沫质量流量，kg/s

D_s=1*10^(-3); % 砂砾直径，m（输入）
rho_s=2000; % 砂砾密度，kg/m^3（输入）
H_s=L_s_b-L_s_t; % 底部砂床高度，m
PHI=0.6; % 砂床充盈度（输入）
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % 井底砂砾总质量，kg

D_nozzle=4/1000; % 喷嘴直径，m（输入）
N_nozzle=3; % 喷嘴个数（输入）
C=0.95; % 喷嘴流量系数，取0.95（输入）

C0=1.2; % 漂移流方程系数（默认）

M_s(1)=0; % 井底进砂量，kg/s
for t=2:1:Nt
    M_s(t)=PHI*rho_s*V5*0.25*pi*D_t_i(t,Nx_Dsp_pen(1)+t-1)^2; % 井底进砂量，kg/s
end

OutPressure=1*10^6; % 井口压力，Pa（输入）

%% 温度设置
T_i=20; % 泡沫注入温度，℃
T_g=0.02; % 地温梯度，℃/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % 连续管内泡沫温度（假设每个时刻都一样），℃
    end
end
T_a=T_ct; % 环空泡沫温度（假设每个时刻都一样），℃

%% 第1个时间节点（初始时刻）相关参数插值计算（环空）
for x=1:1:Nx
    P_a(1,x)=P_a_4(x); % 环空压力，Pa
    rho_g_a(1,x)=rho_g_a_4(x); % 气相密度，kg/m^3
    rho_l_a(1,x)=rho_l_a_4(x); % 环空液相密度，kg/m^3
    rho_f_a(1,x)=rho_f_a_4(x); % 环空泡沫密度，kg/m^3
    mu_g_a(1,x)=mu_g_a_4(x); % 气相粘度，Pa*s
    mu_l_a(1,x)=mu_l_a_4(x); % 环空液体粘度，Pa*s
    mu_f_a(1,x)=mu_f_a_4(x); % 环空泡沫粘度，Pa*s
    mu_s(1,x)=mu_s_4(x); % 固相粘度，Pa*s
    alpha_g_a(1,x)=alpha_g_a_4(x); % 气相含量
    alpha_l_a(1,x)=alpha_l_a_4(x); % 液相含量
    alpha_f_a(1,x)=alpha_f_a_4(x); % 泡沫含量
    alpha_s(1,x)=alpha_s_4(1,x); % 固相含量
    Va_g_a(1,x)=Va_g_a_4(x); % 气相表观流速，m/s
    Va_l_a(1,x)=Va_l_a_4(x); % 环空液相表观流速，m/s
    Va_f_a(1,x)=Va_f_a_4(x); % 环空泡沫表观流速，m/s
    Va_s(1,x)=Va_s_4(x); % 固相表观流速，m/s
    V_g_a(1,x)=V_g_a_4(x); % 气相流速，m/s
    V_l_a(1,x)=V_l_a_4(x); % 环空液相流速，m/s
    V_f_a(1,x)=V_f_a_4(x); % 环空泡沫流速，m/s
    V_s(1,x)=V_s_4(x); % 固相流速，m/s
    Vsr(1,x)=Vsr_4(x); % 固相滑移速度，m/s
    gamma_g_a(1,x)=gamma_g_a_4(x); % 泡沫质量
    gamma_l_a(1,x)=gamma_l_a_4(1,x); % 液体滞留量
    V_m(1,x)=V_m_4(x); % 环空混合物流速，m/s
    rho_m(1,x)=rho_m_4(x); % 环空混合物密度，kg/m^3
    mu_m(1,x)=mu_m_4(x); % 环空混合物粘度，Pa*s
    Ff_a(1,x)=Ff_a_4(x); % 环空单位长度摩擦压降，Pa/m
    f_a(1,x)=f_a_4(x); % 环空摩擦因子
    Re_a(1,x)=Re_a_4(x); % 环空雷诺数
    flow_pattern_a(1,x)=flow_pattern_a_4(x); % 环空流体流态
end

%% 第2～Nt个时间节点相关参数计算（环空）
for t=2:1:Nt
    P_a(t,Nx_Dsp_pen(1)+t-1)=P_a(t-1,Nx_Dsp_pen(1)+t-1);  % 环空管底压力假设值，Pa
    
    err_OutPressure=1; % 出口压力相对误差
    COUNT_OutPressure=0; % 出口压力迭代次数初值
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % 出口压力迭代次数
        
        % 管底（第Nx_Dsp_pen(1)+t-1个空间节点）处相关参数计算
        rho_g_a(t,Nx_Dsp_pen(1)+t-1)=DensityG(T_a(t,Nx_Dsp_pen(1)+t-1),P_a(t,Nx_Dsp_pen(1)+t-1)); % 环空气体密度，kg/m^3
        rho_l_a(t,Nx_Dsp_pen(1)+t-1)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)+t-1),P_a(t,Nx_Dsp_pen(1)+t-1)); % 环空基液密度，kg/m^3
        mu_g_a(t,Nx_Dsp_pen(1)+t-1)=RheologyG(T_a(t,Nx_Dsp_pen(1)+t-1),P_a(t,Nx_Dsp_pen(1)+t-1)); % 环空气体粘度，Pa*s
        mu_l_a(t,Nx_Dsp_pen(1)+t-1)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)+t-1),P_a(t,Nx_Dsp_pen(1)+t-1)); % 环空基液粘度，Pa*s
        gamma_g_a(t,Nx_Dsp_pen(1)+t-1)=(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)+t-1))/(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)+t-1)+Qm_l_0/rho_l_a(t,Nx_Dsp_pen(1)+t-1)); % 泡沫质量
        gamma_l_a(t,Nx_Dsp_pen(1)+t-1)=1-gamma_g_a(t,Nx_Dsp_pen(1)+t-1); % 液体滞留量
        alpha_g_a(t,Nx_Dsp_pen(1)+t-1)=alpha_f_a(t-1,Nx_Dsp_pen(1)+t-1)*gamma_g_a(t,Nx_Dsp_pen(1)+t-1); % 环空气体含量
        alpha_l_a(t,Nx_Dsp_pen(1)+t-1)=alpha_f_a(t-1,Nx_Dsp_pen(1)+t-1)*gamma_l_a(t,Nx_Dsp_pen(1)+t-1); % 环空基液含量
        rho_f_a(t,Nx_Dsp_pen(1)+t-1)=rho_g_a(t,Nx_Dsp_pen(1)+t-1)*gamma_g_a(t,Nx_Dsp_pen(1)+t-1)+rho_l_a(t,Nx_Dsp_pen(1)+t-1)*gamma_l_a(t,Nx_Dsp_pen(1)+t-1); % 环空泡沫密度，kg/m^3
        mu_f_a(t,Nx_Dsp_pen(1)+t-1)=mu_g_a(t,Nx_Dsp_pen(1)+t-1)*gamma_g_a(t,Nx_Dsp_pen(1)+t-1)+mu_l_a(t,Nx_Dsp_pen(1)+t-1)*gamma_l_a(t,Nx_Dsp_pen(1)+t-1); % 环空泡沫粘度，Pa*s
        mu_s(t,Nx_Dsp_pen(1)+t-1)=mu_f_a(t,Nx_Dsp_pen(1)+t-1); % 固相粘度，Pa*s
        Va_f_a(t,Nx_Dsp_pen(1)+t-1)=Qm_f_0/(A_a(t,Nx_Dsp_pen(1)+t-1)*rho_f_a(t,Nx_Dsp_pen(1)+t-1)); % 环空泡沫表观流速，m/s
        Va_g_a(t,Nx_Dsp_pen(1)+t-1)=Va_f_a(t,Nx_Dsp_pen(1)+t-1); % 环空气体表观流速，m/s
        Va_l_a(t,Nx_Dsp_pen(1)+t-1)=Va_f_a(t,Nx_Dsp_pen(1)+t-1); % 环空基液表观流速，m/s
        Va_s(t,Nx_Dsp_pen(1)+t-1)=M_s(t)/(A_a(t,Nx_Dsp_pen(1)+t-1)*rho_s); % 岩屑表观速度，m/s
        Vsr(t,Nx_Dsp_pen(1)+t-1)=12*(mu_f_a(t,Nx_Dsp_pen(1)+t-1)/(rho_f_a(t,Nx_Dsp_pen(1)+t-1)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Dsp_pen(1)+t-1))/rho_f_a(t,Nx_Dsp_pen(1)+t-1))*((rho_f_a(t,Nx_Dsp_pen(1)+t-1)*D_s/mu_f_a(t,Nx_Dsp_pen(1)+t-1))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
        alpha_s(t,Nx_Dsp_pen(1)+t-1)=Va_s(t,Nx_Dsp_pen(1)+t-1)/(C0*(Va_s(t,Nx_Dsp_pen(1)+t-1)+Va_f_a(t,Nx_Dsp_pen(1)+t-1))-Vsr(t,Nx_Dsp_pen(1)+t-1));  % 固相体积分数
        V_s(t,Nx_Dsp_pen(1)+t-1)=Va_s(t,Nx_Dsp_pen(1)+t-1)/alpha_s(t,Nx_Dsp_pen(1)+t-1); % 岩屑速度，m/s
        
        alpha_f_a(t,Nx_Dsp_pen(1)+t-1)=1-alpha_s(t,Nx_Dsp_pen(1)+t-1); % 环空泡沫含量
        alpha_g_a(t,Nx_Dsp_pen(1)+t-1)=alpha_f_a(t,Nx_Dsp_pen(1)+t-1)*gamma_g_a(t,Nx_Dsp_pen(1)+t-1); % 环空气体含量
        alpha_l_a(t,Nx_Dsp_pen(1)+t-1)=alpha_f_a(t,Nx_Dsp_pen(1)+t-1)*gamma_l_a(t,Nx_Dsp_pen(1)+t-1); % 环空基液含量
        V_f_a(t,Nx_Dsp_pen(1)+t-1)=Va_f_a(t,Nx_Dsp_pen(1)+t-1)/alpha_f_a(t,Nx_Dsp_pen(1)+t-1); % 环空泡沫流速，m/s
        V_g_a(t,Nx_Dsp_pen(1)+t-1)=V_f_a(t,Nx_Dsp_pen(1)+t-1); % 环空气体流速，m/s
        V_l_a(t,Nx_Dsp_pen(1)+t-1)=V_f_a(t,Nx_Dsp_pen(1)+t-1); % 环空基液流速，m/s
        
        V_m(t,Nx_Dsp_pen(1)+t-1)=Va_s(t,Nx_Dsp_pen(1)+t-1)+Va_f_a(t,Nx_Dsp_pen(1)+t-1); % 环空混合物速度，m/s
        rho_m(t,Nx_Dsp_pen(1)+t-1)=alpha_s(t,Nx_Dsp_pen(1)+t-1)*rho_s+alpha_f_a(t,Nx_Dsp_pen(1)+t-1)*rho_f_a(t,Nx_Dsp_pen(1)+t-1); % 环空混合物密度，kg/m^3
        mu_m(t,Nx_Dsp_pen(1)+t-1)=alpha_s(t,Nx_Dsp_pen(1)+t-1)*mu_s(t,Nx_Dsp_pen(1)+t-1)+alpha_f_a(t,Nx_Dsp_pen(1)+t-1)*mu_f_a(t,Nx_Dsp_pen(1)+t-1); % 环空混合物粘度，Pa*s
        [Ff_a(t,Nx_Dsp_pen(1)+t-1),f_a(t,Nx_Dsp_pen(1)+t-1),Re_a(t,Nx_Dsp_pen(1)+t-1),flow_pattern_a(t,Nx_Dsp_pen(1)+t-1)]=Friction_annulus(rho_m(t,Nx_Dsp_pen(1)+t-1),V_m(t,Nx_Dsp_pen(1)+t-1),mu_m(t,Nx_Dsp_pen(1)+t-1),D_h(t,Nx_Dsp_pen(1)+t-1),h_a,rho_f_a(t,Nx_Dsp_pen(1)+t-1),V_f_a(t,Nx_Dsp_pen(1)+t-1)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
        
        % 第Nx_Dsp_pen(1)+t-2～1个空间节点处相关参数计算
        for x=Nx_Dsp_pen(1)+t-2:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % 环空压力假设值，Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液粘度，Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % 环空基液含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            
            err_NodePressure=1; % 环空压力相对误差
            COUNT_NodePressure=0; % 环空压力迭代次数初值
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % 环空压力迭代次数
                
                % 割线法求解固相体积分数
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % 固相体积分数假设值1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % 固相体积分数假设值2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % 固相体积分数绝对误差
                COUNT_NodeEg=0; % 迭代次数
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % 固相体积分数为alpha_s_ass1时（假设值）
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % 固相速度，m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % 固相表观流速，m/s
                    alpha_f_ass1=1-alpha_s_ass1; % 泡沫体积分数
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % 泡沫速度，m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % 泡沫表观流速，m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % 环空混合物速度，m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % 固相体积分数计算值
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 固相体积分数为alpha_s_ass2时（假设值）
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % 固相质量守恒方程离散后公式中间值计算
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % 固相速度，m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % 固相表观流速，m/s
                    alpha_f_ass2=1-alpha_s_ass2; % 泡沫体积分数
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % 泡沫速度，m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % 泡沫表观流速，m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % 环空混合物速度，m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % 固相体积分数计算值
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 割线法计算固相体积分数
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % 新的固相体积分数假设值alpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % 固相体积分数绝对误差
                    alpha_s_ass1=alpha_s_ass2; % 新的固相体积分数假设值1
                    alpha_s_ass2=alpha_s_ass3; % 新的固相体积分数假设值2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % 将迭代求解得到的真实固相体积分数值赋给alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % 当固相体积分数小于一定值时，认为固相体积数为0，用于防止出现后续的计算误差
                    alpha_s(t,x)=0; % 固相体积分数
                    V_s(t,x)=0; % 固相速度，m/s
                    Va_s(t,x)=0; % 固相表观速度，m/s
                    mu_s(t,x)=0; % 固相粘度，Pa*s
                    Vsr(t,x)=0; % 岩屑沉降末速，m/s
                    
                    alpha_f_a(t,x)=1; % 泡沫体积分数                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                else
                    V_s(t,x)=V_s_ass2; % 固相速度，m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % 固相表观流速，m/s
                    Vsr(t,x)=Vsr_ass2; % 岩屑沉降末速，m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % 泡沫密度，kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % 泡沫体积分数
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    V_f_a(t,x)=V_f_ass2; % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % 环空混合物速度，m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % 环空压力计算值，Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % 环空压力相对误差
                P_a(t,x)=P_new; % 新的环空压力假设值，Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % 出口压力相对误差        
        if (P_a(t,1)-OutPressure)>0 % 根据出口压力误差的正负，对环空管底压力假设值进行调节
            P_a(t,Nx_Dsp_pen(1)+t-1)=P_a(t,Nx_Dsp_pen(1)+t-1)-(P_a(t,1)-OutPressure)/2; % 新的环空管底压力假设值，Pa
        else
            P_a(t,Nx_Dsp_pen(1)+t-1)=P_a(t,Nx_Dsp_pen(1)+t-1)-(P_a(t,1)-OutPressure)/2*0.3; % 新的环空管底压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(1)+t～Nx个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)+t:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1-PHI; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=PHI; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
end

%% 第1～Nt时间节点钻头压降计算
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Dsp_pen(1)+t-1); % 井底泡沫体积流量，m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Dsp_pen(1)+t-1)); % 计算冲砂工具压降（Pa）及射流喷嘴流速V_nozzle（m/s）
end

%% 第1～Nt时间节点相关参数计算（连续管内）
for t=1:1:Nt
    % 管底（第Nx_Dsp_pen(1)+t-1个空间节点）处相关参数计算
    P_ct(t,Nx_Dsp_pen(1)+t-1)=P_a(t,Nx_Dsp_pen(1)+t-1)+delta_P_SWT(t); % 管内压力，Pa
    rho_g_ct(t,Nx_Dsp_pen(1)+t-1)=DensityG(T_ct(t,Nx_Dsp_pen(1)+t-1),P_ct(t,Nx_Dsp_pen(1)+t-1)); % 管内气体密度，kg/m^3
    rho_l_ct(t,Nx_Dsp_pen(1)+t-1)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)+t-1),P_ct(t,Nx_Dsp_pen(1)+t-1)); % 管内基液密度，kg/m^3
    mu_g_ct(t,Nx_Dsp_pen(1)+t-1)=RheologyG(T_ct(t,Nx_Dsp_pen(1)+t-1),P_ct(t,Nx_Dsp_pen(1)+t-1)); % 管内气体粘度，Pa*s
    mu_l_ct(t,Nx_Dsp_pen(1)+t-1)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)+t-1),P_ct(t,Nx_Dsp_pen(1)+t-1)); % 管内基液粘度，Pa*s
    alpha_f_ct(t,Nx_Dsp_pen(1)+t-1)=1; % 管内泡沫含量
    gamma_g_ct(t,Nx_Dsp_pen(1)+t-1)=(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)+t-1))/(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)+t-1)+Qm_l_0/rho_l_ct(t,Nx_Dsp_pen(1)+t-1)); % 管内泡沫质量
    gamma_l_ct(t,Nx_Dsp_pen(1)+t-1)=1-gamma_g_ct(t,Nx_Dsp_pen(1)+t-1); % 管内液体滞留量
    alpha_g_ct(t,Nx_Dsp_pen(1)+t-1)=alpha_f_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_g_ct(t,Nx_Dsp_pen(1)+t-1); % 管内气体含量
    alpha_l_ct(t,Nx_Dsp_pen(1)+t-1)=alpha_f_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_l_ct(t,Nx_Dsp_pen(1)+t-1); % 管内基液含量
    rho_f_ct(t,Nx_Dsp_pen(1)+t-1)=rho_g_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_g_ct(t,Nx_Dsp_pen(1)+t-1)+rho_l_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_l_ct(t,Nx_Dsp_pen(1)+t-1); % 管内泡沫密度，kg/m^3
    mu_f_ct(t,Nx_Dsp_pen(1)+t-1)=mu_g_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_g_ct(t,Nx_Dsp_pen(1)+t-1)+mu_l_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_l_ct(t,Nx_Dsp_pen(1)+t-1); % 管内泡沫粘度，Pa*s
    V_f_ct(t,Nx_Dsp_pen(1)+t-1)=Qm_f_0/(rho_f_ct(t,Nx_Dsp_pen(1)+t-1)*A_ct(t,Nx_Dsp_pen(1)+t-1)); % 管内气体流速，m/s
    [Ff_ct(t,Nx_Dsp_pen(1)+t-1),f_ct(t,Nx_Dsp_pen(1)+t-1),Re_ct(t,Nx_Dsp_pen(1)+t-1),flow_pattern_ct(t,Nx_Dsp_pen(1)+t-1)]=Friction_coiledtubing(rho_f_ct(t,Nx_Dsp_pen(1)+t-1),V_f_ct(t,Nx_Dsp_pen(1)+t-1),mu_f_ct(t,Nx_Dsp_pen(1)+t-1),D_ct_i(t,Nx_Dsp_pen(1)+t-1),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    % 第Nx_Dsp_pen(1)+t-2～1个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)+t-2:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % 管内压力假设值，Pa
        
        err_DriPipePressure=1; % 管内压力相对误差
        COUNT_DriPipePressure=0; % 管内压力迭代次数初值
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % 管内压力迭代次数
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体密度，kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液密度，kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体粘度，Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液粘度，Pa*s
            alpha_f_ct(t,x)=1; % 管内泡沫含量
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % 管内泡沫质量
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % 管内气体含量
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % 管内基液含量
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫密度，kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫粘度，Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % 管内泡沫流速，m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % 管内压力，Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % 计算管内压力假设值与计算值之间的相对误差
            P_ct_ass(t,x)=P_ct(t,x); % 新的管内压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(1)+t～Nx个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)+t:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % 管内气体含量
        alpha_l_ct(t,x)=alpha_l_a(t,x); % 管内基液含量
        gamma_g_ct(t,x)=gamma_g_a(t,x); % 管内泡沫质量
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
        rho_f_ct(t,x)=rho_m(t,x); % 管内气体密度，kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % 管内气体粘度，Pa*s
        V_f_ct(t,t)=V_m(t,x); % 管内气体流速，sm/s
        Ff_ct(t,x)=Ff_a(t,x); % 管内流体单位长度摩擦压降（Pa/m）
        f_ct(t,x)=f_a(t,x); % 管内流体范宁摩擦因子
        Re_ct(t,x)=Re_a(t,x); % 管内流体雷诺数、流体流态
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % 管内流体流态
        P_ct(t,x)=P_a(t,x); % 管内压力，Pa
    end
end

%% 第1～Nt时间节点盘管段出口压力及泵压计算（连续管内）
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % 地面管内流体流速，m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % 地面管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % 盘管段流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % 导向器内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % 井口到注入头顶部段摩擦压降，Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % 导向器段摩擦压降，Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % 导向器到滚筒段摩擦压降，Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % 盘管段摩擦压降，Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % 盘管段出口压力，Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % 泵压，Pa
end

%% 计算井口累积出砂量
M_w_tem(1)=0; % 井口瞬时出砂量初值，kg
M_w_tot(1)=0; % 井口累积出砂量初值，kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口瞬时出砂量，kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口累积出砂量，kg
end

%% 计算末态井筒最大砂浓度
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(2)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% 计算环空泡沫返速平均值V_l_mre
V_f_mre=0; % 环空返速平均值，m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % 环空返速平均值，m/s
    end
end

%% 计算固相沉降末速平均值Vsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % 从Vsr矩阵中提取值非零的沉降末速，m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % 固相沉降末速平均值，m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % 固相沉降末速平均值，m/s
end

%% 判断是否有效冲砂（环空返速平均值大于2倍砂砾沉降末速平均值）、是否完成冲砂
fprintf("冲洗钻进过程2：\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“有效冲砂”
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % 命令行窗口说明当前条件下“完成冲砂”
    else
        fprintf("Sand Cleanout UnFinished!\n"); % 命令行窗口说明当前条件下“未完成冲砂”
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“无效冲砂”
end

%% 井筒ECD计算
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % 环空ECD，kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % 管内ECD，kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % 环空ECD，kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % 管内ECD，kg/m^3
end

%% 冲洗钻进过程末状态数据
alpha_f_a_5=alpha_f_a(Nt,:); % 冲洗钻进过程末状态泡沫含量
alpha_g_a_5=alpha_g_a(Nt,:); % 冲洗钻进过程末状态气相含量
alpha_l_a_5=alpha_l_a(Nt,:); % 冲洗钻进过程末状态液相含量
alpha_s_5=alpha_s(Nt,:); % 冲洗钻进过程末状态固相含量
f_a_5=f_a(Nt,:); % 冲洗钻进过程末状态环空摩擦因子
Ff_a_5=Ff_a(Nt,:); % 冲洗钻进过程末状态环空单位长度摩擦压降，Pa/m
flow_pattern_a_5=flow_pattern_a(Nt,:); % 冲洗钻进过程末状态环空流体流态
gamma_g_a_5=gamma_g_a(Nt,:); % 冲洗钻进过程末状态泡沫质量
gamma_l_a_5=gamma_l_a(Nt,:); % 冲洗钻进过程末状态液体滞留量
mu_f_a_5=mu_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫粘度，Pa*s
mu_g_a_5=mu_g_a(Nt,:); % 冲洗钻进过程末状态气相粘度，Pa*s
mu_l_a_5=mu_l_a(Nt,:); % 冲洗钻进过程末状态环空液体粘度，Pa*s
mu_m_5=mu_m(Nt,:); % 冲洗钻进过程末状态环空混合物粘度，Pa*s
mu_s_5=mu_s(Nt,:); % 冲洗钻进过程末状态固相粘度，Pa*s
P_a_5=P_a(Nt,:); % 冲洗钻进过程末状态环空压力，Pa
Re_a_5=Re_a(Nt,:); % 冲洗钻进过程末状态环空雷诺数
rho_f_a_5=rho_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫密度，kg/m^3
rho_g_a_5=rho_g_a(Nt,:); % 冲洗钻进过程末状态环空气相密度，kg/m^3
rho_l_a_5=rho_l_a(Nt,:); % 冲洗钻进过程末状态环空液相密度，kg/m^3
rho_m_5=rho_m(Nt,:); % 冲洗钻进过程末状态环空混合物密度，kg/m^3
V_f_a_5=V_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫流速，m/s
V_g_a_5=V_g_a(Nt,:); % 冲洗钻进过程末状态气相流速，m/s
V_l_a_5=V_l_a(Nt,:); % 冲洗钻进过程末状态环空液相流速，m/s
V_m_5=V_m(Nt,:); % 冲洗钻进过程末状态环空混合物流速，m/s
V_s_5=V_s(Nt,:); % 冲洗钻进过程末状态固相流速，m/s
Va_f_a_5=Va_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫表观流速，m/s
Va_g_a_5=Va_g_a(Nt,:); % 冲洗钻进过程末状态气相表观流速，m/s
Va_l_a_5=Va_l_a(Nt,:); % 冲洗钻进过程末状态环空液相表观流速，m/s
Va_s_5=Va_s(Nt,:); % 冲洗钻进过程末状态固相表观流速，m/s
Vsr_5=Vsr(Nt,:); % 冲洗钻进过程末状态固相滑移速度，m/s

%% 数据存储（冲洗钻进过程）
ANS_Nt_5=Nt; % 时间节点数
ANS_Time_5=Time(1:Nt); % 冲洗时间，s
ANS_alpha_g_a_5=alpha_g_a(1:Nt,:); % 环空气相体积分数
ANS_alpha_l_a_5=alpha_l_a(1:Nt,:); % 环空液相体积分数
ANS_alpha_s_5=alpha_s(1:Nt,:); % 固相含量
ANS_P_ct_5=P_ct(1:Nt,:); % 管内压力（Pa）
ANS_delta_P_SWT_5=delta_P_SWT(1:Nt); % 冲砂工具压降（Pa）
ANS_P_a_5=P_a(1:Nt,:); % 环空压力（Pa）
ANS_T_a_5=T_a(1:Nt,:); % 环空温度（℃）
ANS_T_ct_5=T_ct(1:Nt,:); % 管内温度（℃）
ANS_P_coil_5=P_coil(1:Nt); % 盘管段出口压力（Pa）
ANS_P_pump_5=P_pump(1:Nt); % 泵压（Pa）
ANS_M_w_tem_5=M_w_tem(1:Nt); % 井口瞬时出砂量（kg）
ANS_M_w_tot_5=M_w_tot(1:Nt); % 井口累积出砂量（kg）
ANS_Va_s_5=Va_s(1:Nt,:); % 岩屑表观速度（m/s）
ANS_Va_f_a_5=Va_f_a(1:Nt,:); % 环空泡沫表观流速（m/s）
ANS_V_s_5=V_s(1:Nt,:); % 岩屑沉降速度（m/s）
ANS_V_f_a_5=V_f_a(1:Nt,:); % 环空泡沫返速（m/s）
ANS_alpha_g_ct_5=alpha_g_ct(1:Nt,:); % 管内气相体积分数
ANS_alpha_l_ct_5=alpha_l_ct(1:Nt,:); % 管内液相体积分数
ANS_gamma_g_ct_5=gamma_g_ct(1:Nt,:); % 管内泡沫质量
ANS_gamma_g_a_5=gamma_g_a(1:Nt,:); % 环空泡沫质量
ANS_ECD_a_5=ECD_a(1:Nt,:); % 环空ECD（kg/m^3）
ANS_L_coil_5=L_coil(1:Nt); % 连续管下深（m）



%% 回拖短起过程2（国际单位制）
%% 变量数值清零
dt=zeros(); % 时间步长，s
dx=zeros(); % 空间步长，m
Time=zeros(); % 连续管下入总时长，s
L_coil=zeros(); % 连续管下深，m
L_reel=zeros(); % 盘管段长度，m
D_t_i=zeros(); % 油管内径，m
D_ct_o=zeros(); % 连续油管外径，m
D_ct_i=zeros(); % 连续油管内径，m
M_s=zeros(); % 井底进砂量，kg/s
P_a=zeros(); % 环空压力，Pa
rho_g_a=zeros(); % 气相密度，kg/m^3
rho_l_a=zeros(); % 环空液相密度，kg/m^3
rho_f_a=zeros(); % 环空泡沫密度，kg/m^3
mu_g_a=zeros(); % 气相粘度，Pa*s
mu_l_a=zeros(); % 环空液体粘度，Pa*s
mu_f_a=zeros(); % 环空泡沫粘度，Pa*s
mu_s=zeros(); % 固相粘度，Pa*s
alpha_g_a=zeros(); % 气相含量
alpha_l_a=zeros(); % 液相含量
alpha_f_a=zeros(); % 泡沫含量
alpha_s=zeros(); % 固相含量
Va_g_a=zeros(); % 气相表观流速，m/s
Va_l_a=zeros(); % 环空液相表观流速，m/s
Va_f_a=zeros(); % 环空泡沫表观流速，m/s
Va_s=zeros(); % 固相表观流速，m/s
V_g_a=zeros(); % 气相流速，m/s
V_l_a=zeros(); % 环空液相流速，m/s
V_f_a=zeros(); % 环空泡沫流速，m/s
V_s=zeros(); % 固相流速，m/s
Vsr=zeros(); % 固相滑移速度，m/s
gamma_g_a=zeros(); % 泡沫质量
gamma_l_a=zeros(); % 液体滞留量
V_m=zeros(); % 环空混合物流速，m/s
rho_m=zeros(); % 环空混合物密度，kg/m^3
mu_m=zeros(); % 环空混合物粘度，Pa*s
Ff_a=zeros(); % 环空单位长度摩擦压降，Pa/m
f_a=zeros(); % 环空摩擦因子
Re_a=zeros(); % 环空雷诺数
flow_pattern_a=zeros(); % 环空流体流态

%% 计算空间步长dx及相应时间步长dt
V6=(Dsp_pen(2)-Lp)/(3*3600); % 回拖短起速度，m/s（输入）
t_6=(Dsp_pen(2)-Lp)/V6; % 回拖短起总时长，s

Nt=Nx_Dsp_pen(2)-Nx_Lp+1; % 时间节点数
nt=Nt-1; % 时间网格数

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % 空间步长，m
end

Time(1)=0; % 连续管起出总时长初值，s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Dsp_pen(2)-t)/V6; % 每起出一个空间步长所需时间，s
    Time(t+1)=Time(t)+dt(t); % 连续管起出至第(t+1)个空间节点所经历的总时长，s
end

%% 计算不同时刻连续管起出深度L_coil及盘管段连续管长度L_reel
L=10000; % 连续油管总长，m（输入）
L_wg=8; % 井口到注入头顶部段连续管长度，m（输入）
L_goose=3; % 导向器段连续管长度，m（输入）
D_goose=2; % 导向器段半径，m（输入）
H_goose=10; % 导向器距地面高度，m（输入）
L_gr=20; % 导向器到滚筒段连续管长度，m（输入）
theta_gr=acosd(H_goose/L_gr); % 导向器到滚筒段连续管与铅垂线夹角，°
D_r_i=3; % 滚筒内径，m（输入）
D_r_o=5; % 滚筒外径，m（输入）
W_r=5; % 滚筒宽度，m（输入）
D_cable=0.005; % 电缆外径，m

L_coil(1)=Dsp_pen(2); % 初始时刻连续管底部深度，m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % 初始时刻盘管段连续管长度，m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)-dx(Nx_Dsp_pen(2)-t+1); % 连续管起出深度，m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % 盘管段连续管长度，m
end
L_cable=L_coil; % 电缆下入长度，m

%% 数据输入及预处理
D_ct_o_0=0.04445; % 连续油管外径，m（输入）
L1=2000; % 近出口第一段连续管长度，m（输入）
D_ct_i_1=0.03709; % 近出口第一段连续管内径，m（输入）
L2=2000; % 近出口第二段连续管长度，m（输入）
D_ct_i_2=0.03653; % 近出口第二段连续管内径，m（输入）
L3=2000; % 近出口第三段连续管长度，m（输入）
D_ct_i_3=0.03555; % 近出口第三段连续管内径，m（输入）
L4=2000; % 近出口第四段连续管长度，m（输入）
D_ct_i_4=0.03489; % 近出口第四段连续管内径，m（输入）
L5=2000; % 近出口第五段连续管长度，m（输入）
D_ct_i_5=0.03409; % 近出口第五段连续管内径，m（输入）

L_t_1=4000; % 上部油管（或套管或裸眼）长度，m（输入）
D_t_i_1=0.09718; %0.068;%0.09718; % 上部油管（或套管或裸眼）内径，m（输入）
L_t_2=2200; % 下部油管（或套管或裸眼）长度，m（输入）
D_t_i_2=0.09718; %0.13970; % 下部油管（或套管或裸眼）内径，m（输入）

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % 油管（或套管或裸眼）内径，m
        else
            D_t_i(t,x)=D_t_i_2; % 油管（或套管或裸眼）内径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(2)
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % 连续油管外径，m
            else
                D_ct_o(t,x)=0; % 连续油管外径，m
            end
        else
            D_ct_o(t,x)=0; % 连续油管外径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(2)
            if Depth(x)<=L_coil(t) % 决定哪几个点可以有内径
                if (L_coil(t)-Depth(x))<=L1 && (L_coil(t)-Depth(x))>=0
                    D_ct_i(t,x)=D_ct_i_1; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % 连续油管内径，m
                end
            else
                D_ct_i(t,x)=0; % 连续油管内径，m
            end
        else
            D_ct_i(t,x)=0; % 连续油管内径，m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % 决定哪几个点可以有数值
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % 连续油管内截面积，m^2
        else
            A_ct(t,x)=0; % 连续油管内截面积，m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % 环空水力直径，m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % 环空截面积，m^2
    end
end

%% 数据输入及预处理
h_t=25.4*10^(-6); % 油管（或套管或裸眼）绝对粗糙度，m（输入）
h_ct=25.4*10^(-6); % 连续油管绝对粗糙度，m（输入）
h_a=(h_t+h_ct)/2; % 环空平均绝对粗糙度，m
epsilon_e=1*10^(-3); % 迭代求解误差限（输入）
epsilon_t=1*10^3; % 最大迭代次数（输入）
g=9.81; % 重力加速度，m/s^2（默认）

T_0=20; % 基液测试温度，℃（输入）
P_0=0.1*10^6; % 基液测试压力，Pa（输入）
rho_l_0=1150; % T_0、P_0下基液密度，kg/m^3（输入）
mu_l_0=0.03; % T_0、P_0下基液粘度，Pa*s（输入）
Qv_l_0=0.1/60; % 基液体积流量，m^3/s（输入）
Qm_l_0=Qv_l_0*rho_l_0; % 基液质量流量，kg/s
rho_g_0=0.655; % 注入气密度，kg/m^3
mu_g_0=RheologyG(T_0,P_0); % 气体粘度，Pa*s
Qv_g_0=4/60; % 注入气体积流量，m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % 注入气质量流量，kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % 泡沫质量流量，kg/s

D_s=1*10^(-3); % 砂砾直径，m（输入）
rho_s=2000; % 砂砾密度，kg/m^3（输入）
H_s=L_s_b-Nx_Dsp_pen(1); % 底部砂床高度，m
PHI=0.6; % 砂床充盈度（输入）
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % 井底砂砾总质量，kg

D_nozzle=4/1000; % 喷嘴直径，m（输入）
N_nozzle=3; % 喷嘴个数（输入）
C=0.95; % 喷嘴流量系数，取0.95（输入）

C0=1.2; % 漂移流方程系数（默认）

M_s(1)=0; % 井底进砂量，kg/s
for t=2:1:Nt
    M_s(t)=0; % 井底进砂量，kg/s
end

OutPressure=1*10^6; % 井口压力，Pa（输入）

%% 温度设置
T_i=20; % 泡沫注入温度，℃
T_g=0.02; % 地温梯度，℃/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % 连续管内泡沫温度（假设每个时刻都一样），℃
    end
end
T_a=T_ct; % 环空泡沫温度（假设每个时刻都一样），℃

%% 第1个时间节点（初始时刻）相关参数插值计算（环空）
for x=1:1:Nx
    P_a(1,x)=P_a_5(x); % 环空压力，Pa
    rho_g_a(1,x)=rho_g_a_5(x); % 气相密度，kg/m^3
    rho_l_a(1,x)=rho_l_a_5(x); % 环空液相密度，kg/m^3
    rho_f_a(1,x)=rho_f_a_5(x); % 环空泡沫密度，kg/m^3
    mu_g_a(1,x)=mu_g_a_5(x); % 气相粘度，Pa*s
    mu_l_a(1,x)=mu_l_a_5(x); % 环空液体粘度，Pa*s
    mu_f_a(1,x)=mu_f_a_5(x); % 环空泡沫粘度，Pa*s
    mu_s(1,x)=mu_s_5(x); % 固相粘度，Pa*s
    alpha_g_a(1,x)=alpha_g_a_5(x); % 气相含量
    alpha_l_a(1,x)=alpha_l_a_5(x); % 液相含量
    alpha_f_a(1,x)=alpha_f_a_5(x); % 泡沫含量
    alpha_s(1,x)=alpha_s_5(1,x); % 固相含量
    Va_g_a(1,x)=Va_g_a_5(x); % 气相表观流速，m/s
    Va_l_a(1,x)=Va_l_a_5(x); % 环空液相表观流速，m/s
    Va_f_a(1,x)=Va_f_a_5(x); % 环空泡沫表观流速，m/s
    Va_s(1,x)=Va_s_5(x); % 固相表观流速，m/s
    V_g_a(1,x)=V_g_a_5(x); % 气相流速，m/s
    V_l_a(1,x)=V_l_a_5(x); % 环空液相流速，m/s
    V_f_a(1,x)=V_f_a_5(x); % 环空泡沫流速，m/s
    V_s(1,x)=V_s_5(x); % 固相流速，m/s
    Vsr(1,x)=Vsr_5(x); % 固相滑移速度，m/s
    gamma_g_a(1,x)=gamma_g_a_5(x); % 泡沫质量
    gamma_l_a(1,x)=gamma_l_a_5(1,x); % 液体滞留量
    V_m(1,x)=V_m_5(x); % 环空混合物流速，m/s
    rho_m(1,x)=rho_m_5(x); % 环空混合物密度，kg/m^3
    mu_m(1,x)=mu_m_5(x); % 环空混合物粘度，Pa*s
    Ff_a(1,x)=Ff_a_5(x); % 环空单位长度摩擦压降，Pa/m
    f_a(1,x)=f_a_5(x); % 环空摩擦因子
    Re_a(1,x)=Re_a_5(x); % 环空雷诺数
    flow_pattern_a(1,x)=flow_pattern_a_5(x); % 环空流体流态
end

%% 第2～Nt个时间节点相关参数计算（环空）
for t=2:1:Nt
    P_a(t,Nx_Dsp_pen(2)-t+1)=P_a(t-1,Nx_Dsp_pen(2)-t+1);  % 环空管底压力假设值，Pa
    
    err_OutPressure=1; % 出口压力相对误差
    COUNT_OutPressure=0; % 出口压力迭代次数初值
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % 出口压力迭代次数
        
        % 管底（第Nx_Dsp_pen(2)-t+1个空间节点）处相关参数计算
        rho_g_a(t,Nx_Dsp_pen(2)-t+1)=DensityG(T_a(t,Nx_Dsp_pen(2)-t+1),P_a(t,Nx_Dsp_pen(2)-t+1)); % 环空气体密度，kg/m^3
        rho_l_a(t,Nx_Dsp_pen(2)-t+1)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(2)-t+1),P_a(t,Nx_Dsp_pen(2)-t+1)); % 环空基液密度，kg/m^3
        mu_g_a(t,Nx_Dsp_pen(2)-t+1)=RheologyG(T_a(t,Nx_Dsp_pen(2)-t+1),P_a(t,Nx_Dsp_pen(2)-t+1)); % 环空气体粘度，Pa*s
        mu_l_a(t,Nx_Dsp_pen(2)-t+1)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(2)-t+1),P_a(t,Nx_Dsp_pen(2)-t+1)); % 环空基液粘度，Pa*s
        gamma_g_a(t,Nx_Dsp_pen(2)-t+1)=(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(2)-t+1))/(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(2)-t+1)+Qm_l_0/rho_l_a(t,Nx_Dsp_pen(2)-t+1)); % 泡沫质量
        gamma_l_a(t,Nx_Dsp_pen(2)-t+1)=1-gamma_g_a(t,Nx_Dsp_pen(2)-t+1); % 液体滞留量
        alpha_g_a(t,Nx_Dsp_pen(2)-t+1)=alpha_f_a(t-1,Nx_Dsp_pen(2)-t+1)*gamma_g_a(t,Nx_Dsp_pen(2)-t+1); % 环空气体含量
        alpha_l_a(t,Nx_Dsp_pen(2)-t+1)=alpha_f_a(t-1,Nx_Dsp_pen(2)-t+1)*gamma_l_a(t,Nx_Dsp_pen(2)-t+1); % 环空基液含量
        rho_f_a(t,Nx_Dsp_pen(2)-t+1)=rho_g_a(t,Nx_Dsp_pen(2)-t+1)*gamma_g_a(t,Nx_Dsp_pen(2)-t+1)+rho_l_a(t,Nx_Dsp_pen(2)-t+1)*gamma_l_a(t,Nx_Dsp_pen(2)-t+1); % 环空泡沫密度，kg/m^3
        mu_f_a(t,Nx_Dsp_pen(2)-t+1)=mu_g_a(t,Nx_Dsp_pen(2)-t+1)*gamma_g_a(t,Nx_Dsp_pen(2)-t+1)+mu_l_a(t,Nx_Dsp_pen(2)-t+1)*gamma_l_a(t,Nx_Dsp_pen(2)-t+1); % 环空泡沫粘度，Pa*s
        mu_s(t,Nx_Dsp_pen(2)-t+1)=mu_f_a(t,Nx_Dsp_pen(2)-t+1); % 固相粘度，Pa*s
        Va_f_a(t,Nx_Dsp_pen(2)-t+1)=Qm_f_0/(A_a(t,Nx_Dsp_pen(2)-t+1)*rho_f_a(t,Nx_Dsp_pen(2)-t+1)); % 环空泡沫表观流速，m/s
        Va_g_a(t,Nx_Dsp_pen(2)-t+1)=Va_f_a(t,Nx_Dsp_pen(2)-t+1); % 环空气体表观流速，m/s
        Va_l_a(t,Nx_Dsp_pen(2)-t+1)=Va_f_a(t,Nx_Dsp_pen(2)-t+1); % 环空基液表观流速，m/s
        Va_s(t,Nx_Dsp_pen(2)-t+1)=M_s(t)/(A_a(t,Nx_Dsp_pen(2)-t+1)*rho_s); % 岩屑表观速度，m/s
        Vsr(t,Nx_Dsp_pen(2)-t+1)=12*(mu_f_a(t,Nx_Dsp_pen(2)-t+1)/(rho_f_a(t,Nx_Dsp_pen(2)-t+1)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Dsp_pen(2)-t+1))/rho_f_a(t,Nx_Dsp_pen(2)-t+1))*((rho_f_a(t,Nx_Dsp_pen(2)-t+1)*D_s/mu_f_a(t,Nx_Dsp_pen(2)-t+1))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
        alpha_s(t,Nx_Dsp_pen(2)-t+1)=Va_s(t,Nx_Dsp_pen(2)-t+1)/(C0*(Va_s(t,Nx_Dsp_pen(2)-t+1)+Va_f_a(t,Nx_Dsp_pen(2)-t+1))-Vsr(t,Nx_Dsp_pen(2)-t+1));  % 固相体积分数
        V_s(t,Nx_Dsp_pen(2)-t+1)=0;%Va_s(t,Nx_Dsp_pen(2)-t+1)/alpha_s(t,Nx_Dsp_pen(2)-t+1); % 岩屑速度，m/s
        
        alpha_f_a(t,Nx_Dsp_pen(2)-t+1)=1-alpha_s(t,Nx_Dsp_pen(2)-t+1); % 环空泡沫含量
        alpha_g_a(t,Nx_Dsp_pen(2)-t+1)=alpha_f_a(t,Nx_Dsp_pen(2)-t+1)*gamma_g_a(t,Nx_Dsp_pen(2)-t+1); % 环空气体含量
        alpha_l_a(t,Nx_Dsp_pen(2)-t+1)=alpha_f_a(t,Nx_Dsp_pen(2)-t+1)*gamma_l_a(t,Nx_Dsp_pen(2)-t+1); % 环空基液含量
        V_f_a(t,Nx_Dsp_pen(2)-t+1)=Va_f_a(t,Nx_Dsp_pen(2)-t+1)/alpha_f_a(t,Nx_Dsp_pen(2)-t+1); % 环空泡沫流速，m/s
        V_g_a(t,Nx_Dsp_pen(2)-t+1)=V_f_a(t,Nx_Dsp_pen(2)-t+1); % 环空气体流速，m/s
        V_l_a(t,Nx_Dsp_pen(2)-t+1)=V_f_a(t,Nx_Dsp_pen(2)-t+1); % 环空基液流速，m/s
        
        V_m(t,Nx_Dsp_pen(2)-t+1)=Va_s(t,Nx_Dsp_pen(2)-t+1)+Va_f_a(t,Nx_Dsp_pen(2)-t+1); % 环空混合物速度，m/s
        rho_m(t,Nx_Dsp_pen(2)-t+1)=alpha_s(t,Nx_Dsp_pen(2)-t+1)*rho_s+alpha_f_a(t,Nx_Dsp_pen(2)-t+1)*rho_f_a(t,Nx_Dsp_pen(2)-t+1); % 环空混合物密度，kg/m^3
        mu_m(t,Nx_Dsp_pen(2)-t+1)=alpha_s(t,Nx_Dsp_pen(2)-t+1)*mu_s(t,Nx_Dsp_pen(2)-t+1)+alpha_f_a(t,Nx_Dsp_pen(2)-t+1)*mu_f_a(t,Nx_Dsp_pen(2)-t+1); % 环空混合物粘度，Pa*s
        [Ff_a(t,Nx_Dsp_pen(2)-t+1),f_a(t,Nx_Dsp_pen(2)-t+1),Re_a(t,Nx_Dsp_pen(2)-t+1),flow_pattern_a(t,Nx_Dsp_pen(2)-t+1)]=Friction_annulus(rho_m(t,Nx_Dsp_pen(2)-t+1),V_m(t,Nx_Dsp_pen(2)-t+1),mu_m(t,Nx_Dsp_pen(2)-t+1),D_h(t,Nx_Dsp_pen(2)-t+1),h_a,rho_f_a(t,Nx_Dsp_pen(2)-t+1),V_f_a(t,Nx_Dsp_pen(2)-t+1)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
        
        % 第Nx_Dsp_pen(2)-t～1个空间节点处相关参数计算
        for x=Nx_Dsp_pen(2)-t:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % 环空压力假设值，Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液粘度，Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % 环空基液含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            
            err_NodePressure=1; % 环空压力相对误差
            COUNT_NodePressure=0; % 环空压力迭代次数初值
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % 环空压力迭代次数
                
                % 割线法求解固相体积分数
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % 固相体积分数假设值1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % 固相体积分数假设值2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % 固相体积分数绝对误差
                COUNT_NodeEg=0; % 迭代次数
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % 固相体积分数为alpha_s_ass1时（假设值）
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % 固相速度，m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % 固相表观流速，m/s
                    alpha_f_ass1=1-alpha_s_ass1; % 泡沫体积分数
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % 泡沫速度，m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % 泡沫表观流速，m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % 环空混合物速度，m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % 固相体积分数计算值
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 固相体积分数为alpha_s_ass2时（假设值）
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % 固相质量守恒方程离散后公式中间值计算
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % 固相速度，m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % 固相表观流速，m/s
                    alpha_f_ass2=1-alpha_s_ass2; % 泡沫体积分数
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % 泡沫速度，m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % 泡沫表观流速，m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % 环空混合物速度，m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % 固相体积分数计算值
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 割线法计算固相体积分数
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % 新的固相体积分数假设值alpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % 固相体积分数绝对误差
                    alpha_s_ass1=alpha_s_ass2; % 新的固相体积分数假设值1
                    alpha_s_ass2=alpha_s_ass3; % 新的固相体积分数假设值2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % 将迭代求解得到的真实固相体积分数值赋给alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % 当固相体积分数小于一定值时，认为固相体积数为0，用于防止出现后续的计算误差
                    alpha_s(t,x)=0; % 固相体积分数
                    V_s(t,x)=0; % 固相速度，m/s
                    Va_s(t,x)=0; % 固相表观速度，m/s
                    mu_s(t,x)=0; % 固相粘度，Pa*s
                    Vsr(t,x)=0; % 岩屑沉降末速，m/s
                    
                    alpha_f_a(t,x)=1; % 泡沫体积分数                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                else
                    V_s(t,x)=V_s_ass2; % 固相速度，m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % 固相表观流速，m/s
                    Vsr(t,x)=Vsr_ass2; % 岩屑沉降末速，m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % 泡沫密度，kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % 泡沫体积分数
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    V_f_a(t,x)=V_f_ass2; % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % 环空混合物速度，m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % 环空压力计算值，Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % 环空压力相对误差
                P_a(t,x)=P_new; % 新的环空压力假设值，Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % 出口压力相对误差
        if (P_a(t,1)-OutPressure)>0 % 根据出口压力误差的正负，对环空管底压力假设值进行调节
            P_a(t,Nx_Dsp_pen(2)-t+1)=P_a(t,Nx_Dsp_pen(2)-t+1)-(P_a(t,1)-OutPressure)/2; % 新的环空管底压力假设值，Pa
        else
            P_a(t,Nx_Dsp_pen(2)-t+1)=P_a(t,Nx_Dsp_pen(2)-t+1)-(P_a(t,1)-OutPressure)/2*0.3; % 新的环空管底压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(2)-t+2～Nx个空间节点处相关参数计算
    for x=Nx_Dsp_pen(2)-t+2:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=0; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
end

%% 第1～Nt时间节点钻头压降计算
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Dsp_pen(2)-t+1); % 井底泡沫体积流量，m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Dsp_pen(2)-t+1)); % 计算冲砂工具压降（Pa）及射流喷嘴流速V_nozzle（m/s）
end

%% 第1～Nt时间节点相关参数计算（连续管内）
for t=1:1:Nt
    % 管底（第Nx_Dsp_pen(2)-t+1个空间节点）处相关参数计算
    P_ct(t,Nx_Dsp_pen(2)-t+1)=P_a(t,Nx_Dsp_pen(2)-t+1)+delta_P_SWT(t); % 管内压力，Pa
    rho_g_ct(t,Nx_Dsp_pen(2)-t+1)=DensityG(T_ct(t,Nx_Dsp_pen(2)-t+1),P_ct(t,Nx_Dsp_pen(2)-t+1)); % 管内气体密度，kg/m^3
    rho_l_ct(t,Nx_Dsp_pen(2)-t+1)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(2)-t+1),P_ct(t,Nx_Dsp_pen(2)-t+1)); % 管内基液密度，kg/m^3
    mu_g_ct(t,Nx_Dsp_pen(2)-t+1)=RheologyG(T_ct(t,Nx_Dsp_pen(2)-t+1),P_ct(t,Nx_Dsp_pen(2)-t+1)); % 管内气体粘度，Pa*s
    mu_l_ct(t,Nx_Dsp_pen(2)-t+1)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(2)-t+1),P_ct(t,Nx_Dsp_pen(2)-t+1)); % 管内基液粘度，Pa*s
    alpha_f_ct(t,Nx_Dsp_pen(2)-t+1)=1; % 管内泡沫含量
    gamma_g_ct(t,Nx_Dsp_pen(2)-t+1)=(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(2)-t+1))/(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(2)-t+1)+Qm_l_0/rho_l_ct(t,Nx_Dsp_pen(2)-t+1)); % 管内泡沫质量
    gamma_l_ct(t,Nx_Dsp_pen(2)-t+1)=1-gamma_g_ct(t,Nx_Dsp_pen(2)-t+1); % 管内液体滞留量
    alpha_g_ct(t,Nx_Dsp_pen(2)-t+1)=alpha_f_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(2)-t+1); % 管内气体含量
    alpha_l_ct(t,Nx_Dsp_pen(2)-t+1)=alpha_f_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(2)-t+1); % 管内基液含量
    rho_f_ct(t,Nx_Dsp_pen(2)-t+1)=rho_g_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(2)-t+1)+rho_l_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(2)-t+1); % 管内泡沫密度，kg/m^3
    mu_f_ct(t,Nx_Dsp_pen(2)-t+1)=mu_g_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(2)-t+1)+mu_l_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(2)-t+1); % 管内泡沫粘度，Pa*s
    V_f_ct(t,Nx_Dsp_pen(2)-t+1)=Qm_f_0/(rho_f_ct(t,Nx_Dsp_pen(2)-t+1)*A_ct(t,Nx_Dsp_pen(2)-t+1)); % 管内气体流速，m/s
    [Ff_ct(t,Nx_Dsp_pen(2)-t+1),f_ct(t,Nx_Dsp_pen(2)-t+1),Re_ct(t,Nx_Dsp_pen(2)-t+1),flow_pattern_ct(t,Nx_Dsp_pen(2)-t+1)]=Friction_coiledtubing(rho_f_ct(t,Nx_Dsp_pen(2)-t+1),V_f_ct(t,Nx_Dsp_pen(2)-t+1),mu_f_ct(t,Nx_Dsp_pen(2)-t+1),D_ct_i(t,Nx_Dsp_pen(2)-t+1),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    % 第Nx_Dsp_pen(2)-t～1个空间节点处相关参数计算
    for x=Nx_Dsp_pen(2)-t:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % 管内压力假设值，Pa
        
        err_DriPipePressure=1; % 管内压力相对误差
        COUNT_DriPipePressure=0; % 管内压力迭代次数初值
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % 管内压力迭代次数
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体密度，kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液密度，kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体粘度，Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液粘度，Pa*s
            alpha_f_ct(t,x)=1; % 管内泡沫含量
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % 管内泡沫质量
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % 管内气体含量
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % 管内基液含量
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫密度，kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫粘度，Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % 管内泡沫流速，m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % 管内压力，Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % 计算管内压力假设值与计算值之间的相对误差
            P_ct_ass(t,x)=P_ct(t,x); % 新的管内压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(2)-t+2～Nx个空间节点处相关参数计算
    for x=Nx_Dsp_pen(2)-t+2:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % 管内气体含量
        alpha_l_ct(t,x)=alpha_l_a(t,x); % 管内基液含量
        gamma_g_ct(t,x)=gamma_g_a(t,x); % 管内泡沫质量
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
        rho_f_ct(t,x)=rho_m(t,x); % 管内气体密度，kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % 管内气体粘度，Pa*s
        V_f_ct(t,t)=V_m(t,x); % 管内气体流速，sm/s
        Ff_ct(t,x)=Ff_a(t,x); % 管内流体单位长度摩擦压降（Pa/m）
        f_ct(t,x)=f_a(t,x); % 管内流体范宁摩擦因子
        Re_ct(t,x)=Re_a(t,x); % 管内流体雷诺数、流体流态
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % 管内流体流态
        P_ct(t,x)=P_a(t,x); % 管内压力，Pa
    end
end

%% 第1～Nt时间节点盘管段出口压力及泵压计算（连续管内）
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % 地面管内流体流速，m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % 地面管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % 盘管段流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % 导向器内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % 井口到注入头顶部段摩擦压降，Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % 导向器段摩擦压降，Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % 导向器到滚筒段摩擦压降，Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % 盘管段摩擦压降，Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % 盘管段出口压力，Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % 泵压，Pa
end

%% 计算井口累积出砂量
M_w_tem(1)=0; % 井口瞬时出砂量初值，kg
M_w_tot(1)=0; % 井口累积出砂量初值，kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口瞬时出砂量，kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口累积出砂量，kg
end

%% 计算末态井筒最大砂浓度
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(1)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% 计算环空泡沫返速平均值V_l_mre
V_f_mre=0; % 环空返速平均值，m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % 环空返速平均值，m/s
    end
end

%% 计算固相沉降末速平均值Vsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % 从Vsr矩阵中提取值非零的沉降末速，m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % 固相沉降末速平均值，m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % 固相沉降末速平均值，m/s
end

%% 判断是否有效冲砂（环空返速平均值大于2倍砂砾沉降末速平均值）、是否完成冲砂
fprintf("回拖短起过程2：\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“有效冲砂”
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % 命令行窗口说明当前条件下“完成冲砂”
    else
        fprintf("Sand Cleanout UnFinished!\n"); % 命令行窗口说明当前条件下“未完成冲砂”
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“无效冲砂”
end

%% 井筒ECD计算
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % 环空ECD，kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % 管内ECD，kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % 环空ECD，kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % 管内ECD，kg/m^3
end

%% 回拖短起过程末状态数据
alpha_f_a_6=alpha_f_a(Nt,:); % 冲洗钻进过程末状态泡沫含量
alpha_g_a_6=alpha_g_a(Nt,:); % 冲洗钻进过程末状态气相含量
alpha_l_a_6=alpha_l_a(Nt,:); % 冲洗钻进过程末状态液相含量
alpha_s_6=alpha_s(Nt,:); % 冲洗钻进过程末状态固相含量
f_a_6=f_a(Nt,:); % 冲洗钻进过程末状态环空摩擦因子
Ff_a_6=Ff_a(Nt,:); % 冲洗钻进过程末状态环空单位长度摩擦压降，Pa/m
flow_pattern_a_6=flow_pattern_a(Nt,:); % 冲洗钻进过程末状态环空流体流态
gamma_g_a_6=gamma_g_a(Nt,:); % 冲洗钻进过程末状态泡沫质量
gamma_l_a_6=gamma_l_a(Nt,:); % 冲洗钻进过程末状态液体滞留量
mu_f_a_6=mu_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫粘度，Pa*s
mu_g_a_6=mu_g_a(Nt,:); % 冲洗钻进过程末状态气相粘度，Pa*s
mu_l_a_6=mu_l_a(Nt,:); % 冲洗钻进过程末状态环空液体粘度，Pa*s
mu_m_6=mu_m(Nt,:); % 冲洗钻进过程末状态环空混合物粘度，Pa*s
mu_s_6=mu_s(Nt,:); % 冲洗钻进过程末状态固相粘度，Pa*s
P_a_6=P_a(Nt,:); % 冲洗钻进过程末状态环空压力，Pa
Re_a_6=Re_a(Nt,:); % 冲洗钻进过程末状态环空雷诺数
rho_f_a_6=rho_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫密度，kg/m^3
rho_g_a_6=rho_g_a(Nt,:); % 冲洗钻进过程末状态环空气相密度，kg/m^3
rho_l_a_6=rho_l_a(Nt,:); % 冲洗钻进过程末状态环空液相密度，kg/m^3
rho_m_6=rho_m(Nt,:); % 冲洗钻进过程末状态环空混合物密度，kg/m^3
V_f_a_6=V_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫流速，m/s
V_g_a_6=V_g_a(Nt,:); % 冲洗钻进过程末状态气相流速，m/s
V_l_a_6=V_l_a(Nt,:); % 冲洗钻进过程末状态环空液相流速，m/s
V_m_6=V_m(Nt,:); % 冲洗钻进过程末状态环空混合物流速，m/s
V_s_6=V_s(Nt,:); % 冲洗钻进过程末状态固相流速，m/s
Va_f_a_6=Va_f_a(Nt,:); % 冲洗钻进过程末状态环空泡沫表观流速，m/s
Va_g_a_6=Va_g_a(Nt,:); % 冲洗钻进过程末状态气相表观流速，m/s
Va_l_a_6=Va_l_a(Nt,:); % 冲洗钻进过程末状态环空液相表观流速，m/s
Va_s_6=Va_s(Nt,:); % 冲洗钻进过程末状态固相表观流速，m/s
Vsr_6=Vsr(Nt,:); % 冲洗钻进过程末状态固相滑移速度，m/s

%% 数据存储（回拖短起过程）
ANS_Nt_6=Nt; % 时间节点数
ANS_Time_6=Time(1:Nt); % 连续管起出时间，s
ANS_alpha_g_a_6=alpha_g_a(1:Nt,:); % 环空气相体积分数
ANS_alpha_l_a_6=alpha_l_a(1:Nt,:); % 环空液相体积分数
ANS_alpha_s_6=alpha_s(1:Nt,:); % 固相含量
ANS_P_ct_6=P_ct(1:Nt,:); % 管内压力（Pa）
ANS_delta_P_SWT_6=delta_P_SWT(1:Nt); % 冲砂工具压降（Pa）
ANS_P_a_6=P_a(1:Nt,:); % 环空压力（Pa）
ANS_T_a_6=T_a(1:Nt,:); % 环空温度（℃）
ANS_T_ct_6=T_ct(1:Nt,:); % 管内温度（℃）
ANS_P_coil_6=P_coil(1:Nt); % 盘管段出口压力（Pa）
ANS_P_pump_6=P_pump(1:Nt); % 泵压（Pa）
ANS_M_w_tem_6=M_w_tem(1:Nt); % 井口瞬时出砂量（kg）
ANS_M_w_tot_6=M_w_tot(1:Nt); % 井口累积出砂量（kg）
ANS_Va_s_6=Va_s(1:Nt,:); % 岩屑表观速度（m/s）
ANS_Va_f_a_6=Va_f_a(1:Nt,:); % 环空泡沫表观流速（m/s）
ANS_V_s_6=V_s(1:Nt,:); % 岩屑沉降速度（m/s）
ANS_V_f_a_6=V_f_a(1:Nt,:); % 环空泡沫返速（m/s）
ANS_alpha_g_ct_6=alpha_g_ct(1:Nt,:); % 管内气相体积分数
ANS_alpha_l_ct_6=alpha_l_ct(1:Nt,:); % 管内液相体积分数
ANS_gamma_g_ct_6=gamma_g_ct(1:Nt,:); % 管内泡沫质量
ANS_gamma_g_a_6=gamma_g_a(1:Nt,:); % 环空泡沫质量
ANS_ECD_a_6=ECD_a(1:Nt,:); % 环空ECD（kg/m^3）
ANS_L_coil_6=L_coil(1:Nt); % 连续管下深（m）



%% 定点循环过程2（国际单位制）
%% 变量数值清零
dt=zeros(); % 时间步长，s
dx=zeros(); % 空间步长，m
Time=zeros(); % 连续管下入总时长，s
L_coil=zeros(); % 连续管下深，m
L_reel=zeros(); % 盘管段长度，m
D_t_i=zeros(); % 油管内径，m
D_ct_o=zeros(); % 连续油管外径，m
D_ct_i=zeros(); % 连续油管内径，m
M_s=zeros(); % 井底进砂量，kg/s
P_a=zeros(); % 环空压力，Pa
rho_g_a=zeros(); % 气相密度，kg/m^3
rho_l_a=zeros(); % 环空液相密度，kg/m^3
rho_f_a=zeros(); % 环空泡沫密度，kg/m^3
mu_g_a=zeros(); % 气相粘度，Pa*s
mu_l_a=zeros(); % 环空液体粘度，Pa*s
mu_f_a=zeros(); % 环空泡沫粘度，Pa*s
mu_s=zeros(); % 固相粘度，Pa*s
alpha_g_a=zeros(); % 气相含量
alpha_l_a=zeros(); % 液相含量
alpha_f_a=zeros(); % 泡沫含量
alpha_s=zeros(); % 固相含量
Va_g_a=zeros(); % 气相表观流速，m/s
Va_l_a=zeros(); % 环空液相表观流速，m/s
Va_f_a=zeros(); % 环空泡沫表观流速，m/s
Va_s=zeros(); % 固相表观流速，m/s
V_g_a=zeros(); % 气相流速，m/s
V_l_a=zeros(); % 环空液相流速，m/s
V_f_a=zeros(); % 环空泡沫流速，m/s
V_s=zeros(); % 固相流速，m/s
Vsr=zeros(); % 固相滑移速度，m/s
gamma_g_a=zeros(); % 泡沫质量
gamma_l_a=zeros(); % 液体滞留量
V_m=zeros(); % 环空混合物流速，m/s
rho_m=zeros(); % 环空混合物密度，kg/m^3
mu_m=zeros(); % 环空混合物粘度，Pa*s
Ff_a=zeros(); % 环空单位长度摩擦压降，Pa/m
f_a=zeros(); % 环空摩擦因子
Re_a=zeros(); % 环空雷诺数
flow_pattern_a=zeros(); % 环空流体流态

%% 计算空间步长dx及相应时间步长dt
t_7=3600*1; % 定点循环总时长，s

Nt=101; % 时间节点数（输入）
nt=Nt-1; % 时间网格数

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % 空间步长，m
end

for t=1:1:Nt-1
    dt(t)=t_7/nt; % 时间步长
end

Time(1)=0; % 定点循环总时长，s
for t=1:1:Nt-1
    Time(t+1)=Time(t)+dt(t); % 定点循环至第(t+1)个时刻所经历的总时长，s
end

%% 计算不同时刻连续管下入深度L_coil及盘管段连续管长度L_reel
L=10000; % 连续油管总长，m（输入）
L_wg=8; % 井口到注入头顶部段连续管长度，m（输入）
L_goose=3; % 导向器段连续管长度，m（输入）
D_goose=2; % 导向器段半径，m（输入）
H_goose=10; % 导向器距地面高度，m（输入）
L_gr=20; % 导向器到滚筒段连续管长度，m（输入）
theta_gr=acosd(H_goose/L_gr); % 导向器到滚筒段连续管与铅垂线夹角，°
D_r_i=3; % 滚筒内径，m（输入）
D_r_o=5; % 滚筒外径，m（输入）
W_r=5; % 滚筒宽度，m（输入）
D_cable=0.005; % 电缆外径，m

L_coil(1)=Lp; % 连续管下入深度，m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % 初始时刻盘管段连续管长度，m
for t=2:1:Nt
    L_coil(t)=Lp; % 连续管下入深度，m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % 盘管段连续管长度，m
end
L_cable=L_coil; % 电缆下入长度，m

%% 数据输入及预处理
D_ct_o_0=0.04445; % 连续油管外径，m（输入）
L1=2000; % 近出口第一段连续管长度，m（输入）
D_ct_i_1=0.03709; % 近出口第一段连续管内径，m（输入）
L2=2000; % 近出口第二段连续管长度，m（输入）
D_ct_i_2=0.03653; % 近出口第二段连续管内径，m（输入）
L3=2000; % 近出口第三段连续管长度，m（输入）
D_ct_i_3=0.03555; % 近出口第三段连续管内径，m（输入）
L4=2000; % 近出口第四段连续管长度，m（输入）
D_ct_i_4=0.03489; % 近出口第四段连续管内径，m（输入）
L5=2000; % 近出口第五段连续管长度，m（输入）
D_ct_i_5=0.03409; % 近出口第五段连续管内径，m（输入）

L_t_1=4000; % 上部油管（或套管或裸眼）长度，m（输入）
D_t_i_1=0.09718; %0.068;% 上部油管（或套管或裸眼）内径，m（输入）
L_t_2=2200; % 下部油管（或套管或裸眼）长度，m（输入）
D_t_i_2=0.09718; %0.13970; % 下部油管（或套管或裸眼）内径，m（输入）

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % 油管（或套管或裸眼）内径，m
        else
            D_t_i(t,x)=D_t_i_2; % 油管（或套管或裸眼）内径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            D_ct_o(t,x)=D_ct_o_0; % 连续油管外径，m
        else
            D_ct_o(t,x)=0; % 连续油管外径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            if (L_coil(t)-Depth(x))<=L1 && (L_coil(t)-Depth(x))>=0
                D_ct_i(t,x)=D_ct_i_1; % 连续油管内径，m
            elseif (L_coil(t)-Depth(x))<=L1+L2
                D_ct_i(t,x)=D_ct_i_2; % 连续油管内径，m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                D_ct_i(t,x)=D_ct_i_3; % 连续油管内径，m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                D_ct_i(t,x)=D_ct_i_4; % 连续油管内径，m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                D_ct_i(t,x)=D_ct_i_5; % 连续油管内径，m
            end
        else
            D_ct_i(t,x)=0; % 连续油管内径，m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % 决定哪几个点可以有数值
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % 连续油管内截面积，m^2
        else
            A_ct(t,x)=0; % 连续油管内截面积，m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % 环空水力直径，m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % 环空截面积，m^2
    end
end

%% 数据输入及预处理
h_t=25.4*10^(-6); % 油管（或套管或裸眼）绝对粗糙度，m（输入）
h_ct=25.4*10^(-6); % 连续油管绝对粗糙度，m（输入）
h_a=(h_t+h_ct)/2; % 环空平均绝对粗糙度，m
epsilon_e=1*10^(-3); % 迭代求解误差限（输入）
epsilon_t=1*10^3; % 最大迭代次数（输入）
g=9.81; % 重力加速度，m/s^2（默认）

T_0=20; % 基液测试温度，℃（输入）
P_0=0.1*10^6; % 基液测试压力，Pa（输入）
rho_l_0=1150; % T_0、P_0下基液密度，kg/m^3（输入）
mu_l_0=0.03; % T_0、P_0下基液粘度，Pa*s（输入）
Qv_l_0=0.1/60; % 基液体积流量，m^3/s（输入）
Qm_l_0=Qv_l_0*rho_l_0; % 基液质量流量，kg/s
rho_g_0=0.655; % 注入气密度，kg/m^3
mu_g_0=RheologyG(T_0,P_0); % 气体粘度，Pa*s
Qv_g_0=4/60; % 注入气体积流量，m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % 注入气质量流量，kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % 泡沫质量流量，kg/s

D_s=1*10^(-3); % 砂砾直径，m（输入）
rho_s=2000; % 砂砾密度，kg/m^3（输入）
H_s=L_s_b-Dsp_pen(1); % 底部砂床高度，m
PHI=0.6; % 砂床充盈度（输入）
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % 井底砂砾总质量，kg

D_nozzle=4/1000; % 喷嘴直径，m（输入）
N_nozzle=3; % 喷嘴个数（输入）
C=0.95; % 喷嘴流量系数，取0.95（输入）

C0=1.2; % 漂移流方程系数（默认）

M_s(1)=0; % 井底进砂量，kg/s
for t=2:1:Nt
    M_s(t)=0; % 井底进砂量，kg/s
end

OutPressure=1*10^6; % 井口压力，Pa（输入）

%% 温度设置
T_i=20; % 泡沫注入温度，℃
T_g=0.02; % 地温梯度，℃/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % 连续管内泡沫温度（假设每个时刻都一样），℃
    end
end
T_a=T_ct; % 环空泡沫温度（假设每个时刻都一样），℃

%% 第1个时间节点（初始时刻）相关参数插值计算（环空）
for x=1:1:Nx
    P_a(1,x)=P_a_6(x); % 环空压力，Pa
    rho_g_a(1,x)=rho_g_a_6(x); % 气相密度，kg/m^3
    rho_l_a(1,x)=rho_l_a_6(x); % 环空液相密度，kg/m^3
    rho_f_a(1,x)=rho_f_a_6(x); % 环空泡沫密度，kg/m^3
    mu_g_a(1,x)=mu_g_a_6(x); % 气相粘度，Pa*s
    mu_l_a(1,x)=mu_l_a_6(x); % 环空液体粘度，Pa*s
    mu_f_a(1,x)=mu_f_a_6(x); % 环空泡沫粘度，Pa*s
    mu_s(1,x)=mu_s_6(x); % 固相粘度，Pa*s
    alpha_g_a(1,x)=alpha_g_a_6(x); % 气相含量
    alpha_l_a(1,x)=alpha_l_a_6(x); % 液相含量
    alpha_f_a(1,x)=alpha_f_a_6(x); % 泡沫含量
    alpha_s(1,x)=alpha_s_6(1,x); % 固相含量
    Va_g_a(1,x)=Va_g_a_6(x); % 气相表观流速，m/s
    Va_l_a(1,x)=Va_l_a_6(x); % 环空液相表观流速，m/s
    Va_f_a(1,x)=Va_f_a_6(x); % 环空泡沫表观流速，m/s
    Va_s(1,x)=Va_s_6(x); % 固相表观流速，m/s
    V_g_a(1,x)=V_g_a_6(x); % 气相流速，m/s
    V_l_a(1,x)=V_l_a_6(x); % 环空液相流速，m/s
    V_f_a(1,x)=V_f_a_6(x); % 环空泡沫流速，m/s
    V_s(1,x)=V_s_6(x); % 固相流速，m/s
    Vsr(1,x)=Vsr_6(x); % 固相滑移速度，m/s
    gamma_g_a(1,x)=gamma_g_a_6(x); % 泡沫质量
    gamma_l_a(1,x)=gamma_l_a_6(1,x); % 液体滞留量
    V_m(1,x)=V_m_6(x); % 环空混合物流速，m/s
    rho_m(1,x)=rho_m_6(x); % 环空混合物密度，kg/m^3
    mu_m(1,x)=mu_m_6(x); % 环空混合物粘度，Pa*s
    Ff_a(1,x)=Ff_a_6(x); % 环空单位长度摩擦压降，Pa/m
    f_a(1,x)=f_a_6(x); % 环空摩擦因子
    Re_a(1,x)=Re_a_6(x); % 环空雷诺数
    flow_pattern_a(1,x)=flow_pattern_a_6(x); % 环空流体流态
end

%% 第2～Nt个时间节点相关参数计算（环空）
for t=2:1:Nt
    P_a(t,Nx_Lp)=P_a(t-1,Nx_Lp);  % 环空管底压力假设值，Pa
    
    err_OutPressure=1; % 出口压力相对误差
    COUNT_OutPressure=0; % 出口压力迭代次数初值
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % 出口压力迭代次数
        
        % 管底（第Nx_Lp个空间节点）处相关参数计算
        rho_g_a(t,Nx_Lp)=DensityG(T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % 环空气体密度，kg/m^3
        rho_l_a(t,Nx_Lp)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % 环空基液密度，kg/m^3
        mu_g_a(t,Nx_Lp)=RheologyG(T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % 环空气体粘度，Pa*s
        mu_l_a(t,Nx_Lp)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % 环空基液粘度，Pa*s
        gamma_g_a(t,Nx_Lp)=(Qm_g_0/rho_g_a(t,Nx_Lp))/(Qm_g_0/rho_g_a(t,Nx_Lp)+Qm_l_0/rho_l_a(t,Nx_Lp)); % 泡沫质量
        gamma_l_a(t,Nx_Lp)=1-gamma_g_a(t,Nx_Lp); % 液体滞留量
        alpha_g_a(t,Nx_Lp)=alpha_f_a(t-1,Nx_Lp)*gamma_g_a(t,Nx_Lp); % 环空气体含量
        alpha_l_a(t,Nx_Lp)=alpha_f_a(t-1,Nx_Lp)*gamma_l_a(t,Nx_Lp); % 环空基液含量
        rho_f_a(t,Nx_Lp)=rho_g_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp)+rho_l_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % 环空泡沫密度，kg/m^3
        mu_f_a(t,Nx_Lp)=mu_g_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp)+mu_l_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % 环空泡沫粘度，Pa*s
        mu_s(t,Nx_Lp)=mu_f_a(t,Nx_Lp); % 固相粘度，Pa*s
        Va_f_a(t,Nx_Lp)=Qm_f_0/(A_a(t,Nx_Lp)*rho_f_a(t,Nx_Lp)); % 环空泡沫表观流速，m/s
        Va_g_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp); % 环空气体表观流速，m/s
        Va_l_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp); % 环空基液表观流速，m/s
        Va_s(t,Nx_Lp)=M_s(t)/(A_a(t,Nx_Lp)*rho_s); % 岩屑表观速度，m/s
        Vsr(t,Nx_Lp)=12*(mu_f_a(t,Nx_Lp)/(rho_f_a(t,Nx_Lp)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Lp))/rho_f_a(t,Nx_Lp))*((rho_f_a(t,Nx_Lp)*D_s/mu_f_a(t,Nx_Lp))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
        alpha_s(t,Nx_Lp)=Va_s(t,Nx_Lp)/(C0*(Va_s(t,Nx_Lp)+Va_f_a(t,Nx_Lp))-Vsr(t,Nx_Lp));  % 固相体积分数
        V_s(t,Nx_Lp)=0;%Va_s(t,Nx_Lp)/alpha_s(t,Nx_Lp); % 岩屑速度，m/s
        
        alpha_f_a(t,Nx_Lp)=1-alpha_s(t,Nx_Lp); % 环空泡沫含量
        alpha_g_a(t,Nx_Lp)=alpha_f_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp); % 环空气体含量
        alpha_l_a(t,Nx_Lp)=alpha_f_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % 环空基液含量
        V_f_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp)/alpha_f_a(t,Nx_Lp); % 环空泡沫流速，m/s
        V_g_a(t,Nx_Lp)=V_f_a(t,Nx_Lp); % 环空气体流速，m/s
        V_l_a(t,Nx_Lp)=V_f_a(t,Nx_Lp); % 环空基液流速，m/s
        
        V_m(t,Nx_Lp)=Va_s(t,Nx_Lp)+Va_f_a(t,Nx_Lp); % 环空混合物速度，m/s
        rho_m(t,Nx_Lp)=alpha_s(t,Nx_Lp)*rho_s+alpha_f_a(t,Nx_Lp)*rho_f_a(t,Nx_Lp); % 环空混合物密度，kg/m^3
        mu_m(t,Nx_Lp)=alpha_s(t,Nx_Lp)*mu_s(t,Nx_Lp)+alpha_f_a(t,Nx_Lp)*mu_f_a(t,Nx_Lp); % 环空混合物粘度，Pa*s
        [Ff_a(t,Nx_Lp),f_a(t,Nx_Lp),Re_a(t,Nx_Lp),flow_pattern_a(t,Nx_Lp)]=Friction_annulus(rho_m(t,Nx_Lp),V_m(t,Nx_Lp),mu_m(t,Nx_Lp),D_h(t,Nx_Lp),h_a,rho_f_a(t,Nx_Lp),V_f_a(t,Nx_Lp)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
        
        % 第Nx_Lp-1～1个空间节点处相关参数计算
        for x=Nx_Lp-1:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % 环空压力假设值，Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液粘度，Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % 环空基液含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            
            err_NodePressure=1; % 环空压力相对误差
            COUNT_NodePressure=0; % 环空压力迭代次数初值
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % 环空压力迭代次数
                
                % 割线法求解固相体积分数
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % 固相体积分数假设值1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % 固相体积分数假设值2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % 固相体积分数绝对误差
                COUNT_NodeEg=0; % 迭代次数
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % 固相体积分数为alpha_s_ass1时（假设值）
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % 固相速度，m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % 固相表观流速，m/s
                    alpha_f_ass1=1-alpha_s_ass1; % 泡沫体积分数
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % 泡沫速度，m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % 泡沫表观流速，m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % 环空混合物速度，m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % 固相体积分数计算值
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 固相体积分数为alpha_s_ass2时（假设值）
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % 固相质量守恒方程离散后公式中间值计算
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % 固相速度，m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % 固相表观流速，m/s
                    alpha_f_ass2=1-alpha_s_ass2; % 泡沫体积分数
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % 泡沫速度，m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % 泡沫表观流速，m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % 环空混合物速度，m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % 固相体积分数计算值
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 割线法计算固相体积分数
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % 新的固相体积分数假设值alpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % 固相体积分数绝对误差
                    alpha_s_ass1=alpha_s_ass2; % 新的固相体积分数假设值1
                    alpha_s_ass2=alpha_s_ass3; % 新的固相体积分数假设值2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % 将迭代求解得到的真实固相体积分数值赋给alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % 当固相体积分数小于一定值时，认为固相体积数为0，用于防止出现后续的计算误差
                    alpha_s(t,x)=0; % 固相体积分数
                    V_s(t,x)=0; % 固相速度，m/s
                    Va_s(t,x)=0; % 固相表观速度，m/s
                    mu_s(t,x)=0; % 固相粘度，Pa*s
                    Vsr(t,x)=0; % 岩屑沉降末速，m/s
                    
                    alpha_f_a(t,x)=1; % 泡沫体积分数                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                else
                    V_s(t,x)=V_s_ass2; % 固相速度，m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % 固相表观流速，m/s
                    Vsr(t,x)=Vsr_ass2; % 岩屑沉降末速，m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % 泡沫密度，kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % 泡沫体积分数
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    V_f_a(t,x)=V_f_ass2; % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % 环空混合物速度，m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % 环空压力计算值，Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % 环空压力相对误差
                P_a(t,x)=P_new; % 新的环空压力假设值，Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % 出口压力相对误差        
        if (P_a(t,1)-OutPressure)>0 % 根据出口压力误差的正负，对环空管底压力假设值进行调节
            P_a(t,Nx_Lp)=P_a(t,Nx_Lp)-(P_a(t,1)-OutPressure)/2; % 新的环空管底压力假设值，Pa
        else
            P_a(t,Nx_Lp)=P_a(t,Nx_Lp)-(P_a(t,1)-OutPressure)/2*0.3; % 新的环空管底压力假设值，Pa
        end
    end
    
    % 第Nx_Lp+1～Nx_Dsp_pen(1)个空间节点处相关参数计算
    for x=Nx_Lp+1:1:Nx_Dsp_pen(1)
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=0; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
    
    % 第Nx_Dsp_pen(1)+1～Nx个空间节点处相关参数计算
    for x=Nx_Dsp_pen(1)+1:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=0; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
end

%% 第1～Nt时间节点钻头压降计算
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Lp); % 井底泡沫体积流量，m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Lp)); % 计算冲砂工具压降（Pa）及射流喷嘴流速V_nozzle（m/s）
end

%% 第1～Nt时间节点相关参数计算（连续管内）
for t=1:1:Nt
    % 管底（第Nx_Lp个空间节点）处相关参数计算
    P_ct(t,Nx_Lp)=P_a(t,Nx_Lp)+delta_P_SWT(t); % 管内压力，Pa
    rho_g_ct(t,Nx_Lp)=DensityG(T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % 管内气体密度，kg/m^3
    rho_l_ct(t,Nx_Lp)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % 管内基液密度，kg/m^3
    mu_g_ct(t,Nx_Lp)=RheologyG(T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % 管内气体粘度，Pa*s
    mu_l_ct(t,Nx_Lp)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % 管内基液粘度，Pa*s
    alpha_f_ct(t,Nx_Lp)=1; % 管内泡沫含量
    gamma_g_ct(t,Nx_Lp)=(Qm_g_0/rho_g_ct(t,Nx_Lp))/(Qm_g_0/rho_g_ct(t,Nx_Lp)+Qm_l_0/rho_l_ct(t,Nx_Lp)); % 管内泡沫质量
    gamma_l_ct(t,Nx_Lp)=1-gamma_g_ct(t,Nx_Lp); % 管内液体滞留量
    alpha_g_ct(t,Nx_Lp)=alpha_f_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp); % 管内气体含量
    alpha_l_ct(t,Nx_Lp)=alpha_f_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % 管内基液含量
    rho_f_ct(t,Nx_Lp)=rho_g_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp)+rho_l_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % 管内泡沫密度，kg/m^3
    mu_f_ct(t,Nx_Lp)=mu_g_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp)+mu_l_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % 管内泡沫粘度，Pa*s
    V_f_ct(t,Nx_Lp)=Qm_f_0/(rho_f_ct(t,Nx_Lp)*A_ct(t,Nx_Lp)); % 管内气体流速，m/s
    [Ff_ct(t,Nx_Lp),f_ct(t,Nx_Lp),Re_ct(t,Nx_Lp),flow_pattern_ct(t,Nx_Lp)]=Friction_coiledtubing(rho_f_ct(t,Nx_Lp),V_f_ct(t,Nx_Lp),mu_f_ct(t,Nx_Lp),D_ct_i(t,Nx_Lp),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    % 第Nx_Lp-1～1个空间节点处相关参数计算
    for x=Nx_Lp-1:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % 管内压力假设值，Pa
        
        err_DriPipePressure=1; % 管内压力相对误差
        COUNT_DriPipePressure=0; % 管内压力迭代次数初值
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % 管内压力迭代次数
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体密度，kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液密度，kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体粘度，Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液粘度，Pa*s
            alpha_f_ct(t,x)=1; % 管内泡沫含量
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % 管内泡沫质量
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % 管内气体含量
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % 管内基液含量
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫密度，kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫粘度，Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % 管内泡沫流速，m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % 管内压力，Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % 计算管内压力假设值与计算值之间的相对误差
            P_ct_ass(t,x)=P_ct(t,x); % 新的管内压力假设值，Pa
        end
    end
    
    % 第Nx_Lp+1～Nx个空间节点处相关参数计算
    for x=Nx_Lp+1:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % 管内气体含量
        alpha_l_ct(t,x)=alpha_l_a(t,x); % 管内基液含量
        gamma_g_ct(t,x)=gamma_g_a(t,x); % 管内泡沫质量
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
        rho_f_ct(t,x)=rho_m(t,x); % 管内气体密度，kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % 管内气体粘度，Pa*s
        V_f_ct(t,t)=V_m(t,x); % 管内气体流速，sm/s
        Ff_ct(t,x)=Ff_a(t,x); % 管内流体单位长度摩擦压降（Pa/m）
        f_ct(t,x)=f_a(t,x); % 管内流体范宁摩擦因子
        Re_ct(t,x)=Re_a(t,x); % 管内流体雷诺数、流体流态
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % 管内流体流态
        P_ct(t,x)=P_a(t,x); % 管内压力，Pa
    end
end

%% 第1～Nt时间节点盘管段出口压力及泵压计算（连续管内）
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % 地面管内流体流速，m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % 地面管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % 盘管段流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % 导向器内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % 井口到注入头顶部段摩擦压降，Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % 导向器段摩擦压降，Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % 导向器到滚筒段摩擦压降，Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % 盘管段摩擦压降，Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % 盘管段出口压力，Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % 泵压，Pa
end

%% 计算井口累积出砂量
M_w_tem(1)=0; % 井口瞬时出砂量初值，kg
M_w_tot(1)=0; % 井口累积出砂量初值，kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口瞬时出砂量，kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口累积出砂量，kg
end

%% 计算末态井筒最大砂浓度
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(1)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% 计算环空泡沫返速平均值V_l_mre
V_f_mre=0; % 环空返速平均值，m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % 环空返速平均值，m/s
    end
end

%% 计算固相沉降末速平均值Vsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % 从Vsr矩阵中提取值非零的沉降末速，m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % 固相沉降末速平均值，m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % 固相沉降末速平均值，m/s
end

%% 判断是否有效冲砂（环空返速平均值大于2倍砂砾沉降末速平均值）、是否完成冲砂
fprintf("定点循环过程2：\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“有效冲砂”
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % 命令行窗口说明当前条件下“完成冲砂”
    else
        fprintf("Sand Cleanout UnFinished!\n"); % 命令行窗口说明当前条件下“未完成冲砂”
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“无效冲砂”
end

%% 井筒ECD计算
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % 环空ECD，kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % 管内ECD，kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % 环空ECD，kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % 管内ECD，kg/m^3
end

%% 定点循环过程末状态数据
alpha_f_a_7=alpha_f_a(Nt,:); % 定点循环过程末状态泡沫含量
alpha_g_a_7=alpha_g_a(Nt,:); % 定点循环过程末状态气相含量
alpha_l_a_7=alpha_l_a(Nt,:); % 定点循环过程末状态液相含量
alpha_s_7=alpha_s(Nt,:); % 定点循环过程末状态固相含量
f_a_7=f_a(Nt,:); % 定点循环过程末状态环空摩擦因子
Ff_a_7=Ff_a(Nt,:); % 定点循环过程末状态环空单位长度摩擦压降，Pa/m
flow_pattern_a_7=flow_pattern_a(Nt,:); % 定点循环过程末状态环空流体流态
gamma_g_a_7=gamma_g_a(Nt,:); % 定点循环过程末状态泡沫质量
gamma_l_a_7=gamma_l_a(Nt,:); % 定点循环过程末状态液体滞留量
mu_f_a_7=mu_f_a(Nt,:); % 定点循环过程末状态环空泡沫粘度，Pa*s
mu_g_a_7=mu_g_a(Nt,:); % 定点循环过程末状态气相粘度，Pa*s
mu_l_a_7=mu_l_a(Nt,:); % 定点循环过程末状态环空液体粘度，Pa*s
mu_m_7=mu_m(Nt,:); % 定点循环过程末状态环空混合物粘度，Pa*s
mu_s_7=mu_s(Nt,:); % 定点循环过程末状态固相粘度，Pa*s
P_a_7=P_a(Nt,:); % 定点循环过程末状态环空压力，Pa
Re_a_7=Re_a(Nt,:); % 定点循环过程末状态环空雷诺数
rho_f_a_7=rho_f_a(Nt,:); % 定点循环过程末状态环空泡沫密度，kg/m^3
rho_g_a_7=rho_g_a(Nt,:); % 定点循环过程末状态环空气相密度，kg/m^3
rho_l_a_7=rho_l_a(Nt,:); % 定点循环过程末状态环空液相密度，kg/m^3
rho_m_7=rho_m(Nt,:); % 定点循环过程末状态环空混合物密度，kg/m^3
V_f_a_7=V_f_a(Nt,:); % 定点循环过程末状态环空泡沫流速，m/s
V_g_a_7=V_g_a(Nt,:); % 定点循环过程末状态气相流速，m/s
V_l_a_7=V_l_a(Nt,:); % 定点循环过程末状态环空液相流速，m/s
V_m_7=V_m(Nt,:); % 定点循环过程末状态环空混合物流速，m/s
V_s_7=V_s(Nt,:); % 定点循环过程末状态固相流速，m/s
Va_f_a_7=Va_f_a(Nt,:); % 定点循环过程末状态环空泡沫表观流速，m/s
Va_g_a_7=Va_g_a(Nt,:); % 定点循环过程末状态气相表观流速，m/s
Va_l_a_7=Va_l_a(Nt,:); % 定点循环过程末状态环空液相表观流速，m/s
Va_s_7=Va_s(Nt,:); % 定点循环过程末状态固相表观流速，m/s
Vsr_7=Vsr(Nt,:); % 定点循环过程末状态固相滑移速度，m/s

%% 数据存储（定点循环）
ANS_Nt_7=Nt; % 时间节点数
ANS_Time_7=Time(1:Nt); % 定点循环时间，s
ANS_alpha_g_a_7=alpha_g_a(1:Nt,:); % 环空气相体积分数
ANS_alpha_l_a_7=alpha_l_a(1:Nt,:); % 环空液相体积分数
ANS_alpha_s_7=alpha_s(1:Nt,:); % 固相含量
ANS_P_ct_7=P_ct(1:Nt,:); % 管内压力（Pa）
ANS_delta_P_SWT_7=delta_P_SWT(1:Nt); % 冲砂工具压降（Pa）
ANS_P_a_7=P_a(1:Nt,:); % 环空压力（Pa）
ANS_T_a_7=T_a(1:Nt,:); % 环空温度（℃）
ANS_T_ct_7=T_ct(1:Nt,:); % 管内温度（℃）
ANS_P_coil_7=P_coil(1:Nt); % 盘管段出口压力（Pa）
ANS_P_pump_7=P_pump(1:Nt); % 泵压（Pa）
ANS_M_w_tem_7=M_w_tem(1:Nt); % 井口瞬时出砂量（kg）
ANS_M_w_tot_7=M_w_tot(1:Nt); % 井口累积出砂量（kg）
ANS_Va_s_7=Va_s(1:Nt,:); % 岩屑表观速度（m/s）
ANS_Va_f_a_7=Va_f_a(1:Nt,:); % 环空泡沫表观流速（m/s）
ANS_V_s_7=V_s(1:Nt,:); % 岩屑沉降速度（m/s）
ANS_V_f_a_7=V_f_a(1:Nt,:); % 环空泡沫返速（m/s）
ANS_alpha_g_ct_7=alpha_g_ct(1:Nt,:); % 管内气相体积分数
ANS_alpha_l_ct_7=alpha_l_ct(1:Nt,:); % 管内液相体积分数
ANS_gamma_g_ct_7=gamma_g_ct(1:Nt,:); % 管内泡沫质量
ANS_gamma_g_a_7=gamma_g_a(1:Nt,:); % 环空泡沫质量
ANS_ECD_a_7=ECD_a(1:Nt,:); % 环空ECD（kg/m^3）
ANS_L_coil_7=L_coil(1:Nt); % 连续管下深（m）



%% 连续管起出过程（国际单位制）
%% 变量数值清零
dt=zeros(); % 时间步长，s
dx=zeros(); % 空间步长，m
Time=zeros(); % 连续管下入总时长，s
L_coil=zeros(); % 连续管下深，m
L_reel=zeros(); % 盘管段长度，m
D_t_i=zeros(); % 油管内径，m
D_ct_o=zeros(); % 连续油管外径，m
D_ct_i=zeros(); % 连续油管内径，m
M_s=zeros(); % 井底进砂量，kg/s
P_a=zeros(); % 环空压力，Pa
rho_g_a=zeros(); % 气相密度，kg/m^3
rho_l_a=zeros(); % 环空液相密度，kg/m^3
rho_f_a=zeros(); % 环空泡沫密度，kg/m^3
mu_g_a=zeros(); % 气相粘度，Pa*s
mu_l_a=zeros(); % 环空液体粘度，Pa*s
mu_f_a=zeros(); % 环空泡沫粘度，Pa*s
mu_s=zeros(); % 固相粘度，Pa*s
alpha_g_a=zeros(); % 气相含量
alpha_l_a=zeros(); % 液相含量
alpha_f_a=zeros(); % 泡沫含量
alpha_s=zeros(); % 固相含量
Va_g_a=zeros(); % 气相表观流速，m/s
Va_l_a=zeros(); % 环空液相表观流速，m/s
Va_f_a=zeros(); % 环空泡沫表观流速，m/s
Va_s=zeros(); % 固相表观流速，m/s
V_g_a=zeros(); % 气相流速，m/s
V_l_a=zeros(); % 环空液相流速，m/s
V_f_a=zeros(); % 环空泡沫流速，m/s
V_s=zeros(); % 固相流速，m/s
Vsr=zeros(); % 固相滑移速度，m/s
gamma_g_a=zeros(); % 泡沫质量
gamma_l_a=zeros(); % 液体滞留量
V_m=zeros(); % 环空混合物流速，m/s
rho_m=zeros(); % 环空混合物密度，kg/m^3
mu_m=zeros(); % 环空混合物粘度，Pa*s
Ff_a=zeros(); % 环空单位长度摩擦压降，Pa/m
f_a=zeros(); % 环空摩擦因子
Re_a=zeros(); % 环空雷诺数
flow_pattern_a=zeros(); % 环空流体流态

%% 计算空间步长dx及相应时间步长dt
V8=Lp/(1.5*3600); % 连续管起出速度，m/s（输入）
t_8=Lp/V8; % 连续管起出总时长，s

Nt=Nx_Lp; % 时间节点数
nt=Nt-1; % 时间网格数

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % 空间步长，m
end

Time(1)=0; % 连续管起出总时长初值，s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Lp-t)/V8; % 每起出一个空间步长所需时间，s
    Time(t+1)=Time(t)+dt(t); % 连续管起出至第(t+1)个空间节点所经历的总时长，s
end

%% 计算不同时刻连续管起出深度L_coil及盘管段连续管长度L_reel
L=10000; % 连续油管总长，m（输入）
L_wg=8; % 井口到注入头顶部段连续管长度，m（输入）
L_goose=3; % 导向器段连续管长度，m（输入）
D_goose=2; % 导向器段半径，m（输入）
H_goose=10; % 导向器距地面高度，m（输入）
L_gr=20; % 导向器到滚筒段连续管长度，m（输入）
theta_gr=acosd(H_goose/L_gr); % 导向器到滚筒段连续管与铅垂线夹角，°
D_r_i=3; % 滚筒内径，m（输入）
D_r_o=5; % 滚筒外径，m（输入）
W_r=5; % 滚筒宽度，m（输入）
D_cable=0.005; % 电缆外径，m

L_coil(1)=Lp; % 初始时刻连续管底部深度，m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % 初始时刻盘管段连续管长度，m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)-dx(Nx_Lp-t+1); % 连续管起出深度，m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % 盘管段连续管长度，m
end
L_cable=L_coil; % 电缆下入长度，m

%% 数据输入及预处理
D_ct_o_0=0.04445; % 连续油管外径，m（输入）
L1=2000; % 近出口第一段连续管长度，m（输入）
D_ct_i_1=0.03709; % 近出口第一段连续管内径，m（输入）
L2=2000; % 近出口第二段连续管长度，m（输入）
D_ct_i_2=0.03653; % 近出口第二段连续管内径，m（输入）
L3=2000; % 近出口第三段连续管长度，m（输入）
D_ct_i_3=0.03555; % 近出口第三段连续管内径，m（输入）
L4=2000; % 近出口第四段连续管长度，m（输入）
D_ct_i_4=0.03489; % 近出口第四段连续管内径，m（输入）
L5=2000; % 近出口第五段连续管长度，m（输入）
D_ct_i_5=0.03409; % 近出口第五段连续管内径，m（输入）

L_t_1=4000; % 上部油管（或套管或裸眼）长度，m（输入）
D_t_i_1=0.09718; %0.068;%0.09718; % 上部油管（或套管或裸眼）内径，m（输入）
L_t_2=2200; % 下部油管（或套管或裸眼）长度，m（输入）
D_t_i_2=0.09718; %0.13970; % 下部油管（或套管或裸眼）内径，m（输入）

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % 油管（或套管或裸眼）内径，m
        else
            D_t_i(t,x)=D_t_i_2; % 油管（或套管或裸眼）内径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % 连续油管外径，m
            else
                D_ct_o(t,x)=0; % 连续油管外径，m
            end
        else
            D_ct_o(t,x)=0; % 连续油管外径，m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            if Depth(x)<=L_coil(t) % 决定哪几个点可以有内径
                if (L_coil(t)-Depth(x))<=L1 && (L_coil(t)-Depth(x))>=0
                    D_ct_i(t,x)=D_ct_i_1; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % 连续油管内径，m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % 连续油管内径，m
                end
            else
                D_ct_i(t,x)=0; % 连续油管内径，m
            end
        else
            D_ct_i(t,x)=0; % 连续油管内径，m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % 决定哪几个点可以有数值
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % 连续油管内截面积，m^2
        else
            A_ct(t,x)=0; % 连续油管内截面积，m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % 环空水力直径，m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % 环空截面积，m^2
    end
end

%% 数据输入及预处理
h_t=25.4*10^(-6); % 油管（或套管或裸眼）绝对粗糙度，m（输入）
h_ct=25.4*10^(-6); % 连续油管绝对粗糙度，m（输入）
h_a=(h_t+h_ct)/2; % 环空平均绝对粗糙度，m
epsilon_e=1*10^(-3); % 迭代求解误差限（输入）
epsilon_t=1*10^3; % 最大迭代次数（输入）
g=9.81; % 重力加速度，m/s^2（默认）

T_0=20; % 基液测试温度，℃（输入）
P_0=0.1*10^6; % 基液测试压力，Pa（输入）
rho_l_0=1150; % T_0、P_0下基液密度，kg/m^3（输入）
mu_l_0=0.03; % T_0、P_0下基液粘度，Pa*s（输入）
Qv_l_0=0.1/60; % 基液体积流量，m^3/s（输入）
Qm_l_0=Qv_l_0*rho_l_0; % 基液质量流量，kg/s
rho_g_0=0.655; % 注入气密度，kg/m^3
mu_g_0=RheologyG(T_0,P_0); % 气体粘度，Pa*s
Qv_g_0=4/60; % 注入气体积流量，m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % 注入气质量流量，kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % 泡沫质量流量，kg/s

D_s=1*10^(-3); % 砂砾直径，m（输入）
rho_s=2000; % 砂砾密度，kg/m^3（输入）
H_s=L_s_b-L_s_t; % 底部砂床高度，m
PHI=0.6; % 砂床充盈度（输入）
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % 井底砂砾总质量，kg

D_nozzle=4/1000; % 喷嘴直径，m（输入）
N_nozzle=3; % 喷嘴个数（输入）
C=0.95; % 喷嘴流量系数，取0.95（输入）

C0=1.2; % 漂移流方程系数（默认）

M_s(1)=0; % 井底进砂量，kg/s
for t=2:1:Nt
    M_s(t)=0; % 井底进砂量，kg/s
end

OutPressure=1*10^6; % 井口压力，Pa（输入）

%% 温度设置
T_i=20; % 泡沫注入温度，℃
T_g=0.02; % 地温梯度，℃/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % 连续管内泡沫温度（假设每个时刻都一样），℃
    end
end
T_a=T_ct; % 环空泡沫温度（假设每个时刻都一样），℃

%% 第1个时间节点（初始时刻）相关参数插值计算（环空）
for x=1:1:Nx
    P_a(1,x)=P_a_7(x); % 环空压力，Pa
    rho_g_a(1,x)=rho_g_a_7(x); % 气相密度，kg/m^3
    rho_l_a(1,x)=rho_l_a_7(x); % 环空液相密度，kg/m^3
    rho_f_a(1,x)=rho_f_a_7(x); % 环空泡沫密度，kg/m^3
    mu_g_a(1,x)=mu_g_a_7(x); % 气相粘度，Pa*s
    mu_l_a(1,x)=mu_l_a_7(x); % 环空液体粘度，Pa*s
    mu_f_a(1,x)=mu_f_a_7(x); % 环空泡沫粘度，Pa*s
    mu_s(1,x)=mu_s_7(x); % 固相粘度，Pa*s
    alpha_g_a(1,x)=alpha_g_a_7(x); % 气相含量
    alpha_l_a(1,x)=alpha_l_a_7(x); % 液相含量
    alpha_f_a(1,x)=alpha_f_a_7(x); % 泡沫含量
    alpha_s(1,x)=alpha_s_7(1,x); % 固相含量
    Va_g_a(1,x)=Va_g_a_7(x); % 气相表观流速，m/s
    Va_l_a(1,x)=Va_l_a_7(x); % 环空液相表观流速，m/s
    Va_f_a(1,x)=Va_f_a_7(x); % 环空泡沫表观流速，m/s
    Va_s(1,x)=Va_s_7(x); % 固相表观流速，m/s
    V_g_a(1,x)=V_g_a_7(x); % 气相流速，m/s
    V_l_a(1,x)=V_l_a_7(x); % 环空液相流速，m/s
    V_f_a(1,x)=V_f_a_7(x); % 环空泡沫流速，m/s
    V_s(1,x)=V_s_7(x); % 固相流速，m/s
    Vsr(1,x)=Vsr_7(x); % 固相滑移速度，m/s
    gamma_g_a(1,x)=gamma_g_a_7(x); % 泡沫质量
    gamma_l_a(1,x)=gamma_l_a_7(1,x); % 液体滞留量
    V_m(1,x)=V_m_7(x); % 环空混合物流速，m/s
    rho_m(1,x)=rho_m_7(x); % 环空混合物密度，kg/m^3
    mu_m(1,x)=mu_m_7(x); % 环空混合物粘度，Pa*s
    Ff_a(1,x)=Ff_a_7(x); % 环空单位长度摩擦压降，Pa/m
    f_a(1,x)=f_a_7(x); % 环空摩擦因子
    Re_a(1,x)=Re_a_7(x); % 环空雷诺数
    flow_pattern_a(1,x)=flow_pattern_a_7(x); % 环空流体流态
end

%% 第2～Nt个时间节点相关参数计算（环空）
for t=2:1:Nt
    P_a(t,Nx_Lp-t+1)=P_a(t-1,Nx_Lp-t+1);  % 环空管底压力假设值，Pa
    
    err_OutPressure=1; % 出口压力相对误差
    COUNT_OutPressure=0; % 出口压力迭代次数初值
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % 出口压力迭代次数
        
        % 管底（第Nx_Lp-t+1个空间节点）处相关参数计算
        rho_g_a(t,Nx_Lp-t+1)=DensityG(T_a(t,Nx_Lp-t+1),P_a(t,Nx_Lp-t+1)); % 环空气体密度，kg/m^3
        rho_l_a(t,Nx_Lp-t+1)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Lp-t+1),P_a(t,Nx_Lp-t+1)); % 环空基液密度，kg/m^3
        mu_g_a(t,Nx_Lp-t+1)=RheologyG(T_a(t,Nx_Lp-t+1),P_a(t,Nx_Lp-t+1)); % 环空气体粘度，Pa*s
        mu_l_a(t,Nx_Lp-t+1)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Lp-t+1),P_a(t,Nx_Lp-t+1)); % 环空基液粘度，Pa*s
        gamma_g_a(t,Nx_Lp-t+1)=(Qm_g_0/rho_g_a(t,Nx_Lp-t+1))/(Qm_g_0/rho_g_a(t,Nx_Lp-t+1)+Qm_l_0/rho_l_a(t,Nx_Lp-t+1)); % 泡沫质量
        gamma_l_a(t,Nx_Lp-t+1)=1-gamma_g_a(t,Nx_Lp-t+1); % 液体滞留量
        alpha_g_a(t,Nx_Lp-t+1)=alpha_f_a(t-1,Nx_Lp-t+1)*gamma_g_a(t,Nx_Lp-t+1); % 环空气体含量
        alpha_l_a(t,Nx_Lp-t+1)=alpha_f_a(t-1,Nx_Lp-t+1)*gamma_l_a(t,Nx_Lp-t+1); % 环空基液含量
        rho_f_a(t,Nx_Lp-t+1)=rho_g_a(t,Nx_Lp-t+1)*gamma_g_a(t,Nx_Lp-t+1)+rho_l_a(t,Nx_Lp-t+1)*gamma_l_a(t,Nx_Lp-t+1); % 环空泡沫密度，kg/m^3
        mu_f_a(t,Nx_Lp-t+1)=mu_g_a(t,Nx_Lp-t+1)*gamma_g_a(t,Nx_Lp-t+1)+mu_l_a(t,Nx_Lp-t+1)*gamma_l_a(t,Nx_Lp-t+1); % 环空泡沫粘度，Pa*s
        mu_s(t,Nx_Lp-t+1)=mu_f_a(t,Nx_Lp-t+1); % 固相粘度，Pa*s
        Va_f_a(t,Nx_Lp-t+1)=Qm_f_0/(A_a(t,Nx_Lp-t+1)*rho_f_a(t,Nx_Lp-t+1)); % 环空泡沫表观流速，m/s
        Va_g_a(t,Nx_Lp-t+1)=Va_f_a(t,Nx_Lp-t+1); % 环空气体表观流速，m/s
        Va_l_a(t,Nx_Lp-t+1)=Va_f_a(t,Nx_Lp-t+1); % 环空基液表观流速，m/s
        Va_s(t,Nx_Lp-t+1)=M_s(t)/(A_a(t,Nx_Lp-t+1)*rho_s); % 岩屑表观速度，m/s
        Vsr(t,Nx_Lp-t+1)=12*(mu_f_a(t,Nx_Lp-t+1)/(rho_f_a(t,Nx_Lp-t+1)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Lp-t+1))/rho_f_a(t,Nx_Lp-t+1))*((rho_f_a(t,Nx_Lp-t+1)*D_s/mu_f_a(t,Nx_Lp-t+1))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
        alpha_s(t,Nx_Lp-t+1)=Va_s(t,Nx_Lp-t+1)/(C0*(Va_s(t,Nx_Lp-t+1)+Va_f_a(t,Nx_Lp-t+1))-Vsr(t,Nx_Lp-t+1));  % 固相体积分数
        V_s(t,Nx_Lp-t+1)=0;%Va_s(t,Nx_Lp-t+1)/alpha_s(t,Nx_Lp-t+1); % 岩屑速度，m/s
        
        alpha_f_a(t,Nx_Lp-t+1)=1-alpha_s(t,Nx_Lp-t+1); % 环空泡沫含量
        alpha_g_a(t,Nx_Lp-t+1)=alpha_f_a(t,Nx_Lp-t+1)*gamma_g_a(t,Nx_Lp-t+1); % 环空气体含量
        alpha_l_a(t,Nx_Lp-t+1)=alpha_f_a(t,Nx_Lp-t+1)*gamma_l_a(t,Nx_Lp-t+1); % 环空基液含量
        V_f_a(t,Nx_Lp-t+1)=Va_f_a(t,Nx_Lp-t+1)/alpha_f_a(t,Nx_Lp-t+1); % 环空泡沫流速，m/s
        V_g_a(t,Nx_Lp-t+1)=V_f_a(t,Nx_Lp-t+1); % 环空气体流速，m/s
        V_l_a(t,Nx_Lp-t+1)=V_f_a(t,Nx_Lp-t+1); % 环空基液流速，m/s
        
        V_m(t,Nx_Lp-t+1)=Va_s(t,Nx_Lp-t+1)+Va_f_a(t,Nx_Lp-t+1); % 环空混合物速度，m/s
        rho_m(t,Nx_Lp-t+1)=alpha_s(t,Nx_Lp-t+1)*rho_s+alpha_f_a(t,Nx_Lp-t+1)*rho_f_a(t,Nx_Lp-t+1); % 环空混合物密度，kg/m^3
        mu_m(t,Nx_Lp-t+1)=alpha_s(t,Nx_Lp-t+1)*mu_s(t,Nx_Lp-t+1)+alpha_f_a(t,Nx_Lp-t+1)*mu_f_a(t,Nx_Lp-t+1); % 环空混合物粘度，Pa*s
        [Ff_a(t,Nx_Lp-t+1),f_a(t,Nx_Lp-t+1),Re_a(t,Nx_Lp-t+1),flow_pattern_a(t,Nx_Lp-t+1)]=Friction_annulus(rho_m(t,Nx_Lp-t+1),V_m(t,Nx_Lp-t+1),mu_m(t,Nx_Lp-t+1),D_h(t,Nx_Lp-t+1),h_a,rho_f_a(t,Nx_Lp-t+1),V_f_a(t,Nx_Lp-t+1)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
        
        % 第Nx_Lp-t～1个空间节点处相关参数计算
        for x=Nx_Lp-t:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % 环空压力假设值，Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % 环空基液粘度，Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % 环空基液含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            
            err_NodePressure=1; % 环空压力相对误差
            COUNT_NodePressure=0; % 环空压力迭代次数初值
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % 环空压力迭代次数
                
                % 割线法求解固相体积分数
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % 固相体积分数假设值1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % 固相体积分数假设值2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % 固相体积分数绝对误差
                COUNT_NodeEg=0; % 迭代次数
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % 固相体积分数为alpha_s_ass1时（假设值）
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % 固相速度，m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % 固相表观流速，m/s
                    alpha_f_ass1=1-alpha_s_ass1; % 泡沫体积分数
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % 泡沫速度，m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % 泡沫表观流速，m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % 环空混合物速度，m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % 固相体积分数计算值
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 固相体积分数为alpha_s_ass2时（假设值）
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % 固相质量守恒方程离散后公式中间值计算
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % 固相速度，m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % 固相表观流速，m/s
                    alpha_f_ass2=1-alpha_s_ass2; % 泡沫体积分数
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 泡沫密度，kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % 泡沫速度，m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % 泡沫表观流速，m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % 环空混合物速度，m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % 岩屑沉降末速，m/s（Chien）
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % 固相体积分数计算值
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % 构造的函数，它的解就是真实固相体积分数
                    
                    % 割线法计算固相体积分数
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % 新的固相体积分数假设值alpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % 固相体积分数绝对误差
                    alpha_s_ass1=alpha_s_ass2; % 新的固相体积分数假设值1
                    alpha_s_ass2=alpha_s_ass3; % 新的固相体积分数假设值2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % 将迭代求解得到的真实固相体积分数值赋给alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % 当固相体积分数小于一定值时，认为固相体积数为0，用于防止出现后续的计算误差
                    alpha_s(t,x)=0; % 固相体积分数
                    V_s(t,x)=0; % 固相速度，m/s
                    Va_s(t,x)=0; % 固相表观速度，m/s
                    mu_s(t,x)=0; % 固相粘度，Pa*s
                    Vsr(t,x)=0; % 岩屑沉降末速，m/s
                    
                    alpha_f_a(t,x)=1; % 泡沫体积分数                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                else
                    V_s(t,x)=V_s_ass2; % 固相速度，m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % 固相表观流速，m/s
                    Vsr(t,x)=Vsr_ass2; % 岩屑沉降末速，m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % 泡沫密度，kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % 泡沫体积分数
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
                    V_f_a(t,x)=V_f_ass2; % 泡沫速度，m/s
                    V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
                    V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % 泡沫表观流速，m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % 环空混合物速度，m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % 环空压力计算值，Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % 环空压力相对误差
                P_a(t,x)=P_new; % 新的环空压力假设值，Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % 出口压力相对误差
        if (P_a(t,1)-OutPressure)>0 % 根据出口压力误差的正负，对环空管底压力假设值进行调节
            P_a(t,Nx_Lp-t+1)=P_a(t,Nx_Lp-t+1)-(P_a(t,1)-OutPressure)/2; % 新的环空管底压力假设值，Pa
        else
            P_a(t,Nx_Lp-t+1)=P_a(t,Nx_Lp-t+1)-(P_a(t,1)-OutPressure)/2*0.3; % 新的环空管底压力假设值，Pa
        end
    end
    
    % 第Nx_Lp-t+2～Nx个空间节点处相关参数计算
    for x=Nx_Lp-t+2:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % 环空压力假设值，Pa
        
        err_AnnPressure=1; % 环空压力相对误差
        COUNT_AnnPressure=0; % 环空压力迭代次数初值
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % 环空压力迭代次数
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % 环空气体密度，kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液密度，kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % 环空气体粘度，Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % 环空基液粘度，Pa*s
            alpha_f_a(t,x)=1; % 环空泡沫含量
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % 泡沫质量
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % 液体滞留量
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % 环空气体含量
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % 环空基液含量
            alpha_s(t,x)=0; % 固相体积含量
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫密度，kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % 环空泡沫粘度，Pa*s
            V_s(t,x)=0; % 固相速度，m/s
            Va_s(t,x)=0; % 固相表观流速，m/s
            mu_s(t,x)=mu_f_a(t,x); % 固相粘度，Pa*s
            Vsr(t,x)=0; % 砂砾沉降末速，m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % 环空泡沫表观流速，m/s
            Va_g_a(t,x)=Va_f_a(t,x); % 环空气体表观流速，m/s
            Va_l_a(t,x)=Va_f_a(t,x); % 环空基液表观流速，m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % 环空泡沫流速，m/s
            V_g_a(t,x)=V_f_a(t,x); % 环空气体流速，m/s
            V_l_a(t,x)=V_f_a(t,x); % 环空基液流速，m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % 环空混合物速度，m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % 环空混合物密度，kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % 环空混合物粘度，Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % 环空流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % 环空压力，Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % 计算环空压力假设值与计算值之间的相对误差
            P_a_ass(t,x)=P_a(t,x); % 新的环空压力假设值，Pa
        end
    end
end

%% 第1～Nt时间节点钻头压降计算
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Lp-t+1); % 井底泡沫体积流量，m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Lp-t+1)); % 计算冲砂工具压降（Pa）及射流喷嘴流速V_nozzle（m/s）
end

%% 第1～Nt时间节点相关参数计算（连续管内）
for t=1:1:Nt
    % 管底（第Nx_Lp-t+1个空间节点）处相关参数计算
    P_ct(t,Nx_Lp-t+1)=P_a(t,Nx_Lp-t+1)+delta_P_SWT(t); % 管内压力，Pa
    rho_g_ct(t,Nx_Lp-t+1)=DensityG(T_ct(t,Nx_Lp-t+1),P_ct(t,Nx_Lp-t+1)); % 管内气体密度，kg/m^3
    rho_l_ct(t,Nx_Lp-t+1)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Lp-t+1),P_ct(t,Nx_Lp-t+1)); % 管内基液密度，kg/m^3
    mu_g_ct(t,Nx_Lp-t+1)=RheologyG(T_ct(t,Nx_Lp-t+1),P_ct(t,Nx_Lp-t+1)); % 管内气体粘度，Pa*s
    mu_l_ct(t,Nx_Lp-t+1)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Lp-t+1),P_ct(t,Nx_Lp-t+1)); % 管内基液粘度，Pa*s
    alpha_f_ct(t,Nx_Lp-t+1)=1; % 管内泡沫含量
    gamma_g_ct(t,Nx_Lp-t+1)=(Qm_g_0/rho_g_ct(t,Nx_Lp-t+1))/(Qm_g_0/rho_g_ct(t,Nx_Lp-t+1)+Qm_l_0/rho_l_ct(t,Nx_Lp-t+1)); % 管内泡沫质量
    gamma_l_ct(t,Nx_Lp-t+1)=1-gamma_g_ct(t,Nx_Lp-t+1); % 管内液体滞留量
    alpha_g_ct(t,Nx_Lp-t+1)=alpha_f_ct(t,Nx_Lp-t+1)*gamma_g_ct(t,Nx_Lp-t+1); % 管内气体含量
    alpha_l_ct(t,Nx_Lp-t+1)=alpha_f_ct(t,Nx_Lp-t+1)*gamma_l_ct(t,Nx_Lp-t+1); % 管内基液含量
    rho_f_ct(t,Nx_Lp-t+1)=rho_g_ct(t,Nx_Lp-t+1)*gamma_g_ct(t,Nx_Lp-t+1)+rho_l_ct(t,Nx_Lp-t+1)*gamma_l_ct(t,Nx_Lp-t+1); % 管内泡沫密度，kg/m^3
    mu_f_ct(t,Nx_Lp-t+1)=mu_g_ct(t,Nx_Lp-t+1)*gamma_g_ct(t,Nx_Lp-t+1)+mu_l_ct(t,Nx_Lp-t+1)*gamma_l_ct(t,Nx_Lp-t+1); % 管内泡沫粘度，Pa*s
    V_f_ct(t,Nx_Lp-t+1)=Qm_f_0/(rho_f_ct(t,Nx_Lp-t+1)*A_ct(t,Nx_Lp-t+1)); % 管内气体流速，m/s
    [Ff_ct(t,Nx_Lp-t+1),f_ct(t,Nx_Lp-t+1),Re_ct(t,Nx_Lp-t+1),flow_pattern_ct(t,Nx_Lp-t+1)]=Friction_coiledtubing(rho_f_ct(t,Nx_Lp-t+1),V_f_ct(t,Nx_Lp-t+1),mu_f_ct(t,Nx_Lp-t+1),D_ct_i(t,Nx_Lp-t+1),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    % 第Nx_Lp-t～1个空间节点处相关参数计算
    for x=Nx_Lp-t:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % 管内压力假设值，Pa
        
        err_DriPipePressure=1; % 管内压力相对误差
        COUNT_DriPipePressure=0; % 管内压力迭代次数初值
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % 管内压力迭代次数
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体密度，kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液密度，kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % 管内气体粘度，Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % 管内基液粘度，Pa*s
            alpha_f_ct(t,x)=1; % 管内泡沫含量
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % 管内泡沫质量
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % 管内气体含量
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % 管内基液含量
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫密度，kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % 管内泡沫粘度，Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % 管内泡沫流速，m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % 管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % 管内压力，Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % 计算管内压力假设值与计算值之间的相对误差
            P_ct_ass(t,x)=P_ct(t,x); % 新的管内压力假设值，Pa
        end
    end
    
    % 第Nx_Lp-t+2～Nx个空间节点处相关参数计算
    for x=Nx_Lp-t+2:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % 管内气体含量
        alpha_l_ct(t,x)=alpha_l_a(t,x); % 管内基液含量
        gamma_g_ct(t,x)=gamma_g_a(t,x); % 管内泡沫质量
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % 管内液体滞留量
        rho_f_ct(t,x)=rho_m(t,x); % 管内气体密度，kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % 管内气体粘度，Pa*s
        V_f_ct(t,t)=V_m(t,x); % 管内气体流速，sm/s
        Ff_ct(t,x)=Ff_a(t,x); % 管内流体单位长度摩擦压降（Pa/m）
        f_ct(t,x)=f_a(t,x); % 管内流体范宁摩擦因子
        Re_ct(t,x)=Re_a(t,x); % 管内流体雷诺数、流体流态
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % 管内流体流态
        P_ct(t,x)=P_a(t,x); % 管内压力，Pa
    end
end

%% 第1～Nt时间节点盘管段出口压力及泵压计算（连续管内）
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % 地面管内流体流速，m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % 地面管内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % 盘管段流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % 导向器内流体单位长度摩擦压降（Pa/m）、范宁摩擦因子、雷诺数、流体流态
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % 井口到注入头顶部段摩擦压降，Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % 导向器段摩擦压降，Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % 导向器到滚筒段摩擦压降，Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % 盘管段摩擦压降，Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % 盘管段出口压力，Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % 泵压，Pa
end

%% 计算井口累积出砂量
M_w_tem(1)=0; % 井口瞬时出砂量初值，kg
M_w_tot(1)=0; % 井口累积出砂量初值，kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口瞬时出砂量，kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % 第t个时间节点井口累积出砂量，kg
end

%% 计算末态井筒最大砂浓度
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(2)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% 计算环空泡沫返速平均值V_l_mre
V_f_mre=0; % 环空返速平均值，m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % 环空返速平均值，m/s
    end
end

%% 计算固相沉降末速平均值Vsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % 从Vsr矩阵中提取值非零的沉降末速，m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % 固相沉降末速平均值，m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % 固相沉降末速平均值，m/s
end

%% 判断是否有效冲砂（环空返速平均值大于2倍砂砾沉降末速平均值）、是否完成冲砂
fprintf("连续管起出过程：\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“有效冲砂”
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % 命令行窗口说明当前条件下“完成冲砂”
    else
        fprintf("Sand Cleanout UnFinished!\n"); % 命令行窗口说明当前条件下“未完成冲砂”
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % 命令行窗口说明当前条件下为“无效冲砂”
end

%% 井筒ECD计算
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % 环空ECD，kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % 管内ECD，kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % 环空ECD，kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % 管内ECD，kg/m^3
end

%% 数据存储（连续管起出）
ANS_Nt_8=Nt; % 时间节点数
ANS_Time_8=Time(1:Nt); % 连续管起出时间，s
ANS_alpha_g_a_8=alpha_g_a(1:Nt,:); % 环空气相体积分数
ANS_alpha_l_a_8=alpha_l_a(1:Nt,:); % 环空液相体积分数
ANS_alpha_s_8=alpha_s(1:Nt,:); % 固相含量
ANS_P_ct_8=P_ct(1:Nt,:); % 管内压力（Pa）
ANS_delta_P_SWT_8=delta_P_SWT(1:Nt); % 冲砂工具压降（Pa）
ANS_P_a_8=P_a(1:Nt,:); % 环空压力（Pa）
ANS_T_a_8=T_a(1:Nt,:); % 环空温度（℃）
ANS_T_ct_8=T_ct(1:Nt,:); % 管内温度（℃）
ANS_P_coil_8=P_coil(1:Nt); % 盘管段出口压力（Pa）
ANS_P_pump_8=P_pump(1:Nt); % 泵压（Pa）
ANS_M_w_tem_8=M_w_tem(1:Nt); % 井口瞬时出砂量（kg）
ANS_M_w_tot_8=M_w_tot(1:Nt); % 井口累积出砂量（kg）
ANS_Va_s_8=Va_s(1:Nt,:); % 岩屑表观速度（m/s）
ANS_Va_f_a_8=Va_f_a(1:Nt,:); % 环空泡沫表观流速（m/s）
ANS_V_s_8=V_s(1:Nt,:); % 岩屑沉降速度（m/s）
ANS_V_f_a_8=V_f_a(1:Nt,:); % 环空泡沫返速（m/s）
ANS_alpha_g_ct_8=alpha_g_ct(1:Nt,:); % 管内气相体积分数
ANS_alpha_l_ct_8=alpha_l_ct(1:Nt,:); % 管内液相体积分数
ANS_gamma_g_ct_8=gamma_g_ct(1:Nt,:); % 管内泡沫质量
ANS_gamma_g_a_8=gamma_g_a(1:Nt,:); % 环空泡沫质量
ANS_ECD_a_8=ECD_a(1:Nt,:); % 环空ECD（kg/m^3）
ANS_L_coil_8=L_coil(1:Nt); % 连续管下深（m）



%% 数据
Nt=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7+ANS_Nt_8; % 时间节点数
for t=1:1:ANS_Nt_1
    Time(t)=ANS_Time_1(t); % 总时间，s
    M_w_tot(t)=ANS_M_w_tot_1(t); % 井口累积出砂量（kg）
end
for t=ANS_Nt_1+1:1:ANS_Nt_1+ANS_Nt_2
    Time(t)=Time(ANS_Nt_1)+ANS_Time_2(t-ANS_Nt_1); % 总时间，s
    M_w_tot(t)=M_w_tot(ANS_Nt_1)+ANS_M_w_tot_2(t-ANS_Nt_1); % 井口累积出砂量（kg）
end
for t=ANS_Nt_1+ANS_Nt_2+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2)+ANS_Time_3(t-ANS_Nt_1-ANS_Nt_2); % 总时间，s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2)+ANS_M_w_tot_3(t-ANS_Nt_1-ANS_Nt_2); % 井口累积出砂量（kg）
end

for t=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3)+ANS_Time_4(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3)); % 总时间，s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3)+ANS_M_w_tot_4(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3)); % 井口累积出砂量（kg）
end
for t=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4)+ANS_Time_5(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4)); % 总时间，s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4)+ANS_M_w_tot_5(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4)); % 井口累积出砂量（kg）
end
for t=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5)+ANS_Time_6(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5)); % 总时间，s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5)+ANS_M_w_tot_6(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5)); % 井口累积出砂量（kg）
end
for t=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6)+ANS_Time_7(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6)); % 总时间，s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6)+ANS_M_w_tot_7(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6)); % 井口累积出砂量（kg）
end
for t=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7+ANS_Nt_8
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7)+ANS_Time_8(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7)); % 总时间，s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7)+ANS_M_w_tot_8(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7)); % 井口累积出砂量（kg）
end


alpha_g_a=[ANS_alpha_g_a_1;ANS_alpha_g_a_2;ANS_alpha_g_a_3;ANS_alpha_g_a_4;ANS_alpha_g_a_5;ANS_alpha_g_a_6;ANS_alpha_g_a_7;ANS_alpha_g_a_8]; % 环空气相体积分数
alpha_l_a=[ANS_alpha_l_a_1;ANS_alpha_l_a_2;ANS_alpha_l_a_3;ANS_alpha_l_a_4;ANS_alpha_l_a_5;ANS_alpha_l_a_6;ANS_alpha_l_a_7;ANS_alpha_l_a_8]; % 环空液相体积分数
alpha_s=[ANS_alpha_s_1;ANS_alpha_s_2;ANS_alpha_s_3;ANS_alpha_s_4;ANS_alpha_s_5;ANS_alpha_s_6;ANS_alpha_s_7;ANS_alpha_s_8]; % 固相含量
P_ct=[ANS_P_ct_1;ANS_P_ct_2;ANS_P_ct_3;ANS_P_ct_4;ANS_P_ct_5;ANS_P_ct_6;ANS_P_ct_7;ANS_P_ct_8]; % 管内压力（Pa）
delta_P_SWT=[ANS_delta_P_SWT_1,ANS_delta_P_SWT_2,ANS_delta_P_SWT_3,ANS_delta_P_SWT_4,ANS_delta_P_SWT_5,ANS_delta_P_SWT_6,ANS_delta_P_SWT_7,ANS_delta_P_SWT_8]; % 冲砂工具压降（Pa）
P_a=[ANS_P_a_1;ANS_P_a_2;ANS_P_a_3;ANS_P_a_4;ANS_P_a_5;ANS_P_a_6;ANS_P_a_7;ANS_P_a_8]; % 环空压力（Pa）
T_a=[ANS_T_a_1;ANS_T_a_2;ANS_T_a_3;ANS_T_a_4;ANS_T_a_5;ANS_T_a_6;ANS_T_a_7;ANS_T_a_8]; % 环空温度（℃）
T_ct=[ANS_T_ct_1;ANS_T_ct_2;ANS_T_ct_3;ANS_T_ct_4;ANS_T_ct_5;ANS_T_ct_6;ANS_T_ct_7;ANS_T_ct_8]; % 管内温度（℃）
P_coil=[ANS_P_coil_1,ANS_P_coil_2,ANS_P_coil_3,ANS_P_coil_4,ANS_P_coil_5,ANS_P_coil_6,ANS_P_coil_7,ANS_P_coil_8]; % 盘管段出口压力（Pa）
P_pump=[ANS_P_pump_1,ANS_P_pump_2,ANS_P_pump_3,ANS_P_pump_4,ANS_P_pump_5,ANS_P_pump_6,ANS_P_pump_7,ANS_P_pump_8]; % 泵压（Pa）
M_w_tem=[ANS_M_w_tem_1,ANS_M_w_tem_2,ANS_M_w_tem_3,ANS_M_w_tem_4,ANS_M_w_tem_5,ANS_M_w_tem_6,ANS_M_w_tem_7,ANS_M_w_tem_8]; % 井口瞬时出砂量（kg）
Va_s=[ANS_Va_s_1;ANS_Va_s_2;ANS_Va_s_3;ANS_Va_s_4;ANS_Va_s_5;ANS_Va_s_6;ANS_Va_s_7;ANS_Va_s_8]; % 岩屑表观速度（m/s）
Va_f_a=[ANS_Va_f_a_1;ANS_Va_f_a_2;ANS_Va_f_a_3;ANS_Va_f_a_4;ANS_Va_f_a_5;ANS_Va_f_a_6;ANS_Va_f_a_7;ANS_Va_f_a_8]; % 环空泡沫表观流速（m/s）
V_s=[ANS_V_s_1;ANS_V_s_2;ANS_V_s_3;ANS_V_s_4;ANS_V_s_5;ANS_V_s_6;ANS_V_s_7;ANS_V_s_8]; % 岩屑沉降速度（m/s）
V_f_a=[ANS_V_f_a_1;ANS_V_f_a_2;ANS_V_f_a_3;ANS_V_f_a_4;ANS_V_f_a_5;ANS_V_f_a_6;ANS_V_f_a_7;ANS_V_f_a_8]; % 环空泡沫返速（m/s）
alpha_g_ct=[ANS_alpha_g_ct_1;ANS_alpha_g_ct_2;ANS_alpha_g_ct_3;ANS_alpha_g_ct_4;ANS_alpha_g_ct_5;ANS_alpha_g_ct_6;ANS_alpha_g_ct_7;ANS_alpha_g_ct_8]; % 管内气相体积分数
alpha_l_ct=[ANS_alpha_l_ct_1;ANS_alpha_l_ct_2;ANS_alpha_l_ct_3;ANS_alpha_l_ct_4;ANS_alpha_l_ct_5;ANS_alpha_l_ct_6;ANS_alpha_l_ct_7;ANS_alpha_l_ct_8]; % 管内液相体积分数
gamma_g_ct=[ANS_gamma_g_ct_1;ANS_gamma_g_ct_2;ANS_gamma_g_ct_3;ANS_gamma_g_ct_4;ANS_gamma_g_ct_5;ANS_gamma_g_ct_6;ANS_gamma_g_ct_7;ANS_gamma_g_ct_8]; % 管内泡沫质量
gamma_g_a=[ANS_gamma_g_a_1;ANS_gamma_g_a_2;ANS_gamma_g_a_3;ANS_gamma_g_a_4;ANS_gamma_g_a_5;ANS_gamma_g_a_6;ANS_gamma_g_a_7;ANS_gamma_g_a_8]; % 环空泡沫质量
ECD_a=[ANS_ECD_a_1;ANS_ECD_a_2;ANS_ECD_a_3;ANS_ECD_a_4;ANS_ECD_a_5;ANS_ECD_a_6;ANS_ECD_a_7;ANS_ECD_a_8]; % 环空ECD（kg/m^3）
L_coil=[ANS_L_coil_1,ANS_L_coil_2,ANS_L_coil_3,ANS_L_coil_4,ANS_L_coil_5,ANS_L_coil_6,ANS_L_coil_7,ANS_L_coil_8]; % 连续管下深（m）

%% 绘图
%%
fig1 = figure; % 井筒压力VS井深VS时间
for t = 1:1:Nt
    plot(P_a(t,:)/10^6,Depth,P_ct(t,:)/10^6,Depth,[0,60],[L_coil(t),L_coil(t)],'--','LineWidth',2);
    legend('环空压力','管内压力','连续管下深','FontName','黑体','Location','Best');
    xlabel('压力（MPa）','FontName','黑体','FontSize',12);
    ylabel('井深（m）','FontName','黑体','FontSize',12);
    set(gca,'fontsize',12);
    set(gca,'YDir','reverse');
    box on; % 显示坐标轴的边框
    grid on; % 显示坐标轴的主网格线
    grid minor; % 显示坐标轴的次网格线
    
    if t==1
        title({['作业状态：连续管下至砂床顶'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['作业状态：冲洗钻进1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['作业状态：回拖短起1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['作业状态：定点循环1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['作业状态：连续管下入1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['作业状态：冲洗钻进2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['作业状态：回拖短起2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['作业状态：定点循环2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    else
        title({['作业状态：连续管起出'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    end
    
    frame1 = getframe(fig1);
    im{t} = frame2im(frame1);
end

filename1 = '井筒压力VS井深VS时间.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename1,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename1,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig2 = figure; % 泡沫砂砾速度VS井深VS时间
for t = 1:1:Nt
    subplot(1,2,1);
    plot(V_s(t,:),Depth); % 岩屑沉降速度（m/s）
    xlabel('岩屑沉降速度（m/s）','FontName','黑体','FontSize',10);
    ylabel('井深（m）','FontName','黑体','FontSize',10);
    set(gca,'fontsize',12);
    set(gca,'YDir','reverse');
    box on; % 显示坐标轴的边框
    grid on; % 显示坐标轴的主网格线
    grid minor; % 显示坐标轴的次网格线
    
    if t==1
        title({['作业状态：连续管下至砂床顶'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['作业状态：冲洗钻进1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['作业状态：回拖短起1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['作业状态：定点循环1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['作业状态：连续管下入1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['作业状态：冲洗钻进2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['作业状态：回拖短起2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['作业状态：定点循环2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    else
        title({['作业状态：连续管起出'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    end
    
    subplot(1,2,2);
    plot(V_f_a(t,:),Depth); % 环空泡沫返速（m/s）
    xlabel('环空泡沫返速（m/s）','FontName','黑体','FontSize',10);
    ylabel('井深（m）','FontName','黑体','FontSize',10);
    set(gca,'fontsize',12);
    set(gca,'YDir','reverse');
    box on; % 显示坐标轴的边框
    grid on; % 显示坐标轴的主网格线
    grid minor; % 显示坐标轴的次网格线
    
    frame2 = getframe(fig2);
    im{t} = frame2im(frame2);
end

filename2 = '泡沫砂砾速度VS井深VS时间.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename2,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename2,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig3 = figure; % 砂浓度VS井深VS时间
for t = 1:1:Nt
    alpha_s1(t,:)=alpha_s(t,:);
    ind=find(alpha_s(t,:)==0);
    alpha_s1(t,ind)=NaN;
    plot(Depth,alpha_s1(t,:),[L_coil(t),L_coil(t)],[0,0.6],'--','LineWidth',2);
    legend('砂浓度','连续管下深','FontName','黑体','Location','northwest');
    
    ylabel('砂浓度','FontName','黑体','FontSize',12);
    xlabel('井深（m）','FontName','黑体','FontSize',12);
    set(gca,'fontsize',12);
    box on; % 显示坐标轴的边框
    grid on; % 显示坐标轴的主网格线
    grid minor; % 显示坐标轴的次网格线
    axis([0,6200,0,0.6]);
    
    if t==1
        title({['作业状态：连续管下至砂床顶'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['作业状态：冲洗钻进1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['作业状态：回拖短起1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['作业状态：定点循环1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['作业状态：连续管下入1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['作业状态：冲洗钻进2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['作业状态：回拖短起2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['作业状态：定点循环2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    else
        title({['作业状态：连续管起出'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    end
    
    frame3 = getframe(fig3);
    im{t} = frame2im(frame3);
end

filename3 = '砂浓度VS井深VS时间.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename3,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename3,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig4 = figure; % 泵压VS井深VS时间
for t = 1:1:Nt
    plot(Time(t)/60,P_pump(t)/10^6,'-o');
    hold on;
    xlabel('时间（min）','FontName','黑体','FontSize',12);
    ylabel('泵压（MPa）','FontName','黑体','FontSize',12);
    set(gca,'fontsize',12);
    box on; % 显示坐标轴的边框
    grid on; % 显示坐标轴的主网格线
    grid minor; % 显示坐标轴的次网格线
    
    if t==1
        title({['作业状态：连续管下至砂床顶'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['作业状态：冲洗钻进1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['作业状态：回拖短起1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['作业状态：定点循环1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['作业状态：连续管下入1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['作业状态：冲洗钻进2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['作业状态：回拖短起2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['作业状态：定点循环2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    else
        title({['作业状态：连续管起出'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    end
    
    frame4 = getframe(fig4);
    im{t} = frame2im(frame4);
end

filename4 = '泵压VS井深VS时间.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename4,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename4,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig5 = figure; % 盘管段出口压力VS井深VS时间
for t = 1:1:Nt
    plot(Time(t)/60,P_ct(t,1)/10^6,'-o');
    hold on;
    xlabel('时间（min）','FontName','黑体','FontSize',12);
    ylabel('盘管段出口压力（MPa）','FontName','黑体','FontSize',12);
    set(gca,'fontsize',12);
    box on; % 显示坐标轴的边框
    grid on; % 显示坐标轴的主网格线
    grid minor; % 显示坐标轴的次网格线
    
    if t==1
        title({['作业状态：连续管下至砂床顶'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['作业状态：冲洗钻进1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['作业状态：回拖短起1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['作业状态：定点循环1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['作业状态：连续管下入1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['作业状态：冲洗钻进2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['作业状态：回拖短起2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['作业状态：定点循环2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    else
        title({['作业状态：连续管起出'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    end
    
    frame5 = getframe(fig5);
    im{t} = frame2im(frame5);
end

filename5 = '盘管段出口压力VS井深VS时间.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename5,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename5,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig6 = figure; % 环空井底压力VS井深VS时间
for t = 1:1:Nt
    plot(Time(t)/60,P_a(t,Nx)/10^6,'-o');
    hold on;
    xlabel('时间（min）','FontName','黑体','FontSize',12);
    ylabel('环空井底压力（MPa）','FontName','黑体','FontSize',12);
    set(gca,'fontsize',12);
    box on; % 显示坐标轴的边框
    grid on; % 显示坐标轴的主网格线
    grid minor; % 显示坐标轴的次网格线
    
    if t==1
        title({['作业状态：连续管下至砂床顶'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['作业状态：冲洗钻进1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['作业状态：回拖短起1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['作业状态：定点循环1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['作业状态：连续管下入1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['作业状态：冲洗钻进2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['作业状态：回拖短起2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['作业状态：定点循环2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    else
        title({['作业状态：连续管起出'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    end
    
    frame6 = getframe(fig6);
    im{t} = frame2im(frame6);
end

filename6 = '环空井底压力VS井深VS时间.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename6,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename6,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig7 = figure; % 环空ECDVS井深VS时间
for t = 1:1:Nt
    plot(P_a(t,13:Nx)./(g*Depth(13:Nx)),Depth(13:Nx),[max(P_a(t,13:Nx)./(g*Depth(13:Nx)))-30,max(P_a(t,13:Nx)./(g*Depth(13:Nx)))+30],[L_coil(t),L_coil(t)],'--','LineWidth',2);
    legend('环空ECD','连续管下深','FontName','黑体','Location','Best');
    xlabel('环空ECD（kg/m3）','FontName','黑体','FontSize',12);
    ylabel('井深（m）','FontName','黑体','FontSize',12);
    set(gca,'fontsize',12);
    set(gca,'YDir','reverse');
    box on; % 显示坐标轴的边框
    grid on; % 显示坐标轴的主网格线
    grid minor; % 显示坐标轴的次网格线
    
    if t==1
        title({['作业状态：连续管下至砂床顶'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['作业状态：冲洗钻进1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['作业状态：回拖短起1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['作业状态：定点循环1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['作业状态：连续管下入1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['作业状态：冲洗钻进2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['作业状态：回拖短起2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['作业状态：定点循环2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    else
        title({['作业状态：连续管起出'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    end
    
    frame7 = getframe(fig7);
    im{t} = frame2im(frame7);
end

filename7 = '环空ECDVS井深VS时间.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename7,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename7,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig8 = figure; % 井筒泡沫质量VS井深VS时间
for t = 1:1:Nt
    plot(gamma_g_a(t,:),Depth,gamma_g_ct(t,:),Depth,[0,1],[L_coil(t),L_coil(t)],'--','LineWidth',2);
    legend('环空泡沫质量','管内泡沫质量','连续管下深','FontName','黑体','Location','Best');
    xlabel('泡沫质量','FontName','黑体','FontSize',12);
    ylabel('井深（m）','FontName','黑体','FontSize',12);
    set(gca,'fontsize',12);
    set(gca,'YDir','reverse');
    box on; % 显示坐标轴的边框
    grid on; % 显示坐标轴的主网格线
    grid minor; % 显示坐标轴的次网格线
    
    if t==1
        title({['作业状态：连续管下至砂床顶'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['作业状态：冲洗钻进1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['作业状态：回拖短起1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['作业状态：定点循环1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['作业状态：连续管下入1'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['作业状态：冲洗钻进2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['作业状态：回拖短起2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['作业状态：定点循环2'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    else
        title({['作业状态：连续管起出'],['时间：',num2str(round(Time(t)/60)),' min    ','连续管下深：',num2str(round(L_coil(t))),' m'],['泵压：',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','黑体','FontSize',12);
        pause(0.1);
    end
    
    frame8 = getframe(fig8);
    im{t} = frame2im(frame8);
end

filename8 = '井筒泡沫质量VS井深VS时间.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 0
        imwrite(A,map,filename8,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename8,'gif','WriteMode','append','DelayTime',0.11);
    end
end