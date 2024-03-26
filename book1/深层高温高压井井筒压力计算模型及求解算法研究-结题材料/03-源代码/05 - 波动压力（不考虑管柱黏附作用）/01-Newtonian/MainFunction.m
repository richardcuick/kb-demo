%% ��ĭ˲̬��ɰ�����ʵ�λ�ƣ�
close all;clear all;clc;
%% �����б���ݣ�����ԭʼ�ռ�ڵ�DEPTH
DATA=csvread("TrajectoryData.csv",1,0); % ������m������б�ǣ��㣩����б��λ�ǣ��㣩�����룩
Dm=DATA(:,1); % ���m
alpha=DATA(:,2); % ��б�ǣ���
phi=DATA(:,3); % ��б��λ�ǣ���
[r,c]=size(DATA); % �ж��ж����С�������
L_s_t=5800; % ɰ��������ȣ����������͹�������̵�����������ȣ���m�����룩
L_s_b=6200; % ɰ���ײ���ȣ����ڳ�ϴ������̵�����������ȣ���m�����룩
Lp=4000; % ���϶�����ȣ�m�����룩

nx_1=100; % ���������϶�����ȶΣ�0��Lp��ԭʼ�ռ������������룩
Nx_1=nx_1+1; % ���������϶�����ȶΣ�0��Lp���ռ�ڵ���
dx_1=Lp/nx_1; % ���������϶�����ȶΣ�0��Lp���ռ䲽����m

for i=1:1:Nx_1
    DEPTH(i)=dx_1*(i-1); % ԭʼ�ռ�ڵ㣬m
end

nx_2=50; % ���϶��������ɰ��������ȶΣ�Lp��L_s_t��ԭʼ�ռ������������룩
Nx_2=nx_2+1; % ���϶��������ɰ��������ȶΣ�Lp��L_s_t���ռ�ڵ���
dx_2=(L_s_t-Lp)/nx_2; % ���϶��������ɰ��������ȶΣ�Lp��L_s_t���ռ䲽����m

for i=Nx_1+1:1:Nx_1+Nx_2-1
    DEPTH(i)=DEPTH(i-1)+dx_2; % ԭʼ�ռ�ڵ㣬m
end

nx_3=50; % ɰ�����������ɰ���ײ���ȶΣ�L_s_t��L_s_b��ԭʼ�ռ������������룩
Nx_3=nx_3+1; % ɰ�����������ɰ���ײ���ȶΣ�L_s_t��L_s_b���ռ�ڵ���
dx_3=(L_s_b-L_s_t)/nx_3; % ɰ�����������ɰ���ײ���ȶΣ�L_s_t��L_s_b���ռ䲽����m

for i=Nx_1+Nx_2:1:Nx_1+Nx_2+Nx_3-2
    DEPTH(i)=DEPTH(i-1)+dx_3; % ԭʼ�ռ�ڵ㣬m
end

Nx=Nx_1+Nx_2+Nx_3-2; % ԭʼ�ռ�ڵ���

%% ����ڵ�
Dsp_wc=[3980]; % ����㾮��侶������m���꾮���ݲ�����

t_pen=2; % �ܳ�ϴ������������룩
L_pen=(L_s_b-L_s_t)/t_pen; % ÿ�γ�ϴ������ȣ�m
for i=1:1:t_pen
    Dsp_pen(i)=L_s_t+L_pen*i; % ����㾮��侶������m����ϴ���������
end

Dsp=[Dsp_wc,Dsp_pen]; % ����㾮��侶������m����������꾮���ݾ�����
n_sp=length(Dsp); % �ж��ж��ٸ������

%% �������㾮��ڵ㣬����ʵ��ʹ�ÿռ�ڵ�Depth
Depth=DEPTH;
for i=1:1:n_sp
    for j=1:1:Nx-1
        if Dsp(i)>DEPTH(j) && Dsp(i)<DEPTH(j+1) % �����ϴ�ʱ��Ҫ��ӿռ�ڵ�
            Depth(j+1)=Dsp(i); % ������㾮��浽Depth��m
            for k=j+1:1:Nx
                Depth(k+1)=DEPTH(k); % ��DEPTH(j+1)֮��Ŀռ�ڵ㸳��Depth
            end
            Nx=Nx+1; % �ռ�ڵ�����һ
            DEPTH=Depth; % ����DEPTH
        end
    end
end
nx=Nx-1; % �ռ�������

Nx_Lp=0; % ���϶���������ϵĿռ�ڵ���
for x=1:1:Nx
    if Depth(x)<=Lp
        Nx_Lp=Nx_Lp+1; % ���϶���������ϵĿռ�ڵ���
    end
end

Nx_L_s_t=0; % ɰ������������ϵĿռ�ڵ���
for x=1:1:Nx
    if Depth(x)<=L_s_t
        Nx_L_s_t=Nx_L_s_t+1; % ɰ������������ϵĿռ�ڵ���
    end
end

for i=1:1:length(Dsp_pen)
    Nx_Dsp_pen(i)=0; % ��ϴ���������ϵĿռ�ڵ���
    for x=1:1:Nx
        if Depth(x)<=Dsp_pen(i)
            Nx_Dsp_pen(i)=Nx_Dsp_pen(i)+1; % ��ϴ���������ϵĿռ�ڵ���
        end
    end
end

%% ��ֵ����ռ�ڵ�Depth��Ӧ�ľ�б��theta
for x=1:1:Nx
    for i=1:1:r-1
        if Depth(x)>=Dm(i) && Depth(x)<Dm(i+1)
            theta(x)=alpha(i)+(alpha(i+1)-alpha(i))*(Depth(x)-Dm(i))/(Dm(i+1)-Dm(i)); % ���Բ�ֵ�õ����ڵ㾮б�ǣ���
        elseif Depth(x)==Dm(i+1)
            theta(x)=alpha(i+1); % ��б�ǣ���
        end
    end
end

data=zeros(length(Depth),2); % ����ʹ�þ����б��
data(:,1)=Depth; % ���m
data(:,2)=theta; % ��б�ǣ���



%% ��ϴ�������1�����ʵ�λ�ƣ�
%% ����ռ䲽��dx����Ӧʱ�䲽��dt
V1=(Dsp_pen(1)-L_s_t)/(1*3600); % ��ϴ����ٶȣ�m/s�����룩
t_1=(Dsp_pen(1)-L_s_t)/V1; % ��ϴ�����ʱ����s

Nt=Nx_Dsp_pen(1)-Nx_L_s_t+1; % ʱ��ڵ���
nt=Nt-1; % ʱ��������

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % �ռ䲽����m
end

Time(1)=0; % ��ϴ�����ʱ����ֵ��s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Dsp_pen(1)-Nt+t)/V1; % ÿ��ϴ���һ���ռ䲽������ʱ�䣬s
    Time(t+1)=Time(t)+dt(t); % ��ϴ�������(t+1)���ռ�ڵ�����������ʱ����s
end

%% ���㲻ͬʱ�������ܳ�ϴ������L_coil���̹ܶ������ܳ���L_reel
L=10000; % �����͹��ܳ���m�����룩
L_wg=8; % ���ڵ�ע��ͷ�����������ܳ��ȣ�m�����룩
L_goose=3; % �������������ܳ��ȣ�m�����룩
D_goose=2; % �������ΰ뾶��m�����룩
H_goose=10; % �����������߶ȣ�m�����룩
L_gr=20; % ����������Ͳ�������ܳ��ȣ�m�����룩
theta_gr=acosd(H_goose/L_gr); % ����������Ͳ����������Ǧ���߼нǣ���
D_r_i=3; % ��Ͳ�ھ���m�����룩
D_r_o=5; % ��Ͳ�⾶��m�����룩
W_r=5; % ��Ͳ��ȣ�m�����룩
D_cable=0.005; % �����⾶��m

L_coil(1)=L_s_t; % ��ʼʱ�������ܵײ���ȣ�m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % ��ʼʱ���̹ܶ������ܳ��ȣ�m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)+dx(Nx_Dsp_pen(1)-Nt+t-1); % ��ϴ�����ȣ�m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % �̹ܶ������ܳ��ȣ�m
end
L_cable=L_coil; % �������볤�ȣ�m

%% �������뼰Ԥ����
D_ct_o_0=0.04445; % �����͹��⾶��m�����룩
L1=2000; % �����ڵ�һ�������ܳ��ȣ�m�����룩
D_ct_i_1=0.03709; % �����ڵ�һ���������ھ���m�����룩
L2=2000; % �����ڵڶ��������ܳ��ȣ�m�����룩
D_ct_i_2=0.03653; % �����ڵڶ����������ھ���m�����룩
L3=2000; % �����ڵ����������ܳ��ȣ�m�����룩
D_ct_i_3=0.03555; % �����ڵ������������ھ���m�����룩
L4=2000; % �����ڵ��Ķ������ܳ��ȣ�m�����룩
D_ct_i_4=0.03489; % �����ڵ��Ķ��������ھ���m�����룩
L5=2000; % �����ڵ���������ܳ��ȣ�m�����룩
D_ct_i_5=0.03409; % �����ڵ�����������ھ���m�����룩

L_t_1=4000; % �ϲ��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_1=0.09718; %0.068;% �ϲ��͹ܣ����׹ܻ����ۣ��ھ���m�����룩
L_t_2=2200; % �²��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_2=0.09718; %0.13970; % �²��͹ܣ����׹ܻ����ۣ��ھ���m�����룩

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % �͹ܣ����׹ܻ����ۣ��ھ���m
        else
            D_t_i(t,x)=D_t_i_2; % �͹ܣ����׹ܻ����ۣ��ھ���m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(1)
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % �����͹��⾶��m
            else
                D_ct_o(t,x)=0; % �����͹��⾶��m
            end
        else
            D_ct_o(t,x)=0; % �����͹��⾶��m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(1)
            if Depth(x)<=L_coil(t) % �����ļ�����������ھ�
                if (L_coil(t)-Depth(x))<=L1 %&& (L_coil(t)-Depth(x))>=0
                    D_ct_i(t,x)=D_ct_i_1; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % �����͹��ھ���m
                end
            else
                D_ct_i(t,x)=0; % �����͹��ھ���m
            end
        else
            D_ct_i(t,x)=0; % �����͹��ھ���m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % �����ļ������������ֵ
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % �����͹��ڽ������m^2
        else
            A_ct(t,x)=0; % �����͹��ڽ������m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % ����ˮ��ֱ����m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % ���ս������m^2
    end
end

%% �������뼰Ԥ����
h_t=25.4*10^(-6); % �͹ܣ����׹ܻ����ۣ����Դֲڶȣ�m�����룩
h_ct=25.4*10^(-6); % �����͹ܾ��Դֲڶȣ�m�����룩
h_a=(h_t+h_ct)/2; % ����ƽ�����Դֲڶȣ�m
epsilon_e=1*10^(-3); % �����������ޣ����룩
epsilon_t=1*10^3; % ���������������룩
g=9.81; % �������ٶȣ�m/s^2��Ĭ�ϣ�

T_0=20; % ��Һ�����¶ȣ��棨���룩
P_0=0.1*10^6; % ��Һ����ѹ����Pa�����룩
rho_l_0=1150; % T_0��P_0�»�Һ�ܶȣ�kg/m^3�����룩
mu_l_0=0.03; % T_0��P_0�»�Һճ�ȣ�Pa*s�����룩
Qv_l_0=0.1/60; % ��Һ���������m^3/s�����룩
Qm_l_0=Qv_l_0*rho_l_0; % ��Һ����������kg/s
rho_g_0=0.655; % ע�����ܶȣ�kg/m^3
mu_g_0=RheologyG(T_0,P_0); % ����ճ�ȣ�Pa*s
Qv_g_0=4/60; % ע�������������m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % ע��������������kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % ��ĭ����������kg/s

D_s=1*10^(-3); % ɰ��ֱ����m�����룩
rho_s=2000; % ɰ���ܶȣ�kg/m^3�����룩
H_s=L_s_b-L_s_t; % �ײ�ɰ���߶ȣ�m
PHI=0.6; % ɰ����ӯ�ȣ����룩
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % ����ɰ����������kg

D_nozzle=4/1000; % ����ֱ����m�����룩
N_nozzle=3; % ������������룩
C=0.95; % ��������ϵ����ȡ0.95�����룩

C0=1.2; % Ư��������ϵ����Ĭ�ϣ�

M_s(1)=0; % ���׽�ɰ����kg/s
for t=2:1:Nt
    M_s(t)=PHI*rho_s*V1*0.25*pi*D_t_i(t,Nx_Dsp_pen(1)-Nt+t)^2; % ���׽�ɰ����kg/s
end

OutPressure=1*10^6; % ����ѹ����Pa�����룩

%% �¶�����
T_i=20; % ��ĭע���¶ȣ���
T_g=0.02; % �����ݶȣ���/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % ����������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������
    end
end
T_a=T_ct; % ������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������

%% ��1��ʱ��ڵ㣨��ʼʱ�̣���ز������㣨���գ�
% ��1���ռ�ڵ㣨���ڣ�����ز�������
P_a(1,1)=OutPressure; % ����ѹ����Pa
rho_g_a(1,1)=DensityG(T_a(1,1),P_a(1,1)); % ���������ܶȣ�kg/m^3
rho_l_a(1,1)=DensityL(rho_l_0,T_0,P_0,T_a(1,1),P_a(1,1)); % ���ջ�Һ�ܶȣ�kg/m^3
mu_g_a(1,1)=RheologyG(T_a(1,1),P_a(1,1)); % ��������ճ�ȣ�Pa*s
mu_l_a(1,1)=RheologyL(mu_l_0,T_0,P_0,T_a(1,1),P_a(1,1)); % ���ջ�Һճ�ȣ�Pa*s
alpha_f_a(1,1)=1; % ��ĭ����
gamma_g_a(1,1)=(Qm_g_0/rho_g_a(1,1))/(Qm_g_0/rho_g_a(1,1)+Qm_l_0/rho_l_a(1,1)); % ��ĭ����
gamma_l_a(1,1)=1-gamma_g_a(1,1); % Һ��������
alpha_g_a(1,1)=alpha_f_a(1,1)*gamma_g_a(1,1); % �������庬��
alpha_l_a(1,1)=alpha_f_a(1,1)*gamma_l_a(1,1); % ���ջ�Һ����
alpha_s(1,1)=0; % �����������
rho_f_a(1,1)=rho_g_a(1,1)*gamma_g_a(1,1)+rho_l_a(1,1)*gamma_l_a(1,1); % ������ĭ�ܶȣ�kg/m^3
mu_f_a(1,1)=mu_g_a(1,1)*gamma_g_a(1,1)+mu_l_a(1,1)*gamma_l_a(1,1); % ������ĭճ�ȣ�Pa*s
V_s(1,1)=0; % �����ٶȣ�m/s
Va_s(1,1)=0; % ���������٣�m/s
mu_s(1,1)=mu_f_a(1,1); % ����ճ�ȣ�Pa*s
Vsr(1,1)=0; % ɰ������ĩ�٣�m/s
Va_f_a(1,1)=Qm_f_0/(rho_f_a(1,1)*A_a(1,1)); % ������ĭ������٣�m/s
Va_g_a(1,1)=Va_f_a(1,1); % �������������٣�m/s
Va_l_a(1,1)=Va_f_a(1,1); % ���ջ�Һ������٣�m/s
V_f_a(1,1)=Va_f_a(1,1)/alpha_f_a(1,1); % ������ĭ���٣�m/s
V_g_a(1,1)=V_f_a(1,1); % �����������٣�m/s
V_l_a(1,1)=V_f_a(1,1); % ���ջ�Һ���٣�m/s

V_m(1,1)=Va_s(1,1)+Va_f_a(1,1); % ���ջ�����ٶȣ�m/s
rho_m(1,1)=alpha_s(1,1)*rho_s+alpha_f_a(1,1)*rho_f_a(1,1); % ���ջ�����ܶȣ�kg/m^3
mu_m(1,1)=alpha_s(1,1)*mu_s(1,1)+alpha_f_a(1,1)*mu_f_a(1,1); % ���ջ����ճ�ȣ�Pa*s
[Ff_a(1,1),f_a(1,1),Re_a(1,1),flow_pattern_a(1,1)]=Friction_annulus(rho_m(1,1),V_m(1,1),mu_m(1,1),D_h(1,1),h_a,rho_f_a(1,1),V_f_a(1,1)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬

% ��2��Nx_L_s_t���ռ�ڵ㴦��ز�������
for x=2:1:Nx_L_s_t
    P_a_ass(1,x)=P_a(1,x-1)+rho_f_a(1,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
    
    err_AnnPressure=1; % ����ѹ��������
    COUNT_AnnPressure=0; % ����ѹ������������ֵ
    while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
        COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
        
        rho_g_a(1,x)=DensityG(T_a(1,x),P_a_ass(1,x)); % ���������ܶȣ�kg/m^3
        rho_l_a(1,x)=DensityL(rho_l_0,T_0,P_0,T_a(1,x),P_a_ass(1,x)); % ���ջ�Һ�ܶȣ�kg/m^3
        mu_g_a(1,x)=RheologyG(T_a(1,x),P_a_ass(1,x)); % ��������ճ�ȣ�Pa*s
        mu_l_a(1,x)=RheologyL(mu_l_0,T_0,P_0,T_a(1,x),P_a_ass(1,x)); % ���ջ�Һճ�ȣ�Pa*s
        alpha_f_a(1,x)=1; % ������ĭ����
        gamma_g_a(1,x)=(Qm_g_0/rho_g_a(1,x))/(Qm_g_0/rho_g_a(1,x)+Qm_l_0/rho_l_a(1,x)); % ��ĭ����
        gamma_l_a(1,x)=1-gamma_g_a(1,x); % Һ��������
        alpha_g_a(1,x)=alpha_f_a(1,x)*gamma_g_a(1,x); % �������庬��
        alpha_l_a(1,x)=alpha_f_a(1,x)*gamma_l_a(1,x); % ���ջ�Һ����
        alpha_s(1,x)=0; % �����������
        rho_f_a(1,x)=rho_g_a(1,x)*gamma_g_a(1,x)+rho_l_a(1,x)*gamma_l_a(1,x); % ������ĭ�ܶȣ�kg/m^3
        mu_f_a(1,x)=mu_g_a(1,x)*gamma_g_a(1,x)+mu_l_a(1,x)*gamma_l_a(1,x); % ������ĭճ�ȣ�Pa*s
        V_s(1,x)=0; % �����ٶȣ�m/s
        Va_s(1,x)=0; % ���������٣�m/s
        mu_s(1,x)=mu_f_a(1,x); % ����ճ�ȣ�Pa*s
        Vsr(1,x)=0; % ɰ������ĩ�٣�m/s
        Va_f_a(1,x)=Qm_f_0/(rho_f_a(1,x)*A_a(1,x)); % ������ĭ������٣�m/s
        Va_g_a(1,x)=Va_f_a(1,x); % �������������٣�m/s
        Va_l_a(1,x)=Va_f_a(1,x); % ���ջ�Һ������٣�m/s
        V_f_a(1,x)=Va_f_a(1,x)/alpha_f_a(1,x); % ������ĭ���٣�m/s
        V_g_a(1,x)=V_f_a(1,x); % �����������٣�m/s
        V_l_a(1,x)=V_f_a(1,x); % ���ջ�Һ���٣�m/s
        
        V_m(1,x)=Va_s(1,x)+Va_f_a(1,x); % ���ջ�����ٶȣ�m/s
        rho_m(1,x)=alpha_s(1,x)*rho_s+alpha_f_a(1,x)*rho_f_a(1,x); % ���ջ�����ܶȣ�kg/m^3
        mu_m(1,x)=alpha_s(1,x)*mu_s(1,x)+alpha_f_a(1,x)*mu_f_a(1,x); % ���ջ����ճ�ȣ�Pa*s
        [Ff_a(1,x),f_a(1,x),Re_a(1,x),flow_pattern_a(1,x)]=Friction_annulus(rho_m(1,x),V_m(1,x),mu_m(1,x),D_h(1,x),h_a,rho_f_a(1,x),V_f_a(1,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
        P_a(1,x)=-rho_m(1,x)*V_m(1,x)^2+P_a(1,x-1)+rho_m(1,x-1)*V_m(1,x-1)^2+((rho_m(1,x)*g*cosd(theta(x))+Ff_a(1,x)+rho_m(1,x-1)*g*cosd(theta(x-1))+Ff_a(1,x-1))*dx(x-1))/2; % ����ѹ����Pa
        
        err_AnnPressure=abs(P_a(1,x)-P_a_ass(1,x))/P_a_ass(1,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
        P_a_ass(1,x)=P_a(1,x); % �µĻ���ѹ������ֵ��Pa
    end
end

% ��Nx_L_s_t+1��Nx���ռ�ڵ㴦��ز�������
for x=Nx_L_s_t+1:1:Nx
    P_a_ass(1,x)=P_a(1,x-1)+rho_f_a(1,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
    
    err_AnnPressure=1; % ����ѹ��������
    COUNT_AnnPressure=0; % ����ѹ������������ֵ
    while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
        COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
        
        rho_g_a(1,x)=DensityG(T_a(1,x),P_a_ass(1,x)); % ���������ܶȣ�kg/m^3
        rho_l_a(1,x)=DensityL(rho_l_0,T_0,P_0,T_a(1,x),P_a_ass(1,x)); % ���ջ�Һ�ܶȣ�kg/m^3
        mu_g_a(1,x)=RheologyG(T_a(1,x),P_a_ass(1,x)); % ��������ճ�ȣ�Pa*s
        mu_l_a(1,x)=RheologyL(mu_l_0,T_0,P_0,T_a(1,x),P_a_ass(1,x)); % ���ջ�Һճ�ȣ�Pa*s
        alpha_f_a(1,x)=1-PHI; % ������ĭ����
        gamma_g_a(1,x)=(Qm_g_0/rho_g_a(1,x))/(Qm_g_0/rho_g_a(1,x)+Qm_l_0/rho_l_a(1,x)); % ��ĭ����
        gamma_l_a(1,x)=1-gamma_g_a(1,x); % Һ��������
        alpha_g_a(1,x)=alpha_f_a(1,x)*gamma_g_a(1,x); % �������庬��
        alpha_l_a(1,x)=alpha_f_a(1,x)*gamma_l_a(1,x); % ���ջ�Һ����
        alpha_s(1,x)=PHI; % �����������
        rho_f_a(1,x)=rho_g_a(1,x)*gamma_g_a(1,x)+rho_l_a(1,x)*gamma_l_a(1,x); % ������ĭ�ܶȣ�kg/m^3
        mu_f_a(1,x)=mu_g_a(1,x)*gamma_g_a(1,x)+mu_l_a(1,x)*gamma_l_a(1,x); % ������ĭճ�ȣ�Pa*s
        V_s(1,x)=0; % �����ٶȣ�m/s
        Va_s(1,x)=0; % ���������٣�m/s
        mu_s(1,x)=mu_f_a(1,x); % ����ճ�ȣ�Pa*s
        Vsr(1,x)=0; % ɰ������ĩ�٣�m/s
        Va_f_a(1,x)=0;%Qm_f_0/(rho_f_a(1,x)*A_a(1,x)); % ������ĭ������٣�m/s
        Va_g_a(1,x)=Va_f_a(1,x); % �������������٣�m/s
        Va_l_a(1,x)=Va_f_a(1,x); % ���ջ�Һ������٣�m/s
        V_f_a(1,x)=Va_f_a(1,x)/alpha_f_a(1,x); % ������ĭ���٣�m/s
        V_g_a(1,x)=V_f_a(1,x); % �����������٣�m/s
        V_l_a(1,x)=V_f_a(1,x); % ���ջ�Һ���٣�m/s
        
        V_m(1,x)=Va_s(1,x)+Va_f_a(1,x); % ���ջ�����ٶȣ�m/s
        rho_m(1,x)=alpha_s(1,x)*rho_s+alpha_f_a(1,x)*rho_f_a(1,x); % ���ջ�����ܶȣ�kg/m^3
        mu_m(1,x)=alpha_s(1,x)*mu_s(1,x)+alpha_f_a(1,x)*mu_f_a(1,x); % ���ջ����ճ�ȣ�Pa*s
        [Ff_a(1,x),f_a(1,x),Re_a(1,x),flow_pattern_a(1,x)]=Friction_annulus(rho_m(1,x),V_m(1,x),mu_m(1,x),D_h(1,x),h_a,rho_f_a(1,x),V_f_a(1,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
        P_a(1,x)=-rho_m(1,x)*V_m(1,x)^2+P_a(1,x-1)+rho_m(1,x-1)*V_m(1,x-1)^2+((rho_m(1,x)*g*cosd(theta(x))+Ff_a(1,x)+rho_m(1,x-1)*g*cosd(theta(x-1))+Ff_a(1,x-1))*dx(x-1))/2; % ����ѹ����Pa
        
        err_AnnPressure=abs(P_a(1,x)-P_a_ass(1,x))/P_a_ass(1,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
        P_a_ass(1,x)=P_a(1,x); % �µĻ���ѹ������ֵ��Pa
    end
end

%% ��2��Nt��ʱ��ڵ���ز������㣨���գ�
for t=2:1:Nt
    P_a(t,Nx_Dsp_pen(1)-Nt+t)=P_a(t-1,Nx_Dsp_pen(1)-Nt+t);  % ���չܵ�ѹ������ֵ��Pa
    
    err_OutPressure=1; % ����ѹ��������
    COUNT_OutPressure=0; % ����ѹ������������ֵ
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % ����ѹ����������
        
        % �ܵף���Nx_Dsp_pen(1)-Nt+t���ռ�ڵ㣩����ز�������
        rho_g_a(t,Nx_Dsp_pen(1)-Nt+t)=DensityG(T_a(t,Nx_Dsp_pen(1)-Nt+t),P_a(t,Nx_Dsp_pen(1)-Nt+t)); % ���������ܶȣ�kg/m^3
        rho_l_a(t,Nx_Dsp_pen(1)-Nt+t)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)-Nt+t),P_a(t,Nx_Dsp_pen(1)-Nt+t)); % ���ջ�Һ�ܶȣ�kg/m^3
        mu_g_a(t,Nx_Dsp_pen(1)-Nt+t)=RheologyG(T_a(t,Nx_Dsp_pen(1)-Nt+t),P_a(t,Nx_Dsp_pen(1)-Nt+t)); % ��������ճ�ȣ�Pa*s
        mu_l_a(t,Nx_Dsp_pen(1)-Nt+t)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)-Nt+t),P_a(t,Nx_Dsp_pen(1)-Nt+t)); % ���ջ�Һճ�ȣ�Pa*s
        gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t)=(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)-Nt+t))/(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)-Nt+t)+Qm_l_0/rho_l_a(t,Nx_Dsp_pen(1)-Nt+t)); % ��ĭ����
        gamma_l_a(t,Nx_Dsp_pen(1)-Nt+t)=1-gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t); % Һ��������
        alpha_g_a(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_a(t-1,Nx_Dsp_pen(1)-Nt+t)*gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t); % �������庬��
        alpha_l_a(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_a(t-1,Nx_Dsp_pen(1)-Nt+t)*gamma_l_a(t,Nx_Dsp_pen(1)-Nt+t); % ���ջ�Һ����
        rho_f_a(t,Nx_Dsp_pen(1)-Nt+t)=rho_g_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t)+rho_l_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_a(t,Nx_Dsp_pen(1)-Nt+t); % ������ĭ�ܶȣ�kg/m^3
        mu_f_a(t,Nx_Dsp_pen(1)-Nt+t)=mu_g_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t)+mu_l_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_a(t,Nx_Dsp_pen(1)-Nt+t); % ������ĭճ�ȣ�Pa*s
        mu_s(t,Nx_Dsp_pen(1)-Nt+t)=mu_f_a(t,Nx_Dsp_pen(1)-Nt+t); % ����ճ�ȣ�Pa*s
        Va_f_a(t,Nx_Dsp_pen(1)-Nt+t)=Qm_f_0/(A_a(t,Nx_Dsp_pen(1)-Nt+t)*rho_f_a(t,Nx_Dsp_pen(1)-Nt+t)); % ������ĭ������٣�m/s
        Va_g_a(t,Nx_Dsp_pen(1)-Nt+t)=Va_f_a(t,Nx_Dsp_pen(1)-Nt+t); % �������������٣�m/s
        Va_l_a(t,Nx_Dsp_pen(1)-Nt+t)=Va_f_a(t,Nx_Dsp_pen(1)-Nt+t); % ���ջ�Һ������٣�m/s
        Va_s(t,Nx_Dsp_pen(1)-Nt+t)=M_s(t)/(A_a(t,Nx_Dsp_pen(1)-Nt+t)*rho_s); % ��м����ٶȣ�m/s
        Vsr(t,Nx_Dsp_pen(1)-Nt+t)=12*(mu_f_a(t,Nx_Dsp_pen(1)-Nt+t)/(rho_f_a(t,Nx_Dsp_pen(1)-Nt+t)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Dsp_pen(1)-Nt+t))/rho_f_a(t,Nx_Dsp_pen(1)-Nt+t))*((rho_f_a(t,Nx_Dsp_pen(1)-Nt+t)*D_s/mu_f_a(t,Nx_Dsp_pen(1)-Nt+t))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
        alpha_s(t,Nx_Dsp_pen(1)-Nt+t)=Va_s(t,Nx_Dsp_pen(1)-Nt+t)/(C0*(Va_s(t,Nx_Dsp_pen(1)-Nt+t)+Va_f_a(t,Nx_Dsp_pen(1)-Nt+t))-Vsr(t,Nx_Dsp_pen(1)-Nt+t));  % �����������
        V_s(t,Nx_Dsp_pen(1)-Nt+t)=Va_s(t,Nx_Dsp_pen(1)-Nt+t)/alpha_s(t,Nx_Dsp_pen(1)-Nt+t); % ��м�ٶȣ�m/s
        
        alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t)=1-alpha_s(t,Nx_Dsp_pen(1)-Nt+t); % ������ĭ����
        alpha_g_a(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_a(t,Nx_Dsp_pen(1)-Nt+t); % �������庬��
        alpha_l_a(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_a(t,Nx_Dsp_pen(1)-Nt+t); % ���ջ�Һ����
        V_f_a(t,Nx_Dsp_pen(1)-Nt+t)=Va_f_a(t,Nx_Dsp_pen(1)-Nt+t)/alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t); % ������ĭ���٣�m/s
        V_g_a(t,Nx_Dsp_pen(1)-Nt+t)=V_f_a(t,Nx_Dsp_pen(1)-Nt+t); % �����������٣�m/s
        V_l_a(t,Nx_Dsp_pen(1)-Nt+t)=V_f_a(t,Nx_Dsp_pen(1)-Nt+t); % ���ջ�Һ���٣�m/s
        
        V_m(t,Nx_Dsp_pen(1)-Nt+t)=Va_s(t,Nx_Dsp_pen(1)-Nt+t)+Va_f_a(t,Nx_Dsp_pen(1)-Nt+t); % ���ջ�����ٶȣ�m/s
        rho_m(t,Nx_Dsp_pen(1)-Nt+t)=alpha_s(t,Nx_Dsp_pen(1)-Nt+t)*rho_s+alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t)*rho_f_a(t,Nx_Dsp_pen(1)-Nt+t); % ���ջ�����ܶȣ�kg/m^3
        mu_m(t,Nx_Dsp_pen(1)-Nt+t)=alpha_s(t,Nx_Dsp_pen(1)-Nt+t)*mu_s(t,Nx_Dsp_pen(1)-Nt+t)+alpha_f_a(t,Nx_Dsp_pen(1)-Nt+t)*mu_f_a(t,Nx_Dsp_pen(1)-Nt+t); % ���ջ����ճ�ȣ�Pa*s
        [Ff_a(t,Nx_Dsp_pen(1)-Nt+t),f_a(t,Nx_Dsp_pen(1)-Nt+t),Re_a(t,Nx_Dsp_pen(1)-Nt+t),flow_pattern_a(t,Nx_Dsp_pen(1)-Nt+t)]=Friction_annulus(rho_m(t,Nx_Dsp_pen(1)-Nt+t),V_m(t,Nx_Dsp_pen(1)-Nt+t),mu_m(t,Nx_Dsp_pen(1)-Nt+t),D_h(t,Nx_Dsp_pen(1)-Nt+t),h_a,rho_f_a(t,Nx_Dsp_pen(1)-Nt+t),V_f_a(t,Nx_Dsp_pen(1)-Nt+t)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
        
        % ��Nx_Dsp_pen(1)-Nt+t-1��1���ռ�ڵ㴦��ز�������
        for x=Nx_Dsp_pen(1)-Nt+t-1:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % ����ѹ������ֵ��Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % ���ջ�Һ����
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            
            err_NodePressure=1; % ����ѹ��������
            COUNT_NodePressure=0; % ����ѹ������������ֵ
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % ����ѹ����������
                
                % ���߷��������������
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % ���������������ֵ1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % ���������������ֵ2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % ������������������
                COUNT_NodeEg=0; % ��������
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % �����������Ϊalpha_s_ass1ʱ������ֵ��
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % �����ٶȣ�m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % ���������٣�m/s
                    alpha_f_ass1=1-alpha_s_ass1; % ��ĭ�������
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % ��ĭ�ٶȣ�m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % ��ĭ������٣�m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % ���ջ�����ٶȣ�m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % ���������������ֵ
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % ����ĺ��������Ľ������ʵ�����������
                    
                    % �����������Ϊalpha_s_ass2ʱ������ֵ��
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % ���������غ㷽����ɢ��ʽ�м�ֵ����
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % �����ٶȣ�m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % ���������٣�m/s
                    alpha_f_ass2=1-alpha_s_ass2; % ��ĭ�������
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % ��ĭ�ٶȣ�m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % ��ĭ������٣�m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % ���ջ�����ٶȣ�m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % ���������������ֵ
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % ����ĺ��������Ľ������ʵ�����������
                    
                    % ���߷���������������
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % �µĹ��������������ֵalpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % ������������������
                    alpha_s_ass1=alpha_s_ass2; % �µĹ��������������ֵ1
                    alpha_s_ass2=alpha_s_ass3; % �µĹ��������������ֵ2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % ���������õ�����ʵ�����������ֵ����alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % �������������С��һ��ֵʱ����Ϊ���������Ϊ0�����ڷ�ֹ���ֺ����ļ������
                    alpha_s(t,x)=0; % �����������
                    V_s(t,x)=0; % �����ٶȣ�m/s
                    Va_s(t,x)=0; % �������ٶȣ�m/s
                    mu_s(t,x)=0; % ����ճ�ȣ�Pa*s
                    Vsr(t,x)=0; % ��м����ĩ�٣�m/s
                    
                    alpha_f_a(t,x)=1; % ��ĭ�������                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                else
                    V_s(t,x)=V_s_ass2; % �����ٶȣ�m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % ���������٣�m/s
                    Vsr(t,x)=Vsr_ass2; % ��м����ĩ�٣�m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % ��ĭ�ܶȣ�kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % ��ĭ�������
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    V_f_a(t,x)=V_f_ass2; % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % ���ջ�����ٶȣ�m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % ����ѹ������ֵ��Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % ����ѹ��������
                P_a(t,x)=P_new; % �µĻ���ѹ������ֵ��Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % ����ѹ��������        
        if (P_a(t,1)-OutPressure)>0 % ���ݳ���ѹ�������������Ի��չܵ�ѹ������ֵ���е���
            P_a(t,Nx_Dsp_pen(1)-Nt+t)=P_a(t,Nx_Dsp_pen(1)-Nt+t)-(P_a(t,1)-OutPressure)/2; % �µĻ��չܵ�ѹ������ֵ��Pa
        else
            P_a(t,Nx_Dsp_pen(1)-Nt+t)=P_a(t,Nx_Dsp_pen(1)-Nt+t)-(P_a(t,1)-OutPressure)/2*0.3; % �µĻ��չܵ�ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(1)-Nt+t+1��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)-Nt+t+1:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1-PHI; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=PHI; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
end

%% ��1��Ntʱ��ڵ���ͷѹ������
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Dsp_pen(1)-Nt+t); % ������ĭ���������m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Dsp_pen(1)-Nt+t)); % �����ɰ����ѹ����Pa����������������V_nozzle��m/s��
end

%% ��1��Ntʱ��ڵ���ز������㣨�������ڣ�
for t=1:1:Nt
    % �ܵף���Nx_Dsp_pen(1)-Nt+t���ռ�ڵ㣩����ز�������
    P_ct(t,Nx_Dsp_pen(1)-Nt+t)=P_a(t,Nx_Dsp_pen(1)-Nt+t)+delta_P_SWT(t); % ����ѹ����Pa
    rho_g_ct(t,Nx_Dsp_pen(1)-Nt+t)=DensityG(T_ct(t,Nx_Dsp_pen(1)-Nt+t),P_ct(t,Nx_Dsp_pen(1)-Nt+t)); % ���������ܶȣ�kg/m^3
    rho_l_ct(t,Nx_Dsp_pen(1)-Nt+t)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)-Nt+t),P_ct(t,Nx_Dsp_pen(1)-Nt+t)); % ���ڻ�Һ�ܶȣ�kg/m^3
    mu_g_ct(t,Nx_Dsp_pen(1)-Nt+t)=RheologyG(T_ct(t,Nx_Dsp_pen(1)-Nt+t),P_ct(t,Nx_Dsp_pen(1)-Nt+t)); % ��������ճ�ȣ�Pa*s
    mu_l_ct(t,Nx_Dsp_pen(1)-Nt+t)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)-Nt+t),P_ct(t,Nx_Dsp_pen(1)-Nt+t)); % ���ڻ�Һճ�ȣ�Pa*s
    alpha_f_ct(t,Nx_Dsp_pen(1)-Nt+t)=1; % ������ĭ����
    gamma_g_ct(t,Nx_Dsp_pen(1)-Nt+t)=(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)-Nt+t))/(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)-Nt+t)+Qm_l_0/rho_l_ct(t,Nx_Dsp_pen(1)-Nt+t)); % ������ĭ����
    gamma_l_ct(t,Nx_Dsp_pen(1)-Nt+t)=1-gamma_g_ct(t,Nx_Dsp_pen(1)-Nt+t); % ����Һ��������
    alpha_g_ct(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_ct(t,Nx_Dsp_pen(1)-Nt+t); % �������庬��
    alpha_l_ct(t,Nx_Dsp_pen(1)-Nt+t)=alpha_f_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_ct(t,Nx_Dsp_pen(1)-Nt+t); % ���ڻ�Һ����
    rho_f_ct(t,Nx_Dsp_pen(1)-Nt+t)=rho_g_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_ct(t,Nx_Dsp_pen(1)-Nt+t)+rho_l_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_ct(t,Nx_Dsp_pen(1)-Nt+t); % ������ĭ�ܶȣ�kg/m^3
    mu_f_ct(t,Nx_Dsp_pen(1)-Nt+t)=mu_g_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_g_ct(t,Nx_Dsp_pen(1)-Nt+t)+mu_l_ct(t,Nx_Dsp_pen(1)-Nt+t)*gamma_l_ct(t,Nx_Dsp_pen(1)-Nt+t); % ������ĭճ�ȣ�Pa*s
    V_f_ct(t,Nx_Dsp_pen(1)-Nt+t)=Qm_f_0/(rho_f_ct(t,Nx_Dsp_pen(1)-Nt+t)*A_ct(t,Nx_Dsp_pen(1)-Nt+t)); % �����������٣�m/s
    [Ff_ct(t,Nx_Dsp_pen(1)-Nt+t),f_ct(t,Nx_Dsp_pen(1)-Nt+t),Re_ct(t,Nx_Dsp_pen(1)-Nt+t),flow_pattern_ct(t,Nx_Dsp_pen(1)-Nt+t)]=Friction_coiledtubing(rho_f_ct(t,Nx_Dsp_pen(1)-Nt+t),V_f_ct(t,Nx_Dsp_pen(1)-Nt+t),mu_f_ct(t,Nx_Dsp_pen(1)-Nt+t),D_ct_i(t,Nx_Dsp_pen(1)-Nt+t),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    % ��Nx_Dsp_pen(1)-Nt+t-1��1���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)-Nt+t-1:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % ����ѹ������ֵ��Pa
        
        err_DriPipePressure=1; % ����ѹ��������
        COUNT_DriPipePressure=0; % ����ѹ������������ֵ
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % ����ѹ����������
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һ�ܶȣ�kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һճ�ȣ�Pa*s
            alpha_f_ct(t,x)=1; % ������ĭ����
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % ������ĭ����
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % �������庬��
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % ���ڻ�Һ����
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭճ�ȣ�Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % ������ĭ���٣�m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % ����ѹ����Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % �������ѹ������ֵ�����ֵ֮���������
            P_ct_ass(t,x)=P_ct(t,x); % �µĹ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(1)-Nt+t+1��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)-Nt+t+1:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % �������庬��
        alpha_l_ct(t,x)=alpha_l_a(t,x); % ���ڻ�Һ����
        gamma_g_ct(t,x)=gamma_g_a(t,x); % ������ĭ����
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
        rho_f_ct(t,x)=rho_m(t,x); % ���������ܶȣ�kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % ��������ճ�ȣ�Pa*s
        V_f_ct(t,t)=V_m(t,x); % �����������٣�sm/s
        Ff_ct(t,x)=Ff_a(t,x); % �������嵥λ����Ħ��ѹ����Pa/m��
        f_ct(t,x)=f_a(t,x); % �������巶��Ħ������
        Re_ct(t,x)=Re_a(t,x); % ����������ŵ����������̬
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % ����������̬
        P_ct(t,x)=P_a(t,x); % ����ѹ����Pa
    end
end

%% ��1��Ntʱ��ڵ��̹ܶγ���ѹ������ѹ���㣨�������ڣ�
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % ��������������٣�m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % ����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % �̹ܶ����嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % �����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % ���ڵ�ע��ͷ������Ħ��ѹ����Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % ��������Ħ��ѹ����Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % ����������Ͳ��Ħ��ѹ����Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % �̹ܶ�Ħ��ѹ����Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % �̹ܶγ���ѹ����Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % ��ѹ��Pa
end

%% ���㾮���ۻ���ɰ��
M_w_tem(1)=0; % ����˲ʱ��ɰ����ֵ��kg
M_w_tot(1)=0; % �����ۻ���ɰ����ֵ��kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮��˲ʱ��ɰ����kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮���ۻ���ɰ����kg
end

%% ����ĩ̬��Ͳ���ɰŨ��
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(1)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% ���㻷����ĭ����ƽ��ֵV_l_mre
V_f_mre=0; % ���շ���ƽ��ֵ��m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % ���շ���ƽ��ֵ��m/s
    end
end

%% ����������ĩ��ƽ��ֵVsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % ��Vsr��������ȡֵ����ĳ���ĩ�٣�m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % �������ĩ��ƽ��ֵ��m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % �������ĩ��ƽ��ֵ��m/s
end

%% �ж��Ƿ���Ч��ɰ�����շ���ƽ��ֵ����2��ɰ������ĩ��ƽ��ֵ�����Ƿ���ɳ�ɰ
fprintf("��ϴ�������1��\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % �����д���˵����ǰ�����¡���ɳ�ɰ��
    else
        fprintf("Sand Cleanout UnFinished!\n"); % �����д���˵����ǰ�����¡�δ��ɳ�ɰ��
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
end

%% ��ͲECD����
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % ����ECD��kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % ����ECD��kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
end

%% ��ϴ�������ĩ״̬����
alpha_f_a_1=alpha_f_a(Nt,:); % ��ϴ�������ĩ״̬��ĭ����
alpha_g_a_1=alpha_g_a(Nt,:); % ��ϴ�������ĩ״̬���ຬ��
alpha_l_a_1=alpha_l_a(Nt,:); % ��ϴ�������ĩ״̬Һ�ຬ��
alpha_s_1=alpha_s(Nt,:); % ��ϴ�������ĩ״̬���ຬ��
f_a_1=f_a(Nt,:); % ��ϴ�������ĩ״̬����Ħ������
Ff_a_1=Ff_a(Nt,:); % ��ϴ�������ĩ״̬���յ�λ����Ħ��ѹ����Pa/m
flow_pattern_a_1=flow_pattern_a(Nt,:); % ��ϴ�������ĩ״̬����������̬
gamma_g_a_1=gamma_g_a(Nt,:); % ��ϴ�������ĩ״̬��ĭ����
gamma_l_a_1=gamma_l_a(Nt,:); % ��ϴ�������ĩ״̬Һ��������
mu_f_a_1=mu_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭճ�ȣ�Pa*s
mu_g_a_1=mu_g_a(Nt,:); % ��ϴ�������ĩ״̬����ճ�ȣ�Pa*s
mu_l_a_1=mu_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ��ճ�ȣ�Pa*s
mu_m_1=mu_m(Nt,:); % ��ϴ�������ĩ״̬���ջ����ճ�ȣ�Pa*s
mu_s_1=mu_s(Nt,:); % ��ϴ�������ĩ״̬����ճ�ȣ�Pa*s
P_a_1=P_a(Nt,:); % ��ϴ�������ĩ״̬����ѹ����Pa
Re_a_1=Re_a(Nt,:); % ��ϴ�������ĩ״̬������ŵ��
rho_f_a_1=rho_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ�ܶȣ�kg/m^3
rho_g_a_1=rho_g_a(Nt,:); % ��ϴ�������ĩ״̬���������ܶȣ�kg/m^3
rho_l_a_1=rho_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ���ܶȣ�kg/m^3
rho_m_1=rho_m(Nt,:); % ��ϴ�������ĩ״̬���ջ�����ܶȣ�kg/m^3
V_f_a_1=V_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ���٣�m/s
V_g_a_1=V_g_a(Nt,:); % ��ϴ�������ĩ״̬�������٣�m/s
V_l_a_1=V_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ�����٣�m/s
V_m_1=V_m(Nt,:); % ��ϴ�������ĩ״̬���ջ�������٣�m/s
V_s_1=V_s(Nt,:); % ��ϴ�������ĩ״̬�������٣�m/s
Va_f_a_1=Va_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ������٣�m/s
Va_g_a_1=Va_g_a(Nt,:); % ��ϴ�������ĩ״̬���������٣�m/s
Va_l_a_1=Va_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ�������٣�m/s
Va_s_1=Va_s(Nt,:); % ��ϴ�������ĩ״̬���������٣�m/s
Vsr_1=Vsr(Nt,:); % ��ϴ�������ĩ״̬���໬���ٶȣ�m/s

%% ���ݴ洢����ϴ������̣�
ANS_Nt_1=Nt; % ʱ��ڵ���
ANS_Time_1=Time(1:Nt); % ��ϴʱ�䣬s
ANS_alpha_g_a_1=alpha_g_a(1:Nt,:); % ���������������
ANS_alpha_l_a_1=alpha_l_a(1:Nt,:); % ����Һ���������
ANS_alpha_s_1=alpha_s(1:Nt,:); % ���ຬ��
ANS_P_ct_1=P_ct(1:Nt,:); % ����ѹ����Pa��
ANS_delta_P_SWT_1=delta_P_SWT(1:Nt); % ��ɰ����ѹ����Pa��
ANS_P_a_1=P_a(1:Nt,:); % ����ѹ����Pa��
ANS_T_a_1=T_a(1:Nt,:); % �����¶ȣ��棩
ANS_T_ct_1=T_ct(1:Nt,:); % �����¶ȣ��棩
ANS_P_coil_1=P_coil(1:Nt); % �̹ܶγ���ѹ����Pa��
ANS_P_pump_1=P_pump(1:Nt); % ��ѹ��Pa��
ANS_M_w_tem_1=M_w_tem(1:Nt); % ����˲ʱ��ɰ����kg��
ANS_M_w_tot_1=M_w_tot(1:Nt); % �����ۻ���ɰ����kg��
ANS_Va_s_1=Va_s(1:Nt,:); % ��м����ٶȣ�m/s��
ANS_Va_f_a_1=Va_f_a(1:Nt,:); % ������ĭ������٣�m/s��
ANS_V_s_1=V_s(1:Nt,:); % ��м�����ٶȣ�m/s��
ANS_V_f_a_1=V_f_a(1:Nt,:); % ������ĭ���٣�m/s��
ANS_alpha_g_ct_1=alpha_g_ct(1:Nt,:); % ���������������
ANS_alpha_l_ct_1=alpha_l_ct(1:Nt,:); % ����Һ���������
ANS_gamma_g_ct_1=gamma_g_ct(1:Nt,:); % ������ĭ����
ANS_gamma_g_a_1=gamma_g_a(1:Nt,:); % ������ĭ����
ANS_ECD_a_1=ECD_a(1:Nt,:); % ����ECD��kg/m^3��
ANS_L_coil_1=L_coil(1:Nt); % ���������m��



%% ���϶������1�����ʵ�λ�ƣ�
%% ������ֵ����
dt=zeros(); % ʱ�䲽����s
dx=zeros(); % �ռ䲽����m
Time=zeros(); % ������������ʱ����s
L_coil=zeros(); % ���������m
L_reel=zeros(); % �̹ܶγ��ȣ�m
D_t_i=zeros(); % �͹��ھ���m
D_ct_o=zeros(); % �����͹��⾶��m
D_ct_i=zeros(); % �����͹��ھ���m
M_s=zeros(); % ���׽�ɰ����kg/s
P_a=zeros(); % ����ѹ����Pa
rho_g_a=zeros(); % �����ܶȣ�kg/m^3
rho_l_a=zeros(); % ����Һ���ܶȣ�kg/m^3
rho_f_a=zeros(); % ������ĭ�ܶȣ�kg/m^3
mu_g_a=zeros(); % ����ճ�ȣ�Pa*s
mu_l_a=zeros(); % ����Һ��ճ�ȣ�Pa*s
mu_f_a=zeros(); % ������ĭճ�ȣ�Pa*s
mu_s=zeros(); % ����ճ�ȣ�Pa*s
alpha_g_a=zeros(); % ���ຬ��
alpha_l_a=zeros(); % Һ�ຬ��
alpha_f_a=zeros(); % ��ĭ����
alpha_s=zeros(); % ���ຬ��
Va_g_a=zeros(); % ���������٣�m/s
Va_l_a=zeros(); % ����Һ�������٣�m/s
Va_f_a=zeros(); % ������ĭ������٣�m/s
Va_s=zeros(); % ���������٣�m/s
V_g_a=zeros(); % �������٣�m/s
V_l_a=zeros(); % ����Һ�����٣�m/s
V_f_a=zeros(); % ������ĭ���٣�m/s
V_s=zeros(); % �������٣�m/s
Vsr=zeros(); % ���໬���ٶȣ�m/s
gamma_g_a=zeros(); % ��ĭ����
gamma_l_a=zeros(); % Һ��������
V_m=zeros(); % ���ջ�������٣�m/s
rho_m=zeros(); % ���ջ�����ܶȣ�kg/m^3
mu_m=zeros(); % ���ջ����ճ�ȣ�Pa*s
Ff_a=zeros(); % ���յ�λ����Ħ��ѹ����Pa/m
f_a=zeros(); % ����Ħ������
Re_a=zeros(); % ������ŵ��
flow_pattern_a=zeros(); % ����������̬

%% ����ռ䲽��dx����Ӧʱ�䲽��dt
V2=(Dsp_pen(1)-Lp)/(3*3600); % ���϶����ٶȣ�m/s�����룩
t_2=(Dsp_pen(1)-Lp)/V2; % ���϶�����ʱ����s

Nt=Nx_Dsp_pen(1)-Nx_Lp+1; % ʱ��ڵ���
nt=Nt-1; % ʱ��������

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % �ռ䲽����m
end

Time(1)=0; % �����������ʱ����ֵ��s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Dsp_pen(1)-t)/V2; % ÿ���һ���ռ䲽������ʱ�䣬s
    Time(t+1)=Time(t)+dt(t); % �������������(t+1)���ռ�ڵ�����������ʱ����s
end

%% ���㲻ͬʱ��������������L_coil���̹ܶ������ܳ���L_reel
L=10000; % �����͹��ܳ���m�����룩
L_wg=8; % ���ڵ�ע��ͷ�����������ܳ��ȣ�m�����룩
L_goose=3; % �������������ܳ��ȣ�m�����룩
D_goose=2; % �������ΰ뾶��m�����룩
H_goose=10; % �����������߶ȣ�m�����룩
L_gr=20; % ����������Ͳ�������ܳ��ȣ�m�����룩
theta_gr=acosd(H_goose/L_gr); % ����������Ͳ����������Ǧ���߼нǣ���
D_r_i=3; % ��Ͳ�ھ���m�����룩
D_r_o=5; % ��Ͳ�⾶��m�����룩
W_r=5; % ��Ͳ��ȣ�m�����룩
D_cable=0.005; % �����⾶��m

L_coil(1)=Dsp_pen(1); % ��ʼʱ�������ܵײ���ȣ�m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % ��ʼʱ���̹ܶ������ܳ��ȣ�m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)-dx(Nx_Dsp_pen(1)-t+1); % �����������ȣ�m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % �̹ܶ������ܳ��ȣ�m
end
L_cable=L_coil; % �������볤�ȣ�m

%% �������뼰Ԥ����
D_ct_o_0=0.04445; % �����͹��⾶��m�����룩
L1=2000; % �����ڵ�һ�������ܳ��ȣ�m�����룩
D_ct_i_1=0.03709; % �����ڵ�һ���������ھ���m�����룩
L2=2000; % �����ڵڶ��������ܳ��ȣ�m�����룩
D_ct_i_2=0.03653; % �����ڵڶ����������ھ���m�����룩
L3=2000; % �����ڵ����������ܳ��ȣ�m�����룩
D_ct_i_3=0.03555; % �����ڵ������������ھ���m�����룩
L4=2000; % �����ڵ��Ķ������ܳ��ȣ�m�����룩
D_ct_i_4=0.03489; % �����ڵ��Ķ��������ھ���m�����룩
L5=2000; % �����ڵ���������ܳ��ȣ�m�����룩
D_ct_i_5=0.03409; % �����ڵ�����������ھ���m�����룩

L_t_1=4000; % �ϲ��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_1=0.09718; %0.068;%0.09718; % �ϲ��͹ܣ����׹ܻ����ۣ��ھ���m�����룩
L_t_2=2200; % �²��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_2=0.09718; %0.13970; % �²��͹ܣ����׹ܻ����ۣ��ھ���m�����룩

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % �͹ܣ����׹ܻ����ۣ��ھ���m
        else
            D_t_i(t,x)=D_t_i_2; % �͹ܣ����׹ܻ����ۣ��ھ���m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(1)
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % �����͹��⾶��m
            else
                D_ct_o(t,x)=0; % �����͹��⾶��m
            end
        else
            D_ct_o(t,x)=0; % �����͹��⾶��m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(1)
            if Depth(x)<=L_coil(t) % �����ļ�����������ھ�
                if (L_coil(t)-Depth(x))<=L1 && (L_coil(t)-Depth(x))>=0
                    D_ct_i(t,x)=D_ct_i_1; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % �����͹��ھ���m
                end
            else
                D_ct_i(t,x)=0; % �����͹��ھ���m
            end
        else
            D_ct_i(t,x)=0; % �����͹��ھ���m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % �����ļ������������ֵ
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % �����͹��ڽ������m^2
        else
            A_ct(t,x)=0; % �����͹��ڽ������m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % ����ˮ��ֱ����m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % ���ս������m^2
    end
end

%% �������뼰Ԥ����
h_t=25.4*10^(-6); % �͹ܣ����׹ܻ����ۣ����Դֲڶȣ�m�����룩
h_ct=25.4*10^(-6); % �����͹ܾ��Դֲڶȣ�m�����룩
h_a=(h_t+h_ct)/2; % ����ƽ�����Դֲڶȣ�m
epsilon_e=1*10^(-3); % �����������ޣ����룩
epsilon_t=1*10^3; % ���������������룩
g=9.81; % �������ٶȣ�m/s^2��Ĭ�ϣ�

T_0=20; % ��Һ�����¶ȣ��棨���룩
P_0=0.1*10^6; % ��Һ����ѹ����Pa�����룩
rho_l_0=1150; % T_0��P_0�»�Һ�ܶȣ�kg/m^3�����룩
mu_l_0=0.03; % T_0��P_0�»�Һճ�ȣ�Pa*s�����룩
Qv_l_0=0.1/60; % ��Һ���������m^3/s�����룩
Qm_l_0=Qv_l_0*rho_l_0; % ��Һ����������kg/s
rho_g_0=0.655; % ע�����ܶȣ�kg/m^3
mu_g_0=RheologyG(T_0,P_0); % ����ճ�ȣ�Pa*s
Qv_g_0=4/60; % ע�������������m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % ע��������������kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % ��ĭ����������kg/s

D_s=1*10^(-3); % ɰ��ֱ����m�����룩
rho_s=2000; % ɰ���ܶȣ�kg/m^3�����룩
H_s=L_s_b-Nx_Dsp_pen(1); % �ײ�ɰ���߶ȣ�m
PHI=0.6; % ɰ����ӯ�ȣ����룩
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % ����ɰ����������kg

D_nozzle=4/1000; % ����ֱ����m�����룩
N_nozzle=3; % ������������룩
C=0.95; % ��������ϵ����ȡ0.95�����룩

C0=1.2; % Ư��������ϵ����Ĭ�ϣ�

M_s(1)=0; % ���׽�ɰ����kg/s
for t=2:1:Nt
    M_s(t)=0; % ���׽�ɰ����kg/s
end

OutPressure=1*10^6; % ����ѹ����Pa�����룩

%% �¶�����
T_i=20; % ��ĭע���¶ȣ���
T_g=0.02; % �����ݶȣ���/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % ����������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������
    end
end
T_a=T_ct; % ������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������

%% ��1��ʱ��ڵ㣨��ʼʱ�̣���ز�����ֵ���㣨���գ�
for x=1:1:Nx
    P_a(1,x)=P_a_1(x); % ����ѹ����Pa
    rho_g_a(1,x)=rho_g_a_1(x); % �����ܶȣ�kg/m^3
    rho_l_a(1,x)=rho_l_a_1(x); % ����Һ���ܶȣ�kg/m^3
    rho_f_a(1,x)=rho_f_a_1(x); % ������ĭ�ܶȣ�kg/m^3
    mu_g_a(1,x)=mu_g_a_1(x); % ����ճ�ȣ�Pa*s
    mu_l_a(1,x)=mu_l_a_1(x); % ����Һ��ճ�ȣ�Pa*s
    mu_f_a(1,x)=mu_f_a_1(x); % ������ĭճ�ȣ�Pa*s
    mu_s(1,x)=mu_s_1(x); % ����ճ�ȣ�Pa*s
    alpha_g_a(1,x)=alpha_g_a_1(x); % ���ຬ��
    alpha_l_a(1,x)=alpha_l_a_1(x); % Һ�ຬ��
    alpha_f_a(1,x)=alpha_f_a_1(x); % ��ĭ����
    alpha_s(1,x)=alpha_s_1(1,x); % ���ຬ��
    Va_g_a(1,x)=Va_g_a_1(x); % ���������٣�m/s
    Va_l_a(1,x)=Va_l_a_1(x); % ����Һ�������٣�m/s
    Va_f_a(1,x)=Va_f_a_1(x); % ������ĭ������٣�m/s
    Va_s(1,x)=Va_s_1(x); % ���������٣�m/s
    V_g_a(1,x)=V_g_a_1(x); % �������٣�m/s
    V_l_a(1,x)=V_l_a_1(x); % ����Һ�����٣�m/s
    V_f_a(1,x)=V_f_a_1(x); % ������ĭ���٣�m/s
    V_s(1,x)=V_s_1(x); % �������٣�m/s
    Vsr(1,x)=Vsr_1(x); % ���໬���ٶȣ�m/s
    gamma_g_a(1,x)=gamma_g_a_1(x); % ��ĭ����
    gamma_l_a(1,x)=gamma_l_a_1(1,x); % Һ��������
    V_m(1,x)=V_m_1(x); % ���ջ�������٣�m/s
    rho_m(1,x)=rho_m_1(x); % ���ջ�����ܶȣ�kg/m^3
    mu_m(1,x)=mu_m_1(x); % ���ջ����ճ�ȣ�Pa*s
    Ff_a(1,x)=Ff_a_1(x); % ���յ�λ����Ħ��ѹ����Pa/m
    f_a(1,x)=f_a_1(x); % ����Ħ������
    Re_a(1,x)=Re_a_1(x); % ������ŵ��
    flow_pattern_a(1,x)=flow_pattern_a_1(x); % ����������̬
end

%% ��2��Nt��ʱ��ڵ���ز������㣨���գ�
for t=2:1:Nt
    P_a(t,Nx_Dsp_pen(1)-t+1)=P_a(t-1,Nx_Dsp_pen(1)-t+1);  % ���չܵ�ѹ������ֵ��Pa
    
    err_OutPressure=1; % ����ѹ��������
    COUNT_OutPressure=0; % ����ѹ������������ֵ
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % ����ѹ����������
        
        % �ܵף���Nx_Dsp_pen(1)-t+1���ռ�ڵ㣩����ز�������
        rho_g_a(t,Nx_Dsp_pen(1)-t+1)=DensityG(T_a(t,Nx_Dsp_pen(1)-t+1),P_a(t,Nx_Dsp_pen(1)-t+1)); % ���������ܶȣ�kg/m^3
        rho_l_a(t,Nx_Dsp_pen(1)-t+1)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)-t+1),P_a(t,Nx_Dsp_pen(1)-t+1)); % ���ջ�Һ�ܶȣ�kg/m^3
        mu_g_a(t,Nx_Dsp_pen(1)-t+1)=RheologyG(T_a(t,Nx_Dsp_pen(1)-t+1),P_a(t,Nx_Dsp_pen(1)-t+1)); % ��������ճ�ȣ�Pa*s
        mu_l_a(t,Nx_Dsp_pen(1)-t+1)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)-t+1),P_a(t,Nx_Dsp_pen(1)-t+1)); % ���ջ�Һճ�ȣ�Pa*s
        gamma_g_a(t,Nx_Dsp_pen(1)-t+1)=(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)-t+1))/(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)-t+1)+Qm_l_0/rho_l_a(t,Nx_Dsp_pen(1)-t+1)); % ��ĭ����
        gamma_l_a(t,Nx_Dsp_pen(1)-t+1)=1-gamma_g_a(t,Nx_Dsp_pen(1)-t+1); % Һ��������
        alpha_g_a(t,Nx_Dsp_pen(1)-t+1)=alpha_f_a(t-1,Nx_Dsp_pen(1)-t+1)*gamma_g_a(t,Nx_Dsp_pen(1)-t+1); % �������庬��
        alpha_l_a(t,Nx_Dsp_pen(1)-t+1)=alpha_f_a(t-1,Nx_Dsp_pen(1)-t+1)*gamma_l_a(t,Nx_Dsp_pen(1)-t+1); % ���ջ�Һ����
        rho_f_a(t,Nx_Dsp_pen(1)-t+1)=rho_g_a(t,Nx_Dsp_pen(1)-t+1)*gamma_g_a(t,Nx_Dsp_pen(1)-t+1)+rho_l_a(t,Nx_Dsp_pen(1)-t+1)*gamma_l_a(t,Nx_Dsp_pen(1)-t+1); % ������ĭ�ܶȣ�kg/m^3
        mu_f_a(t,Nx_Dsp_pen(1)-t+1)=mu_g_a(t,Nx_Dsp_pen(1)-t+1)*gamma_g_a(t,Nx_Dsp_pen(1)-t+1)+mu_l_a(t,Nx_Dsp_pen(1)-t+1)*gamma_l_a(t,Nx_Dsp_pen(1)-t+1); % ������ĭճ�ȣ�Pa*s
        mu_s(t,Nx_Dsp_pen(1)-t+1)=mu_f_a(t,Nx_Dsp_pen(1)-t+1); % ����ճ�ȣ�Pa*s
        Va_f_a(t,Nx_Dsp_pen(1)-t+1)=Qm_f_0/(A_a(t,Nx_Dsp_pen(1)-t+1)*rho_f_a(t,Nx_Dsp_pen(1)-t+1)); % ������ĭ������٣�m/s
        Va_g_a(t,Nx_Dsp_pen(1)-t+1)=Va_f_a(t,Nx_Dsp_pen(1)-t+1); % �������������٣�m/s
        Va_l_a(t,Nx_Dsp_pen(1)-t+1)=Va_f_a(t,Nx_Dsp_pen(1)-t+1); % ���ջ�Һ������٣�m/s
        Va_s(t,Nx_Dsp_pen(1)-t+1)=M_s(t)/(A_a(t,Nx_Dsp_pen(1)-t+1)*rho_s); % ��м����ٶȣ�m/s
        Vsr(t,Nx_Dsp_pen(1)-t+1)=12*(mu_f_a(t,Nx_Dsp_pen(1)-t+1)/(rho_f_a(t,Nx_Dsp_pen(1)-t+1)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Dsp_pen(1)-t+1))/rho_f_a(t,Nx_Dsp_pen(1)-t+1))*((rho_f_a(t,Nx_Dsp_pen(1)-t+1)*D_s/mu_f_a(t,Nx_Dsp_pen(1)-t+1))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
        alpha_s(t,Nx_Dsp_pen(1)-t+1)=Va_s(t,Nx_Dsp_pen(1)-t+1)/(C0*(Va_s(t,Nx_Dsp_pen(1)-t+1)+Va_f_a(t,Nx_Dsp_pen(1)-t+1))-Vsr(t,Nx_Dsp_pen(1)-t+1));  % �����������
        V_s(t,Nx_Dsp_pen(1)-t+1)=0;%Va_s(t,Nx_Dsp_pen(1)-t+1)/alpha_s(t,Nx_Dsp_pen(1)-t+1); % ��м�ٶȣ�m/s
        
        alpha_f_a(t,Nx_Dsp_pen(1)-t+1)=1-alpha_s(t,Nx_Dsp_pen(1)-t+1); % ������ĭ����
        alpha_g_a(t,Nx_Dsp_pen(1)-t+1)=alpha_f_a(t,Nx_Dsp_pen(1)-t+1)*gamma_g_a(t,Nx_Dsp_pen(1)-t+1); % �������庬��
        alpha_l_a(t,Nx_Dsp_pen(1)-t+1)=alpha_f_a(t,Nx_Dsp_pen(1)-t+1)*gamma_l_a(t,Nx_Dsp_pen(1)-t+1); % ���ջ�Һ����
        V_f_a(t,Nx_Dsp_pen(1)-t+1)=Va_f_a(t,Nx_Dsp_pen(1)-t+1)/alpha_f_a(t,Nx_Dsp_pen(1)-t+1); % ������ĭ���٣�m/s
        V_g_a(t,Nx_Dsp_pen(1)-t+1)=V_f_a(t,Nx_Dsp_pen(1)-t+1); % �����������٣�m/s
        V_l_a(t,Nx_Dsp_pen(1)-t+1)=V_f_a(t,Nx_Dsp_pen(1)-t+1); % ���ջ�Һ���٣�m/s
        
        V_m(t,Nx_Dsp_pen(1)-t+1)=Va_s(t,Nx_Dsp_pen(1)-t+1)+Va_f_a(t,Nx_Dsp_pen(1)-t+1); % ���ջ�����ٶȣ�m/s
        rho_m(t,Nx_Dsp_pen(1)-t+1)=alpha_s(t,Nx_Dsp_pen(1)-t+1)*rho_s+alpha_f_a(t,Nx_Dsp_pen(1)-t+1)*rho_f_a(t,Nx_Dsp_pen(1)-t+1); % ���ջ�����ܶȣ�kg/m^3
        mu_m(t,Nx_Dsp_pen(1)-t+1)=alpha_s(t,Nx_Dsp_pen(1)-t+1)*mu_s(t,Nx_Dsp_pen(1)-t+1)+alpha_f_a(t,Nx_Dsp_pen(1)-t+1)*mu_f_a(t,Nx_Dsp_pen(1)-t+1); % ���ջ����ճ�ȣ�Pa*s
        [Ff_a(t,Nx_Dsp_pen(1)-t+1),f_a(t,Nx_Dsp_pen(1)-t+1),Re_a(t,Nx_Dsp_pen(1)-t+1),flow_pattern_a(t,Nx_Dsp_pen(1)-t+1)]=Friction_annulus(rho_m(t,Nx_Dsp_pen(1)-t+1),V_m(t,Nx_Dsp_pen(1)-t+1),mu_m(t,Nx_Dsp_pen(1)-t+1),D_h(t,Nx_Dsp_pen(1)-t+1),h_a,rho_f_a(t,Nx_Dsp_pen(1)-t+1),V_f_a(t,Nx_Dsp_pen(1)-t+1)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
        
        % ��Nx_Dsp_pen(1)-t��1���ռ�ڵ㴦��ز�������
        for x=Nx_Dsp_pen(1)-t:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % ����ѹ������ֵ��Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % ���ջ�Һ����
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            
            err_NodePressure=1; % ����ѹ��������
            COUNT_NodePressure=0; % ����ѹ������������ֵ
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % ����ѹ����������
                
                % ���߷��������������
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % ���������������ֵ1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % ���������������ֵ2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % ������������������
                COUNT_NodeEg=0; % ��������
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % �����������Ϊalpha_s_ass1ʱ������ֵ��
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % �����ٶȣ�m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % ���������٣�m/s
                    alpha_f_ass1=1-alpha_s_ass1; % ��ĭ�������
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % ��ĭ�ٶȣ�m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % ��ĭ������٣�m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % ���ջ�����ٶȣ�m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % ���������������ֵ
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % ����ĺ��������Ľ������ʵ�����������
                    
                    % �����������Ϊalpha_s_ass2ʱ������ֵ��
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % ���������غ㷽����ɢ��ʽ�м�ֵ����
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % �����ٶȣ�m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % ���������٣�m/s
                    alpha_f_ass2=1-alpha_s_ass2; % ��ĭ�������
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % ��ĭ�ٶȣ�m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % ��ĭ������٣�m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % ���ջ�����ٶȣ�m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % ���������������ֵ
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % ����ĺ��������Ľ������ʵ�����������
                    
                    % ���߷���������������
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % �µĹ��������������ֵalpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % ������������������
                    alpha_s_ass1=alpha_s_ass2; % �µĹ��������������ֵ1
                    alpha_s_ass2=alpha_s_ass3; % �µĹ��������������ֵ2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % ���������õ�����ʵ�����������ֵ����alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % �������������С��һ��ֵʱ����Ϊ���������Ϊ0�����ڷ�ֹ���ֺ����ļ������
                    alpha_s(t,x)=0; % �����������
                    V_s(t,x)=0; % �����ٶȣ�m/s
                    Va_s(t,x)=0; % �������ٶȣ�m/s
                    mu_s(t,x)=0; % ����ճ�ȣ�Pa*s
                    Vsr(t,x)=0; % ��м����ĩ�٣�m/s
                    
                    alpha_f_a(t,x)=1; % ��ĭ�������                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                else
                    V_s(t,x)=V_s_ass2; % �����ٶȣ�m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % ���������٣�m/s
                    Vsr(t,x)=Vsr_ass2; % ��м����ĩ�٣�m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % ��ĭ�ܶȣ�kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % ��ĭ�������
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    V_f_a(t,x)=V_f_ass2; % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % ���ջ�����ٶȣ�m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % ����ѹ������ֵ��Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % ����ѹ��������
                P_a(t,x)=P_new; % �µĻ���ѹ������ֵ��Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % ����ѹ��������
        if (P_a(t,1)-OutPressure)>0 % ���ݳ���ѹ�������������Ի��չܵ�ѹ������ֵ���е���
            P_a(t,Nx_Dsp_pen(1)-t+1)=P_a(t,Nx_Dsp_pen(1)-t+1)-(P_a(t,1)-OutPressure)/2; % �µĻ��չܵ�ѹ������ֵ��Pa
        else
            P_a(t,Nx_Dsp_pen(1)-t+1)=P_a(t,Nx_Dsp_pen(1)-t+1)-(P_a(t,1)-OutPressure)/2*0.3; % �µĻ��չܵ�ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(1)-t+2��Nx_Dsp_pen(1)���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)-t+2:1:Nx_Dsp_pen(1)
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=0; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(1)+1��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)+1:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1-PHI; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=PHI; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
end

%% ��1��Ntʱ��ڵ���ͷѹ������
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Dsp_pen(1)-t+1); % ������ĭ���������m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Dsp_pen(1)-t+1)); % �����ɰ����ѹ����Pa����������������V_nozzle��m/s��
end

%% ��1��Ntʱ��ڵ���ز������㣨�������ڣ�
for t=1:1:Nt
    % �ܵף���Nx_Dsp_pen(1)-t+1���ռ�ڵ㣩����ز�������
    P_ct(t,Nx_Dsp_pen(1)-t+1)=P_a(t,Nx_Dsp_pen(1)-t+1)+delta_P_SWT(t); % ����ѹ����Pa
    rho_g_ct(t,Nx_Dsp_pen(1)-t+1)=DensityG(T_ct(t,Nx_Dsp_pen(1)-t+1),P_ct(t,Nx_Dsp_pen(1)-t+1)); % ���������ܶȣ�kg/m^3
    rho_l_ct(t,Nx_Dsp_pen(1)-t+1)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)-t+1),P_ct(t,Nx_Dsp_pen(1)-t+1)); % ���ڻ�Һ�ܶȣ�kg/m^3
    mu_g_ct(t,Nx_Dsp_pen(1)-t+1)=RheologyG(T_ct(t,Nx_Dsp_pen(1)-t+1),P_ct(t,Nx_Dsp_pen(1)-t+1)); % ��������ճ�ȣ�Pa*s
    mu_l_ct(t,Nx_Dsp_pen(1)-t+1)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)-t+1),P_ct(t,Nx_Dsp_pen(1)-t+1)); % ���ڻ�Һճ�ȣ�Pa*s
    alpha_f_ct(t,Nx_Dsp_pen(1)-t+1)=1; % ������ĭ����
    gamma_g_ct(t,Nx_Dsp_pen(1)-t+1)=(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)-t+1))/(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)-t+1)+Qm_l_0/rho_l_ct(t,Nx_Dsp_pen(1)-t+1)); % ������ĭ����
    gamma_l_ct(t,Nx_Dsp_pen(1)-t+1)=1-gamma_g_ct(t,Nx_Dsp_pen(1)-t+1); % ����Һ��������
    alpha_g_ct(t,Nx_Dsp_pen(1)-t+1)=alpha_f_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(1)-t+1); % �������庬��
    alpha_l_ct(t,Nx_Dsp_pen(1)-t+1)=alpha_f_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(1)-t+1); % ���ڻ�Һ����
    rho_f_ct(t,Nx_Dsp_pen(1)-t+1)=rho_g_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(1)-t+1)+rho_l_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(1)-t+1); % ������ĭ�ܶȣ�kg/m^3
    mu_f_ct(t,Nx_Dsp_pen(1)-t+1)=mu_g_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(1)-t+1)+mu_l_ct(t,Nx_Dsp_pen(1)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(1)-t+1); % ������ĭճ�ȣ�Pa*s
    V_f_ct(t,Nx_Dsp_pen(1)-t+1)=Qm_f_0/(rho_f_ct(t,Nx_Dsp_pen(1)-t+1)*A_ct(t,Nx_Dsp_pen(1)-t+1)); % �����������٣�m/s
    [Ff_ct(t,Nx_Dsp_pen(1)-t+1),f_ct(t,Nx_Dsp_pen(1)-t+1),Re_ct(t,Nx_Dsp_pen(1)-t+1),flow_pattern_ct(t,Nx_Dsp_pen(1)-t+1)]=Friction_coiledtubing(rho_f_ct(t,Nx_Dsp_pen(1)-t+1),V_f_ct(t,Nx_Dsp_pen(1)-t+1),mu_f_ct(t,Nx_Dsp_pen(1)-t+1),D_ct_i(t,Nx_Dsp_pen(1)-t+1),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    % ��Nx_Dsp_pen(1)-t��1���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)-t:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % ����ѹ������ֵ��Pa
        
        err_DriPipePressure=1; % ����ѹ��������
        COUNT_DriPipePressure=0; % ����ѹ������������ֵ
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % ����ѹ����������
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һ�ܶȣ�kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һճ�ȣ�Pa*s
            alpha_f_ct(t,x)=1; % ������ĭ����
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % ������ĭ����
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % �������庬��
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % ���ڻ�Һ����
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭճ�ȣ�Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % ������ĭ���٣�m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % ����ѹ����Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % �������ѹ������ֵ�����ֵ֮���������
            P_ct_ass(t,x)=P_ct(t,x); % �µĹ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(1)-t+2��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)-t+2:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % �������庬��
        alpha_l_ct(t,x)=alpha_l_a(t,x); % ���ڻ�Һ����
        gamma_g_ct(t,x)=gamma_g_a(t,x); % ������ĭ����
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
        rho_f_ct(t,x)=rho_m(t,x); % ���������ܶȣ�kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % ��������ճ�ȣ�Pa*s
        V_f_ct(t,t)=V_m(t,x); % �����������٣�sm/s
        Ff_ct(t,x)=Ff_a(t,x); % �������嵥λ����Ħ��ѹ����Pa/m��
        f_ct(t,x)=f_a(t,x); % �������巶��Ħ������
        Re_ct(t,x)=Re_a(t,x); % ����������ŵ����������̬
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % ����������̬
        P_ct(t,x)=P_a(t,x); % ����ѹ����Pa
    end
end

%% ��1��Ntʱ��ڵ��̹ܶγ���ѹ������ѹ���㣨�������ڣ�
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % ��������������٣�m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % ����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % �̹ܶ����嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % �����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % ���ڵ�ע��ͷ������Ħ��ѹ����Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % ��������Ħ��ѹ����Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % ����������Ͳ��Ħ��ѹ����Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % �̹ܶ�Ħ��ѹ����Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % �̹ܶγ���ѹ����Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % ��ѹ��Pa
end

%% ���㾮���ۻ���ɰ��
M_w_tem(1)=0; % ����˲ʱ��ɰ����ֵ��kg
M_w_tot(1)=0; % �����ۻ���ɰ����ֵ��kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮��˲ʱ��ɰ����kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮���ۻ���ɰ����kg
end

%% ����ĩ̬��Ͳ���ɰŨ��
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(1)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% ���㻷����ĭ����ƽ��ֵV_l_mre
V_f_mre=0; % ���շ���ƽ��ֵ��m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % ���շ���ƽ��ֵ��m/s
    end
end

%% ����������ĩ��ƽ��ֵVsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % ��Vsr��������ȡֵ����ĳ���ĩ�٣�m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % �������ĩ��ƽ��ֵ��m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % �������ĩ��ƽ��ֵ��m/s
end

%% �ж��Ƿ���Ч��ɰ�����շ���ƽ��ֵ����2��ɰ������ĩ��ƽ��ֵ�����Ƿ���ɳ�ɰ
fprintf("���϶������1��\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % �����д���˵����ǰ�����¡���ɳ�ɰ��
    else
        fprintf("Sand Cleanout UnFinished!\n"); % �����д���˵����ǰ�����¡�δ��ɳ�ɰ��
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
end

%% ��ͲECD����
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % ����ECD��kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % ����ECD��kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
end

%% ���϶������ĩ״̬����
alpha_f_a_2=alpha_f_a(Nt,:); % ��ϴ�������ĩ״̬��ĭ����
alpha_g_a_2=alpha_g_a(Nt,:); % ��ϴ�������ĩ״̬���ຬ��
alpha_l_a_2=alpha_l_a(Nt,:); % ��ϴ�������ĩ״̬Һ�ຬ��
alpha_s_2=alpha_s(Nt,:); % ��ϴ�������ĩ״̬���ຬ��
f_a_2=f_a(Nt,:); % ��ϴ�������ĩ״̬����Ħ������
Ff_a_2=Ff_a(Nt,:); % ��ϴ�������ĩ״̬���յ�λ����Ħ��ѹ����Pa/m
flow_pattern_a_2=flow_pattern_a(Nt,:); % ��ϴ�������ĩ״̬����������̬
gamma_g_a_2=gamma_g_a(Nt,:); % ��ϴ�������ĩ״̬��ĭ����
gamma_l_a_2=gamma_l_a(Nt,:); % ��ϴ�������ĩ״̬Һ��������
mu_f_a_2=mu_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭճ�ȣ�Pa*s
mu_g_a_2=mu_g_a(Nt,:); % ��ϴ�������ĩ״̬����ճ�ȣ�Pa*s
mu_l_a_2=mu_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ��ճ�ȣ�Pa*s
mu_m_2=mu_m(Nt,:); % ��ϴ�������ĩ״̬���ջ����ճ�ȣ�Pa*s
mu_s_2=mu_s(Nt,:); % ��ϴ�������ĩ״̬����ճ�ȣ�Pa*s
P_a_2=P_a(Nt,:); % ��ϴ�������ĩ״̬����ѹ����Pa
Re_a_2=Re_a(Nt,:); % ��ϴ�������ĩ״̬������ŵ��
rho_f_a_2=rho_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ�ܶȣ�kg/m^3
rho_g_a_2=rho_g_a(Nt,:); % ��ϴ�������ĩ״̬���������ܶȣ�kg/m^3
rho_l_a_2=rho_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ���ܶȣ�kg/m^3
rho_m_2=rho_m(Nt,:); % ��ϴ�������ĩ״̬���ջ�����ܶȣ�kg/m^3
V_f_a_2=V_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ���٣�m/s
V_g_a_2=V_g_a(Nt,:); % ��ϴ�������ĩ״̬�������٣�m/s
V_l_a_2=V_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ�����٣�m/s
V_m_2=V_m(Nt,:); % ��ϴ�������ĩ״̬���ջ�������٣�m/s
V_s_2=V_s(Nt,:); % ��ϴ�������ĩ״̬�������٣�m/s
Va_f_a_2=Va_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ������٣�m/s
Va_g_a_2=Va_g_a(Nt,:); % ��ϴ�������ĩ״̬���������٣�m/s
Va_l_a_2=Va_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ�������٣�m/s
Va_s_2=Va_s(Nt,:); % ��ϴ�������ĩ״̬���������٣�m/s
Vsr_2=Vsr(Nt,:); % ��ϴ�������ĩ״̬���໬���ٶȣ�m/s

%% ���ݴ洢�����϶�����̣�
ANS_Nt_2=Nt; % ʱ��ڵ���
ANS_Time_2=Time(1:Nt); % ���������ʱ�䣬s
ANS_alpha_g_a_2=alpha_g_a(1:Nt,:); % ���������������
ANS_alpha_l_a_2=alpha_l_a(1:Nt,:); % ����Һ���������
ANS_alpha_s_2=alpha_s(1:Nt,:); % ���ຬ��
ANS_P_ct_2=P_ct(1:Nt,:); % ����ѹ����Pa��
ANS_delta_P_SWT_2=delta_P_SWT(1:Nt); % ��ɰ����ѹ����Pa��
ANS_P_a_2=P_a(1:Nt,:); % ����ѹ����Pa��
ANS_T_a_2=T_a(1:Nt,:); % �����¶ȣ��棩
ANS_T_ct_2=T_ct(1:Nt,:); % �����¶ȣ��棩
ANS_P_coil_2=P_coil(1:Nt); % �̹ܶγ���ѹ����Pa��
ANS_P_pump_2=P_pump(1:Nt); % ��ѹ��Pa��
ANS_M_w_tem_2=M_w_tem(1:Nt); % ����˲ʱ��ɰ����kg��
ANS_M_w_tot_2=M_w_tot(1:Nt); % �����ۻ���ɰ����kg��
ANS_Va_s_2=Va_s(1:Nt,:); % ��м����ٶȣ�m/s��
ANS_Va_f_a_2=Va_f_a(1:Nt,:); % ������ĭ������٣�m/s��
ANS_V_s_2=V_s(1:Nt,:); % ��м�����ٶȣ�m/s��
ANS_V_f_a_2=V_f_a(1:Nt,:); % ������ĭ���٣�m/s��
ANS_alpha_g_ct_2=alpha_g_ct(1:Nt,:); % ���������������
ANS_alpha_l_ct_2=alpha_l_ct(1:Nt,:); % ����Һ���������
ANS_gamma_g_ct_2=gamma_g_ct(1:Nt,:); % ������ĭ����
ANS_gamma_g_a_2=gamma_g_a(1:Nt,:); % ������ĭ����
ANS_ECD_a_2=ECD_a(1:Nt,:); % ����ECD��kg/m^3��
ANS_L_coil_2=L_coil(1:Nt); % ���������m��



%% ����ѭ������1�����ʵ�λ�ƣ�
%% ������ֵ����
dt=zeros(); % ʱ�䲽����s
dx=zeros(); % �ռ䲽����m
Time=zeros(); % ������������ʱ����s
L_coil=zeros(); % ���������m
L_reel=zeros(); % �̹ܶγ��ȣ�m
D_t_i=zeros(); % �͹��ھ���m
D_ct_o=zeros(); % �����͹��⾶��m
D_ct_i=zeros(); % �����͹��ھ���m
M_s=zeros(); % ���׽�ɰ����kg/s
P_a=zeros(); % ����ѹ����Pa
rho_g_a=zeros(); % �����ܶȣ�kg/m^3
rho_l_a=zeros(); % ����Һ���ܶȣ�kg/m^3
rho_f_a=zeros(); % ������ĭ�ܶȣ�kg/m^3
mu_g_a=zeros(); % ����ճ�ȣ�Pa*s
mu_l_a=zeros(); % ����Һ��ճ�ȣ�Pa*s
mu_f_a=zeros(); % ������ĭճ�ȣ�Pa*s
mu_s=zeros(); % ����ճ�ȣ�Pa*s
alpha_g_a=zeros(); % ���ຬ��
alpha_l_a=zeros(); % Һ�ຬ��
alpha_f_a=zeros(); % ��ĭ����
alpha_s=zeros(); % ���ຬ��
Va_g_a=zeros(); % ���������٣�m/s
Va_l_a=zeros(); % ����Һ�������٣�m/s
Va_f_a=zeros(); % ������ĭ������٣�m/s
Va_s=zeros(); % ���������٣�m/s
V_g_a=zeros(); % �������٣�m/s
V_l_a=zeros(); % ����Һ�����٣�m/s
V_f_a=zeros(); % ������ĭ���٣�m/s
V_s=zeros(); % �������٣�m/s
Vsr=zeros(); % ���໬���ٶȣ�m/s
gamma_g_a=zeros(); % ��ĭ����
gamma_l_a=zeros(); % Һ��������
V_m=zeros(); % ���ջ�������٣�m/s
rho_m=zeros(); % ���ջ�����ܶȣ�kg/m^3
mu_m=zeros(); % ���ջ����ճ�ȣ�Pa*s
Ff_a=zeros(); % ���յ�λ����Ħ��ѹ����Pa/m
f_a=zeros(); % ����Ħ������
Re_a=zeros(); % ������ŵ��
flow_pattern_a=zeros(); % ����������̬

%% ����ռ䲽��dx����Ӧʱ�䲽��dt
t_3=0.5*3600; % ����ѭ����ʱ����s

Nt=101; % ʱ��ڵ��������룩
nt=Nt-1; % ʱ��������

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % �ռ䲽����m
end

for t=1:1:Nt-1
    dt(t)=t_3/nt; % ʱ�䲽��
end

Time(1)=0; % ����ѭ����ʱ����s
for t=1:1:Nt-1
    Time(t+1)=Time(t)+dt(t); % ����ѭ������(t+1)��ʱ������������ʱ����s
end

%% ���㲻ͬʱ���������������L_coil���̹ܶ������ܳ���L_reel
L=10000; % �����͹��ܳ���m�����룩
L_wg=8; % ���ڵ�ע��ͷ�����������ܳ��ȣ�m�����룩
L_goose=3; % �������������ܳ��ȣ�m�����룩
D_goose=2; % �������ΰ뾶��m�����룩
H_goose=10; % �����������߶ȣ�m�����룩
L_gr=20; % ����������Ͳ�������ܳ��ȣ�m�����룩
theta_gr=acosd(H_goose/L_gr); % ����������Ͳ����������Ǧ���߼нǣ���
D_r_i=3; % ��Ͳ�ھ���m�����룩
D_r_o=5; % ��Ͳ�⾶��m�����룩
W_r=5; % ��Ͳ��ȣ�m�����룩
D_cable=0.005; % �����⾶��m

L_coil(1)=Lp; % ������������ȣ�m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % ��ʼʱ���̹ܶ������ܳ��ȣ�m
for t=2:1:Nt
    L_coil(t)=Lp; % ������������ȣ�m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % �̹ܶ������ܳ��ȣ�m
end
L_cable=L_coil; % �������볤�ȣ�m

%% �������뼰Ԥ����
D_ct_o_0=0.04445; % �����͹��⾶��m�����룩
L1=2000; % �����ڵ�һ�������ܳ��ȣ�m�����룩
D_ct_i_1=0.03709; % �����ڵ�һ���������ھ���m�����룩
L2=2000; % �����ڵڶ��������ܳ��ȣ�m�����룩
D_ct_i_2=0.03653; % �����ڵڶ����������ھ���m�����룩
L3=2000; % �����ڵ����������ܳ��ȣ�m�����룩
D_ct_i_3=0.03555; % �����ڵ������������ھ���m�����룩
L4=2000; % �����ڵ��Ķ������ܳ��ȣ�m�����룩
D_ct_i_4=0.03489; % �����ڵ��Ķ��������ھ���m�����룩
L5=2000; % �����ڵ���������ܳ��ȣ�m�����룩
D_ct_i_5=0.03409; % �����ڵ�����������ھ���m�����룩

L_t_1=4000; % �ϲ��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_1=0.09718; %0.068;% �ϲ��͹ܣ����׹ܻ����ۣ��ھ���m�����룩
L_t_2=2200; % �²��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_2=0.09718; %0.13970; % �²��͹ܣ����׹ܻ����ۣ��ھ���m�����룩

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % �͹ܣ����׹ܻ����ۣ��ھ���m
        else
            D_t_i(t,x)=D_t_i_2; % �͹ܣ����׹ܻ����ۣ��ھ���m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            D_ct_o(t,x)=D_ct_o_0; % �����͹��⾶��m
        else
            D_ct_o(t,x)=0; % �����͹��⾶��m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            if (L_coil(t)-Depth(x))<=L1 && (L_coil(t)-Depth(x))>=0
                D_ct_i(t,x)=D_ct_i_1; % �����͹��ھ���m
            elseif (L_coil(t)-Depth(x))<=L1+L2
                D_ct_i(t,x)=D_ct_i_2; % �����͹��ھ���m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                D_ct_i(t,x)=D_ct_i_3; % �����͹��ھ���m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                D_ct_i(t,x)=D_ct_i_4; % �����͹��ھ���m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                D_ct_i(t,x)=D_ct_i_5; % �����͹��ھ���m
            end
        else
            D_ct_i(t,x)=0; % �����͹��ھ���m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % �����ļ������������ֵ
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % �����͹��ڽ������m^2
        else
            A_ct(t,x)=0; % �����͹��ڽ������m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % ����ˮ��ֱ����m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % ���ս������m^2
    end
end

%% �������뼰Ԥ����
h_t=25.4*10^(-6); % �͹ܣ����׹ܻ����ۣ����Դֲڶȣ�m�����룩
h_ct=25.4*10^(-6); % �����͹ܾ��Դֲڶȣ�m�����룩
h_a=(h_t+h_ct)/2; % ����ƽ�����Դֲڶȣ�m
epsilon_e=1*10^(-3); % �����������ޣ����룩
epsilon_t=1*10^3; % ���������������룩
g=9.81; % �������ٶȣ�m/s^2��Ĭ�ϣ�

T_0=20; % ��Һ�����¶ȣ��棨���룩
P_0=0.1*10^6; % ��Һ����ѹ����Pa�����룩
rho_l_0=1150; % T_0��P_0�»�Һ�ܶȣ�kg/m^3�����룩
mu_l_0=0.03; % T_0��P_0�»�Һճ�ȣ�Pa*s�����룩
Qv_l_0=0.1/60; % ��Һ���������m^3/s�����룩
Qm_l_0=Qv_l_0*rho_l_0; % ��Һ����������kg/s
rho_g_0=0.655; % ע�����ܶȣ�kg/m^3
mu_g_0=RheologyG(T_0,P_0); % ����ճ�ȣ�Pa*s
Qv_g_0=4/60; % ע�������������m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % ע��������������kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % ��ĭ����������kg/s

D_s=1*10^(-3); % ɰ��ֱ����m�����룩
rho_s=2000; % ɰ���ܶȣ�kg/m^3�����룩
H_s=L_s_b-Dsp_pen(1); % �ײ�ɰ���߶ȣ�m
PHI=0.6; % ɰ����ӯ�ȣ����룩
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % ����ɰ����������kg

D_nozzle=4/1000; % ����ֱ����m�����룩
N_nozzle=3; % ������������룩
C=0.95; % ��������ϵ����ȡ0.95�����룩

C0=1.2; % Ư��������ϵ����Ĭ�ϣ�

M_s(1)=0; % ���׽�ɰ����kg/s
for t=2:1:Nt
    M_s(t)=0; % ���׽�ɰ����kg/s
end

OutPressure=1*10^6; % ����ѹ����Pa�����룩

%% �¶�����
T_i=20; % ��ĭע���¶ȣ���
T_g=0.02; % �����ݶȣ���/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % ����������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������
    end
end
T_a=T_ct; % ������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������

%% ��1��ʱ��ڵ㣨��ʼʱ�̣���ز�����ֵ���㣨���գ�
for x=1:1:Nx
    P_a(1,x)=P_a_2(x); % ����ѹ����Pa
    rho_g_a(1,x)=rho_g_a_2(x); % �����ܶȣ�kg/m^3
    rho_l_a(1,x)=rho_l_a_2(x); % ����Һ���ܶȣ�kg/m^3
    rho_f_a(1,x)=rho_f_a_2(x); % ������ĭ�ܶȣ�kg/m^3
    mu_g_a(1,x)=mu_g_a_2(x); % ����ճ�ȣ�Pa*s
    mu_l_a(1,x)=mu_l_a_2(x); % ����Һ��ճ�ȣ�Pa*s
    mu_f_a(1,x)=mu_f_a_2(x); % ������ĭճ�ȣ�Pa*s
    mu_s(1,x)=mu_s_2(x); % ����ճ�ȣ�Pa*s
    alpha_g_a(1,x)=alpha_g_a_2(x); % ���ຬ��
    alpha_l_a(1,x)=alpha_l_a_2(x); % Һ�ຬ��
    alpha_f_a(1,x)=alpha_f_a_2(x); % ��ĭ����
    alpha_s(1,x)=alpha_s_2(1,x); % ���ຬ��
    Va_g_a(1,x)=Va_g_a_2(x); % ���������٣�m/s
    Va_l_a(1,x)=Va_l_a_2(x); % ����Һ�������٣�m/s
    Va_f_a(1,x)=Va_f_a_2(x); % ������ĭ������٣�m/s
    Va_s(1,x)=Va_s_2(x); % ���������٣�m/s
    V_g_a(1,x)=V_g_a_2(x); % �������٣�m/s
    V_l_a(1,x)=V_l_a_2(x); % ����Һ�����٣�m/s
    V_f_a(1,x)=V_f_a_2(x); % ������ĭ���٣�m/s
    V_s(1,x)=V_s_2(x); % �������٣�m/s
    Vsr(1,x)=Vsr_2(x); % ���໬���ٶȣ�m/s
    gamma_g_a(1,x)=gamma_g_a_2(x); % ��ĭ����
    gamma_l_a(1,x)=gamma_l_a_2(1,x); % Һ��������
    V_m(1,x)=V_m_2(x); % ���ջ�������٣�m/s
    rho_m(1,x)=rho_m_2(x); % ���ջ�����ܶȣ�kg/m^3
    mu_m(1,x)=mu_m_2(x); % ���ջ����ճ�ȣ�Pa*s
    Ff_a(1,x)=Ff_a_2(x); % ���յ�λ����Ħ��ѹ����Pa/m
    f_a(1,x)=f_a_2(x); % ����Ħ������
    Re_a(1,x)=Re_a_2(x); % ������ŵ��
    flow_pattern_a(1,x)=flow_pattern_a_2(x); % ����������̬
end

%% ��2��Nt��ʱ��ڵ���ز������㣨���գ�
for t=2:1:Nt
    P_a(t,Nx_Lp)=P_a(t-1,Nx_Lp);  % ���չܵ�ѹ������ֵ��Pa
    
    err_OutPressure=1; % ����ѹ��������
    COUNT_OutPressure=0; % ����ѹ������������ֵ
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % ����ѹ����������
        
        % �ܵף���Nx_Lp���ռ�ڵ㣩����ز�������
        rho_g_a(t,Nx_Lp)=DensityG(T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % ���������ܶȣ�kg/m^3
        rho_l_a(t,Nx_Lp)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % ���ջ�Һ�ܶȣ�kg/m^3
        mu_g_a(t,Nx_Lp)=RheologyG(T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % ��������ճ�ȣ�Pa*s
        mu_l_a(t,Nx_Lp)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % ���ջ�Һճ�ȣ�Pa*s
        gamma_g_a(t,Nx_Lp)=(Qm_g_0/rho_g_a(t,Nx_Lp))/(Qm_g_0/rho_g_a(t,Nx_Lp)+Qm_l_0/rho_l_a(t,Nx_Lp)); % ��ĭ����
        gamma_l_a(t,Nx_Lp)=1-gamma_g_a(t,Nx_Lp); % Һ��������
        alpha_g_a(t,Nx_Lp)=alpha_f_a(t-1,Nx_Lp)*gamma_g_a(t,Nx_Lp); % �������庬��
        alpha_l_a(t,Nx_Lp)=alpha_f_a(t-1,Nx_Lp)*gamma_l_a(t,Nx_Lp); % ���ջ�Һ����
        rho_f_a(t,Nx_Lp)=rho_g_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp)+rho_l_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % ������ĭ�ܶȣ�kg/m^3
        mu_f_a(t,Nx_Lp)=mu_g_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp)+mu_l_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % ������ĭճ�ȣ�Pa*s
        mu_s(t,Nx_Lp)=mu_f_a(t,Nx_Lp); % ����ճ�ȣ�Pa*s
        Va_f_a(t,Nx_Lp)=Qm_f_0/(A_a(t,Nx_Lp)*rho_f_a(t,Nx_Lp)); % ������ĭ������٣�m/s
        Va_g_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp); % �������������٣�m/s
        Va_l_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp); % ���ջ�Һ������٣�m/s
        Va_s(t,Nx_Lp)=M_s(t)/(A_a(t,Nx_Lp)*rho_s); % ��м����ٶȣ�m/s
        Vsr(t,Nx_Lp)=12*(mu_f_a(t,Nx_Lp)/(rho_f_a(t,Nx_Lp)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Lp))/rho_f_a(t,Nx_Lp))*((rho_f_a(t,Nx_Lp)*D_s/mu_f_a(t,Nx_Lp))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
        alpha_s(t,Nx_Lp)=Va_s(t,Nx_Lp)/(C0*(Va_s(t,Nx_Lp)+Va_f_a(t,Nx_Lp))-Vsr(t,Nx_Lp));  % �����������
        V_s(t,Nx_Lp)=0;%Va_s(t,Nx_Lp)/alpha_s(t,Nx_Lp); % ��м�ٶȣ�m/s
        
        alpha_f_a(t,Nx_Lp)=1-alpha_s(t,Nx_Lp); % ������ĭ����
        alpha_g_a(t,Nx_Lp)=alpha_f_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp); % �������庬��
        alpha_l_a(t,Nx_Lp)=alpha_f_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % ���ջ�Һ����
        V_f_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp)/alpha_f_a(t,Nx_Lp); % ������ĭ���٣�m/s
        V_g_a(t,Nx_Lp)=V_f_a(t,Nx_Lp); % �����������٣�m/s
        V_l_a(t,Nx_Lp)=V_f_a(t,Nx_Lp); % ���ջ�Һ���٣�m/s
        
        V_m(t,Nx_Lp)=Va_s(t,Nx_Lp)+Va_f_a(t,Nx_Lp); % ���ջ�����ٶȣ�m/s
        rho_m(t,Nx_Lp)=alpha_s(t,Nx_Lp)*rho_s+alpha_f_a(t,Nx_Lp)*rho_f_a(t,Nx_Lp); % ���ջ�����ܶȣ�kg/m^3
        mu_m(t,Nx_Lp)=alpha_s(t,Nx_Lp)*mu_s(t,Nx_Lp)+alpha_f_a(t,Nx_Lp)*mu_f_a(t,Nx_Lp); % ���ջ����ճ�ȣ�Pa*s
        [Ff_a(t,Nx_Lp),f_a(t,Nx_Lp),Re_a(t,Nx_Lp),flow_pattern_a(t,Nx_Lp)]=Friction_annulus(rho_m(t,Nx_Lp),V_m(t,Nx_Lp),mu_m(t,Nx_Lp),D_h(t,Nx_Lp),h_a,rho_f_a(t,Nx_Lp),V_f_a(t,Nx_Lp)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
        
        % ��Nx_Lp-1��1���ռ�ڵ㴦��ز�������
        for x=Nx_Lp-1:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % ����ѹ������ֵ��Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % ���ջ�Һ����
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            
            err_NodePressure=1; % ����ѹ��������
            COUNT_NodePressure=0; % ����ѹ������������ֵ
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % ����ѹ����������
                
                % ���߷��������������
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % ���������������ֵ1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % ���������������ֵ2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % ������������������
                COUNT_NodeEg=0; % ��������
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % �����������Ϊalpha_s_ass1ʱ������ֵ��
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % �����ٶȣ�m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % ���������٣�m/s
                    alpha_f_ass1=1-alpha_s_ass1; % ��ĭ�������
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % ��ĭ�ٶȣ�m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % ��ĭ������٣�m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % ���ջ�����ٶȣ�m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % ���������������ֵ
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % ����ĺ��������Ľ������ʵ�����������
                    
                    % �����������Ϊalpha_s_ass2ʱ������ֵ��
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % ���������غ㷽����ɢ��ʽ�м�ֵ����
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % �����ٶȣ�m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % ���������٣�m/s
                    alpha_f_ass2=1-alpha_s_ass2; % ��ĭ�������
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % ��ĭ�ٶȣ�m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % ��ĭ������٣�m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % ���ջ�����ٶȣ�m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % ���������������ֵ
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % ����ĺ��������Ľ������ʵ�����������
                    
                    % ���߷���������������
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % �µĹ��������������ֵalpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % ������������������
                    alpha_s_ass1=alpha_s_ass2; % �µĹ��������������ֵ1
                    alpha_s_ass2=alpha_s_ass3; % �µĹ��������������ֵ2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % ���������õ�����ʵ�����������ֵ����alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % �������������С��һ��ֵʱ����Ϊ���������Ϊ0�����ڷ�ֹ���ֺ����ļ������
                    alpha_s(t,x)=0; % �����������
                    V_s(t,x)=0; % �����ٶȣ�m/s
                    Va_s(t,x)=0; % �������ٶȣ�m/s
                    mu_s(t,x)=0; % ����ճ�ȣ�Pa*s
                    Vsr(t,x)=0; % ��м����ĩ�٣�m/s
                    
                    alpha_f_a(t,x)=1; % ��ĭ�������                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                else
                    V_s(t,x)=V_s_ass2; % �����ٶȣ�m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % ���������٣�m/s
                    Vsr(t,x)=Vsr_ass2; % ��м����ĩ�٣�m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % ��ĭ�ܶȣ�kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % ��ĭ�������
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    V_f_a(t,x)=V_f_ass2; % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % ���ջ�����ٶȣ�m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % ����ѹ������ֵ��Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % ����ѹ��������
                P_a(t,x)=P_new; % �µĻ���ѹ������ֵ��Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % ����ѹ��������        
        if (P_a(t,1)-OutPressure)>0 % ���ݳ���ѹ�������������Ի��չܵ�ѹ������ֵ���е���
            P_a(t,Nx_Lp)=P_a(t,Nx_Lp)-(P_a(t,1)-OutPressure)/2; % �µĻ��չܵ�ѹ������ֵ��Pa
        else
            P_a(t,Nx_Lp)=P_a(t,Nx_Lp)-(P_a(t,1)-OutPressure)/2*0.3; % �µĻ��չܵ�ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Lp+1��Nx_Dsp_pen(1)���ռ�ڵ㴦��ز�������
    for x=Nx_Lp+1:1:Nx_Dsp_pen(1)
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=0; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(1)+1��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)+1:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1-PHI; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=PHI; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
end

%% ��1��Ntʱ��ڵ���ͷѹ������
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Lp); % ������ĭ���������m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Lp)); % �����ɰ����ѹ����Pa����������������V_nozzle��m/s��
end

%% ��1��Ntʱ��ڵ���ز������㣨�������ڣ�
for t=1:1:Nt
    % �ܵף���Nx_Lp���ռ�ڵ㣩����ز�������
    P_ct(t,Nx_Lp)=P_a(t,Nx_Lp)+delta_P_SWT(t); % ����ѹ����Pa
    rho_g_ct(t,Nx_Lp)=DensityG(T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % ���������ܶȣ�kg/m^3
    rho_l_ct(t,Nx_Lp)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % ���ڻ�Һ�ܶȣ�kg/m^3
    mu_g_ct(t,Nx_Lp)=RheologyG(T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % ��������ճ�ȣ�Pa*s
    mu_l_ct(t,Nx_Lp)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % ���ڻ�Һճ�ȣ�Pa*s
    alpha_f_ct(t,Nx_Lp)=1; % ������ĭ����
    gamma_g_ct(t,Nx_Lp)=(Qm_g_0/rho_g_ct(t,Nx_Lp))/(Qm_g_0/rho_g_ct(t,Nx_Lp)+Qm_l_0/rho_l_ct(t,Nx_Lp)); % ������ĭ����
    gamma_l_ct(t,Nx_Lp)=1-gamma_g_ct(t,Nx_Lp); % ����Һ��������
    alpha_g_ct(t,Nx_Lp)=alpha_f_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp); % �������庬��
    alpha_l_ct(t,Nx_Lp)=alpha_f_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % ���ڻ�Һ����
    rho_f_ct(t,Nx_Lp)=rho_g_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp)+rho_l_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % ������ĭ�ܶȣ�kg/m^3
    mu_f_ct(t,Nx_Lp)=mu_g_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp)+mu_l_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % ������ĭճ�ȣ�Pa*s
    V_f_ct(t,Nx_Lp)=Qm_f_0/(rho_f_ct(t,Nx_Lp)*A_ct(t,Nx_Lp)); % �����������٣�m/s
    [Ff_ct(t,Nx_Lp),f_ct(t,Nx_Lp),Re_ct(t,Nx_Lp),flow_pattern_ct(t,Nx_Lp)]=Friction_coiledtubing(rho_f_ct(t,Nx_Lp),V_f_ct(t,Nx_Lp),mu_f_ct(t,Nx_Lp),D_ct_i(t,Nx_Lp),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    % ��Nx_Lp-1��1���ռ�ڵ㴦��ز�������
    for x=Nx_Lp-1:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % ����ѹ������ֵ��Pa
        
        err_DriPipePressure=1; % ����ѹ��������
        COUNT_DriPipePressure=0; % ����ѹ������������ֵ
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % ����ѹ����������
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һ�ܶȣ�kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һճ�ȣ�Pa*s
            alpha_f_ct(t,x)=1; % ������ĭ����
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % ������ĭ����
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % �������庬��
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % ���ڻ�Һ����
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭճ�ȣ�Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % ������ĭ���٣�m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % ����ѹ����Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % �������ѹ������ֵ�����ֵ֮���������
            P_ct_ass(t,x)=P_ct(t,x); % �µĹ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Lp+1��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Lp+1:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % �������庬��
        alpha_l_ct(t,x)=alpha_l_a(t,x); % ���ڻ�Һ����
        gamma_g_ct(t,x)=gamma_g_a(t,x); % ������ĭ����
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
        rho_f_ct(t,x)=rho_m(t,x); % ���������ܶȣ�kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % ��������ճ�ȣ�Pa*s
        V_f_ct(t,t)=V_m(t,x); % �����������٣�sm/s
        Ff_ct(t,x)=Ff_a(t,x); % �������嵥λ����Ħ��ѹ����Pa/m��
        f_ct(t,x)=f_a(t,x); % �������巶��Ħ������
        Re_ct(t,x)=Re_a(t,x); % ����������ŵ����������̬
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % ����������̬
        P_ct(t,x)=P_a(t,x); % ����ѹ����Pa
    end
end

%% ��1��Ntʱ��ڵ��̹ܶγ���ѹ������ѹ���㣨�������ڣ�
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % ��������������٣�m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % ����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % �̹ܶ����嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % �����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % ���ڵ�ע��ͷ������Ħ��ѹ����Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % ��������Ħ��ѹ����Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % ����������Ͳ��Ħ��ѹ����Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % �̹ܶ�Ħ��ѹ����Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % �̹ܶγ���ѹ����Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % ��ѹ��Pa
end

%% ���㾮���ۻ���ɰ��
M_w_tem(1)=0; % ����˲ʱ��ɰ����ֵ��kg
M_w_tot(1)=0; % �����ۻ���ɰ����ֵ��kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮��˲ʱ��ɰ����kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮���ۻ���ɰ����kg
end

%% ����ĩ̬��Ͳ���ɰŨ��
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(1)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% ���㻷����ĭ����ƽ��ֵV_l_mre
V_f_mre=0; % ���շ���ƽ��ֵ��m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % ���շ���ƽ��ֵ��m/s
    end
end

%% ����������ĩ��ƽ��ֵVsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % ��Vsr��������ȡֵ����ĳ���ĩ�٣�m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % �������ĩ��ƽ��ֵ��m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % �������ĩ��ƽ��ֵ��m/s
end

%% �ж��Ƿ���Ч��ɰ�����շ���ƽ��ֵ����2��ɰ������ĩ��ƽ��ֵ�����Ƿ���ɳ�ɰ
fprintf("����ѭ������1��\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % �����д���˵����ǰ�����¡���ɳ�ɰ��
    else
        fprintf("Sand Cleanout UnFinished!\n"); % �����д���˵����ǰ�����¡�δ��ɳ�ɰ��
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
end

%% ��ͲECD����
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % ����ECD��kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % ����ECD��kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
end

%% ����ѭ������ĩ״̬����
alpha_f_a_3=alpha_f_a(Nt,:); % ����ѭ������ĩ״̬��ĭ����
alpha_g_a_3=alpha_g_a(Nt,:); % ����ѭ������ĩ״̬���ຬ��
alpha_l_a_3=alpha_l_a(Nt,:); % ����ѭ������ĩ״̬Һ�ຬ��
alpha_s_3=alpha_s(Nt,:); % ����ѭ������ĩ״̬���ຬ��
f_a_3=f_a(Nt,:); % ����ѭ������ĩ״̬����Ħ������
Ff_a_3=Ff_a(Nt,:); % ����ѭ������ĩ״̬���յ�λ����Ħ��ѹ����Pa/m
flow_pattern_a_3=flow_pattern_a(Nt,:); % ����ѭ������ĩ״̬����������̬
gamma_g_a_3=gamma_g_a(Nt,:); % ����ѭ������ĩ״̬��ĭ����
gamma_l_a_3=gamma_l_a(Nt,:); % ����ѭ������ĩ״̬Һ��������
mu_f_a_3=mu_f_a(Nt,:); % ����ѭ������ĩ״̬������ĭճ�ȣ�Pa*s
mu_g_a_3=mu_g_a(Nt,:); % ����ѭ������ĩ״̬����ճ�ȣ�Pa*s
mu_l_a_3=mu_l_a(Nt,:); % ����ѭ������ĩ״̬����Һ��ճ�ȣ�Pa*s
mu_m_3=mu_m(Nt,:); % ����ѭ������ĩ״̬���ջ����ճ�ȣ�Pa*s
mu_s_3=mu_s(Nt,:); % ����ѭ������ĩ״̬����ճ�ȣ�Pa*s
P_a_3=P_a(Nt,:); % ����ѭ������ĩ״̬����ѹ����Pa
Re_a_3=Re_a(Nt,:); % ����ѭ������ĩ״̬������ŵ��
rho_f_a_3=rho_f_a(Nt,:); % ����ѭ������ĩ״̬������ĭ�ܶȣ�kg/m^3
rho_g_a_3=rho_g_a(Nt,:); % ����ѭ������ĩ״̬���������ܶȣ�kg/m^3
rho_l_a_3=rho_l_a(Nt,:); % ����ѭ������ĩ״̬����Һ���ܶȣ�kg/m^3
rho_m_3=rho_m(Nt,:); % ����ѭ������ĩ״̬���ջ�����ܶȣ�kg/m^3
V_f_a_3=V_f_a(Nt,:); % ����ѭ������ĩ״̬������ĭ���٣�m/s
V_g_a_3=V_g_a(Nt,:); % ����ѭ������ĩ״̬�������٣�m/s
V_l_a_3=V_l_a(Nt,:); % ����ѭ������ĩ״̬����Һ�����٣�m/s
V_m_3=V_m(Nt,:); % ����ѭ������ĩ״̬���ջ�������٣�m/s
V_s_3=V_s(Nt,:); % ����ѭ������ĩ״̬�������٣�m/s
Va_f_a_3=Va_f_a(Nt,:); % ����ѭ������ĩ״̬������ĭ������٣�m/s
Va_g_a_3=Va_g_a(Nt,:); % ����ѭ������ĩ״̬���������٣�m/s
Va_l_a_3=Va_l_a(Nt,:); % ����ѭ������ĩ״̬����Һ�������٣�m/s
Va_s_3=Va_s(Nt,:); % ����ѭ������ĩ״̬���������٣�m/s
Vsr_3=Vsr(Nt,:); % ����ѭ������ĩ״̬���໬���ٶȣ�m/s

%% ���ݴ洢������ѭ����
ANS_Nt_3=Nt; % ʱ��ڵ���
ANS_Time_3=Time(1:Nt); % ����ѭ��ʱ�䣬s
ANS_alpha_g_a_3=alpha_g_a(1:Nt,:); % ���������������
ANS_alpha_l_a_3=alpha_l_a(1:Nt,:); % ����Һ���������
ANS_alpha_s_3=alpha_s(1:Nt,:); % ���ຬ��
ANS_P_ct_3=P_ct(1:Nt,:); % ����ѹ����Pa��
ANS_delta_P_SWT_3=delta_P_SWT(1:Nt); % ��ɰ����ѹ����Pa��
ANS_P_a_3=P_a(1:Nt,:); % ����ѹ����Pa��
ANS_T_a_3=T_a(1:Nt,:); % �����¶ȣ��棩
ANS_T_ct_3=T_ct(1:Nt,:); % �����¶ȣ��棩
ANS_P_coil_3=P_coil(1:Nt); % �̹ܶγ���ѹ����Pa��
ANS_P_pump_3=P_pump(1:Nt); % ��ѹ��Pa��
ANS_M_w_tem_3=M_w_tem(1:Nt); % ����˲ʱ��ɰ����kg��
ANS_M_w_tot_3=M_w_tot(1:Nt); % �����ۻ���ɰ����kg��
ANS_Va_s_3=Va_s(1:Nt,:); % ��м����ٶȣ�m/s��
ANS_Va_f_a_3=Va_f_a(1:Nt,:); % ������ĭ������٣�m/s��
ANS_V_s_3=V_s(1:Nt,:); % ��м�����ٶȣ�m/s��
ANS_V_f_a_3=V_f_a(1:Nt,:); % ������ĭ���٣�m/s��
ANS_alpha_g_ct_3=alpha_g_ct(1:Nt,:); % ���������������
ANS_alpha_l_ct_3=alpha_l_ct(1:Nt,:); % ����Һ���������
ANS_gamma_g_ct_3=gamma_g_ct(1:Nt,:); % ������ĭ����
ANS_gamma_g_a_3=gamma_g_a(1:Nt,:); % ������ĭ����
ANS_ECD_a_3=ECD_a(1:Nt,:); % ����ECD��kg/m^3��
ANS_L_coil_3=L_coil(1:Nt); % ���������m��



%% �������������1�����ʵ�λ�ƣ�
%% ������ֵ����
dt=zeros(); % ʱ�䲽����s
dx=zeros(); % �ռ䲽����m
Time=zeros(); % ������������ʱ����s
L_coil=zeros(); % ���������m
L_reel=zeros(); % �̹ܶγ��ȣ�m
D_t_i=zeros(); % �͹��ھ���m
D_ct_o=zeros(); % �����͹��⾶��m
D_ct_i=zeros(); % �����͹��ھ���m
M_s=zeros(); % ���׽�ɰ����kg/s
P_a=zeros(); % ����ѹ����Pa
rho_g_a=zeros(); % �����ܶȣ�kg/m^3
rho_l_a=zeros(); % ����Һ���ܶȣ�kg/m^3
rho_f_a=zeros(); % ������ĭ�ܶȣ�kg/m^3
mu_g_a=zeros(); % ����ճ�ȣ�Pa*s
mu_l_a=zeros(); % ����Һ��ճ�ȣ�Pa*s
mu_f_a=zeros(); % ������ĭճ�ȣ�Pa*s
mu_s=zeros(); % ����ճ�ȣ�Pa*s
alpha_g_a=zeros(); % ���ຬ��
alpha_l_a=zeros(); % Һ�ຬ��
alpha_f_a=zeros(); % ��ĭ����
alpha_s=zeros(); % ���ຬ��
Va_g_a=zeros(); % ���������٣�m/s
Va_l_a=zeros(); % ����Һ�������٣�m/s
Va_f_a=zeros(); % ������ĭ������٣�m/s
Va_s=zeros(); % ���������٣�m/s
V_g_a=zeros(); % �������٣�m/s
V_l_a=zeros(); % ����Һ�����٣�m/s
V_f_a=zeros(); % ������ĭ���٣�m/s
V_s=zeros(); % �������٣�m/s
Vsr=zeros(); % ���໬���ٶȣ�m/s
gamma_g_a=zeros(); % ��ĭ����
gamma_l_a=zeros(); % Һ��������
V_m=zeros(); % ���ջ�������٣�m/s
rho_m=zeros(); % ���ջ�����ܶȣ�kg/m^3
mu_m=zeros(); % ���ջ����ճ�ȣ�Pa*s
Ff_a=zeros(); % ���յ�λ����Ħ��ѹ����Pa/m
f_a=zeros(); % ����Ħ������
Re_a=zeros(); % ������ŵ��
flow_pattern_a=zeros(); % ����������̬

%% ����ռ䲽��dx����Ӧʱ�䲽��dt
V4=(Dsp_pen(1)-Lp)/(1*3600); % �����͹������ٶȣ�m/s�����룩
t_4=(Dsp_pen(1)-Lp)/V4; % �����͹�������ʱ����s

Nt=Nx_Dsp_pen(1)-Nx_Lp+1; % ʱ��ڵ���
nt=Nt-1; % ʱ��������

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % �ռ䲽����m
end

Time(1)=0; % ������������ʱ����ֵ��s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Lp+t-1)/V4; % ÿ����һ���ռ䲽������ʱ�䣬s
    Time(t+1)=Time(t)+dt(t); % ��������������(t+1)���ռ�ڵ�����������ʱ����s
end

%% ���㲻ͬʱ���������������L_coil���̹ܶ������ܳ���L_reel
L=10000; % �����͹��ܳ���m�����룩
L_wg=8; % ���ڵ�ע��ͷ�����������ܳ��ȣ�m�����룩
L_goose=3; % �������������ܳ��ȣ�m�����룩
D_goose=2; % �������ΰ뾶��m�����룩
H_goose=10; % �����������߶ȣ�m�����룩
L_gr=20; % ����������Ͳ�������ܳ��ȣ�m�����룩
theta_gr=acosd(H_goose/L_gr); % ����������Ͳ����������Ǧ���߼нǣ���
D_r_i=3; % ��Ͳ�ھ���m�����룩
D_r_o=5; % ��Ͳ�⾶��m�����룩
W_r=5; % ��Ͳ��ȣ�m�����룩
D_cable=0; % �����⾶��m

L_coil(1)=Lp; % ��ʼʱ��������������ȣ�m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % ��ʼʱ���̹ܶ������ܳ��ȣ�m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)+dx(Nx_Lp+t-2); % ������������ȣ�m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % �̹ܶ������ܳ��ȣ�m
end
L_cable=L_coil; % �������볤�ȣ�m

%% �������뼰Ԥ����
D_ct_o_0=0.04445; % �����͹��⾶��m�����룩
L1=2000; % �����ڵ�һ�������ܳ��ȣ�m�����룩
D_ct_i_1=0.03709; % �����ڵ�һ���������ھ���m�����룩
L2=2000; % �����ڵڶ��������ܳ��ȣ�m�����룩
D_ct_i_2=0.03653; % �����ڵڶ����������ھ���m�����룩
L3=2000; % �����ڵ����������ܳ��ȣ�m�����룩
D_ct_i_3=0.03555; % �����ڵ������������ھ���m�����룩
L4=2000; % �����ڵ��Ķ������ܳ��ȣ�m�����룩
D_ct_i_4=0.03489; % �����ڵ��Ķ��������ھ���m�����룩
L5=2000; % �����ڵ���������ܳ��ȣ�m�����룩
D_ct_i_5=0.03409; % �����ڵ�����������ھ���m�����룩

L_t_1=4000; % �ϲ��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_1=0.09718; %0.068;%0.09718; % �ϲ��͹ܣ����׹ܻ����ۣ��ھ���m�����룩
L_t_2=2200; % �²��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_2=0.09718; %0.13970; % �²��͹ܣ����׹ܻ����ۣ��ھ���m�����룩

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % �͹ܣ����׹ܻ����ۣ��ھ���m
        else
            D_t_i(t,x)=D_t_i_2; % �͹ܣ����׹ܻ����ۣ��ھ���m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Dsp_pen(1)
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % �����͹��⾶��m
            else
                D_ct_o(t,x)=0; % �����͹��⾶��m
            end
        else
            D_ct_o(t,x)=0; % �����͹��⾶��m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Dsp_pen(1)
            if Depth(x)<=L_coil(t) % �����ļ�����������ھ�
                if (L_coil(t)-Depth(x))<=L1
                    D_ct_i(t,x)=D_ct_i_1; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % �����͹��ھ���m
                end
            else
                D_ct_i(t,x)=0; % �����͹��ھ���m
            end
        else
            D_ct_i(t,x)=0; % �����͹��ھ���m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % �����ļ������������ֵ
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % �����͹��ڽ������m^2
        else
            A_ct(t,x)=0; % �����͹��ڽ������m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % ����ˮ��ֱ����m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % ���ս������m^2
    end
end

%% �������뼰Ԥ����
h_t=25.4*10^(-6); % �͹ܣ����׹ܻ����ۣ����Դֲڶȣ�m�����룩
h_ct=25.4*10^(-6); % �����͹ܾ��Դֲڶȣ�m�����룩
h_a=(h_t+h_ct)/2; % ����ƽ�����Դֲڶȣ�m
epsilon_e=1*10^(-3); % �����������ޣ����룩
epsilon_t=1*10^3; % ���������������룩
g=9.81; % �������ٶȣ�m/s^2��Ĭ�ϣ�

T_0=20; % ��Һ�����¶ȣ��棨���룩
P_0=0.1*10^6; % ��Һ����ѹ����Pa�����룩
rho_l_0=1150; % T_0��P_0�»�Һ�ܶȣ�kg/m^3�����룩
mu_l_0=0.03; % T_0��P_0�»�Һճ�ȣ�Pa*s�����룩
Qv_l_0=0.1/60; % ��Һ���������m^3/s�����룩
Qm_l_0=Qv_l_0*rho_l_0; % ��Һ����������kg/s
rho_g_0=0.655; % ע�����ܶȣ�kg/m^3
mu_g_0=RheologyG(T_0,P_0); % ����ճ�ȣ�Pa*s
Qv_g_0=4/60; % ע�������������m^3/s�����룩
Qm_g_0=Qv_g_0*rho_g_0;  % ע��������������kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % ��ĭ����������kg/s

D_s=1*10^(-3); % ɰ��ֱ����m�����룩
rho_s=2000; % ɰ���ܶȣ�kg/m^3�����룩
H_s=L_s_b-L_s_t; % �ײ�ɰ���߶ȣ�m
PHI=0.6; % ɰ����ӯ�ȣ����룩
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % ����ɰ����������kg

D_nozzle=4/1000; % ����ֱ����m�����룩
N_nozzle=3; % ������������룩
C=0.95; % ��������ϵ����ȡ0.95�����룩

C0=1.2; % Ư��������ϵ����Ĭ�ϣ�

M_s(1)=0; % ���׽�ɰ����kg/s
for t=2:1:Nt
    M_s(t)=0; % ���׽�ɰ����kg/s
end

OutPressure=1*10^6; % ����ѹ����Pa�����룩

%% �¶�����
T_i=20; % ��ĭע���¶ȣ���
T_g=0.02; % �����ݶȣ���/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % ����������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������
    end
end
T_a=T_ct; % ������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������

%% ��1��ʱ��ڵ㣨��ʼʱ�̣���ز�����ֵ���㣨���գ�
for x=1:1:Nx
    P_a(1,x)=P_a_3(x); % ����ѹ����Pa
    rho_g_a(1,x)=rho_g_a_3(x); % �����ܶȣ�kg/m^3
    rho_l_a(1,x)=rho_l_a_3(x); % ����Һ���ܶȣ�kg/m^3
    rho_f_a(1,x)=rho_f_a_3(x); % ������ĭ�ܶȣ�kg/m^3
    mu_g_a(1,x)=mu_g_a_3(x); % ����ճ�ȣ�Pa*s
    mu_l_a(1,x)=mu_l_a_3(x); % ����Һ��ճ�ȣ�Pa*s
    mu_f_a(1,x)=mu_f_a_3(x); % ������ĭճ�ȣ�Pa*s
    mu_s(1,x)=mu_s_3(x); % ����ճ�ȣ�Pa*s
    alpha_g_a(1,x)=alpha_g_a_3(x); % ���ຬ��
    alpha_l_a(1,x)=alpha_l_a_3(x); % Һ�ຬ��
    alpha_f_a(1,x)=alpha_f_a_3(x); % ��ĭ����
    alpha_s(1,x)=alpha_s_3(1,x); % ���ຬ��
    Va_g_a(1,x)=Va_g_a_3(x); % ���������٣�m/s
    Va_l_a(1,x)=Va_l_a_3(x); % ����Һ�������٣�m/s
    Va_f_a(1,x)=Va_f_a_3(x); % ������ĭ������٣�m/s
    Va_s(1,x)=Va_s_3(x); % ���������٣�m/s
    V_g_a(1,x)=V_g_a_3(x); % �������٣�m/s
    V_l_a(1,x)=V_l_a_3(x); % ����Һ�����٣�m/s
    V_f_a(1,x)=V_f_a_3(x); % ������ĭ���٣�m/s
    V_s(1,x)=V_s_3(x); % �������٣�m/s
    Vsr(1,x)=Vsr_3(x); % ���໬���ٶȣ�m/s
    gamma_g_a(1,x)=gamma_g_a_3(x); % ��ĭ����
    gamma_l_a(1,x)=gamma_l_a_3(1,x); % Һ��������
    V_m(1,x)=V_m_3(x); % ���ջ�������٣�m/s
    rho_m(1,x)=rho_m_3(x); % ���ջ�����ܶȣ�kg/m^3
    mu_m(1,x)=mu_m_3(x); % ���ջ����ճ�ȣ�Pa*s
    Ff_a(1,x)=Ff_a_3(x); % ���յ�λ����Ħ��ѹ����Pa/m
    f_a(1,x)=f_a_3(x); % ����Ħ������
    Re_a(1,x)=Re_a_3(x); % ������ŵ��
    flow_pattern_a(1,x)=flow_pattern_a_3(x); % ����������̬
end

%% ��2��Nt��ʱ��ڵ���ز������㣨���գ�
for t=2:1:Nt
    P_a(t,Nx_Lp+t-1)=P_a(t-1,Nx_Lp+t-1);  % ���չܵ�ѹ������ֵ��Pa
    
    err_OutPressure=1; % ����ѹ��������
    COUNT_OutPressure=0; % ����ѹ������������ֵ
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % ����ѹ����������
        
        % �ܵף���Nx_Lp+t-1���ռ�ڵ㣩����ز�������
        rho_g_a(t,Nx_Lp+t-1)=DensityG(T_a(t,Nx_Lp+t-1),P_a(t,Nx_Lp+t-1)); % ���������ܶȣ�kg/m^3
        rho_l_a(t,Nx_Lp+t-1)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Lp+t-1),P_a(t,Nx_Lp+t-1)); % ���ջ�Һ�ܶȣ�kg/m^3
        mu_g_a(t,Nx_Lp+t-1)=RheologyG(T_a(t,Nx_Lp+t-1),P_a(t,Nx_Lp+t-1)); % ��������ճ�ȣ�Pa*s
        mu_l_a(t,Nx_Lp+t-1)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Lp+t-1),P_a(t,Nx_Lp+t-1)); % ���ջ�Һճ�ȣ�Pa*s
        gamma_g_a(t,Nx_Lp+t-1)=(Qm_g_0/rho_g_a(t,Nx_Lp+t-1))/(Qm_g_0/rho_g_a(t,Nx_Lp+t-1)+Qm_l_0/rho_l_a(t,Nx_Lp+t-1)); % ��ĭ����
        gamma_l_a(t,Nx_Lp+t-1)=1-gamma_g_a(t,Nx_Lp+t-1); % Һ��������
        alpha_g_a(t,Nx_Lp+t-1)=alpha_f_a(t-1,Nx_Lp+t-1)*gamma_g_a(t,Nx_Lp+t-1); % �������庬��
        alpha_l_a(t,Nx_Lp+t-1)=alpha_f_a(t-1,Nx_Lp+t-1)*gamma_l_a(t,Nx_Lp+t-1); % ���ջ�Һ����
        rho_f_a(t,Nx_Lp+t-1)=rho_g_a(t,Nx_Lp+t-1)*gamma_g_a(t,Nx_Lp+t-1)+rho_l_a(t,Nx_Lp+t-1)*gamma_l_a(t,Nx_Lp+t-1); % ������ĭ�ܶȣ�kg/m^3
        mu_f_a(t,Nx_Lp+t-1)=mu_g_a(t,Nx_Lp+t-1)*gamma_g_a(t,Nx_Lp+t-1)+mu_l_a(t,Nx_Lp+t-1)*gamma_l_a(t,Nx_Lp+t-1); % ������ĭճ�ȣ�Pa*s
        mu_s(t,Nx_Lp+t-1)=mu_f_a(t,Nx_Lp+t-1); % ����ճ�ȣ�Pa*s
        Va_f_a(t,Nx_Lp+t-1)=Qm_f_0/(A_a(t,Nx_Lp+t-1)*rho_f_a(t,Nx_Lp+t-1)); % ������ĭ������٣�m/s
        Va_g_a(t,Nx_Lp+t-1)=Va_f_a(t,Nx_Lp+t-1); % �������������٣�m/s
        Va_l_a(t,Nx_Lp+t-1)=Va_f_a(t,Nx_Lp+t-1); % ���ջ�Һ������٣�m/s
        Va_s(t,Nx_Lp+t-1)=M_s(t)/(A_a(t,Nx_Lp+t-1)*rho_s); % ��м����ٶȣ�m/s
        Vsr(t,Nx_Lp+t-1)=12*(mu_f_a(t,Nx_Lp+t-1)/(rho_f_a(t,Nx_Lp+t-1)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Lp+t-1))/rho_f_a(t,Nx_Lp+t-1))*((rho_f_a(t,Nx_Lp+t-1)*D_s/mu_f_a(t,Nx_Lp+t-1))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
        alpha_s(t,Nx_Lp+t-1)=Va_s(t,Nx_Lp+t-1)/(C0*(Va_s(t,Nx_Lp+t-1)+Va_f_a(t,Nx_Lp+t-1))-Vsr(t,Nx_Lp+t-1));  % �����������
        V_s(t,Nx_Lp+t-1)=0;%Va_s(t,Nx_Lp+t-1)/alpha_s(t,Nx_Lp+t-1); % ��м�ٶȣ�m/s
        
        alpha_f_a(t,Nx_Lp+t-1)=1-alpha_s(t,Nx_Lp+t-1); % ������ĭ����
        alpha_g_a(t,Nx_Lp+t-1)=alpha_f_a(t,Nx_Lp+t-1)*gamma_g_a(t,Nx_Lp+t-1); % �������庬��
        alpha_l_a(t,Nx_Lp+t-1)=alpha_f_a(t,Nx_Lp+t-1)*gamma_l_a(t,Nx_Lp+t-1); % ���ջ�Һ����
        V_f_a(t,Nx_Lp+t-1)=Va_f_a(t,Nx_Lp+t-1)/alpha_f_a(t,Nx_Lp+t-1); % ������ĭ���٣�m/s
        V_g_a(t,Nx_Lp+t-1)=V_f_a(t,Nx_Lp+t-1); % �����������٣�m/s
        V_l_a(t,Nx_Lp+t-1)=V_f_a(t,Nx_Lp+t-1); % ���ջ�Һ���٣�m/s
        
        V_m(t,Nx_Lp+t-1)=Va_s(t,Nx_Lp+t-1)+Va_f_a(t,Nx_Lp+t-1); % ���ջ�����ٶȣ�m/s
        rho_m(t,Nx_Lp+t-1)=alpha_s(t,Nx_Lp+t-1)*rho_s+alpha_f_a(t,Nx_Lp+t-1)*rho_f_a(t,Nx_Lp+t-1); % ���ջ�����ܶȣ�kg/m^3
        mu_m(t,Nx_Lp+t-1)=alpha_s(t,Nx_Lp+t-1)*mu_s(t,Nx_Lp+t-1)+alpha_f_a(t,Nx_Lp+t-1)*mu_f_a(t,Nx_Lp+t-1); % ���ջ����ճ�ȣ�Pa*s
        [Ff_a(t,Nx_Lp+t-1),f_a(t,Nx_Lp+t-1),Re_a(t,Nx_Lp+t-1),flow_pattern_a(t,Nx_Lp+t-1)]=Friction_annulus(rho_m(t,Nx_Lp+t-1),V_m(t,Nx_Lp+t-1),mu_m(t,Nx_Lp+t-1),D_h(t,Nx_Lp+t-1),h_a,rho_f_a(t,Nx_Lp+t-1),V_f_a(t,Nx_Lp+t-1)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
        
        % ��Nx_Lp+t-2��1���ռ�ڵ㴦��ز�������
        for x=Nx_Lp+t-2:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % ����ѹ������ֵ��Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % ���ջ�Һ����
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            
            err_NodePressure=1; % ����ѹ��������
            COUNT_NodePressure=0; % ����ѹ������������ֵ
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % ����ѹ����������
                
                % ���߷��������������
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % ���������������ֵ1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % ���������������ֵ2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % ������������������
                COUNT_NodeEg=0; % ��������
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % �����������Ϊalpha_s_ass1ʱ������ֵ��
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % �����ٶȣ�m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % ���������٣�m/s
                    alpha_f_ass1=1-alpha_s_ass1; % ��ĭ�������
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % ��ĭ�ٶȣ�m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % ��ĭ������٣�m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % ���ջ�����ٶȣ�m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % ���������������ֵ
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % ����ĺ��������Ľ������ʵ�����������
                    
                    % �����������Ϊalpha_s_ass2ʱ������ֵ��
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % ���������غ㷽����ɢ��ʽ�м�ֵ����
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % �����ٶȣ�m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % ���������٣�m/s
                    alpha_f_ass2=1-alpha_s_ass2; % ��ĭ�������
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % ��ĭ�ٶȣ�m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % ��ĭ������٣�m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % ���ջ�����ٶȣ�m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % ���������������ֵ
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % ����ĺ��������Ľ������ʵ�����������
                    
                    % ���߷���������������
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % �µĹ��������������ֵalpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % ������������������
                    alpha_s_ass1=alpha_s_ass2; % �µĹ��������������ֵ1
                    alpha_s_ass2=alpha_s_ass3; % �µĹ��������������ֵ2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % ���������õ�����ʵ�����������ֵ����alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % �������������С��һ��ֵʱ����Ϊ���������Ϊ0�����ڷ�ֹ���ֺ����ļ������
                    alpha_s(t,x)=0; % �����������
                    V_s(t,x)=0; % �����ٶȣ�m/s
                    Va_s(t,x)=0; % �������ٶȣ�m/s
                    mu_s(t,x)=0; % ����ճ�ȣ�Pa*s
                    Vsr(t,x)=0; % ��м����ĩ�٣�m/s
                    
                    alpha_f_a(t,x)=1; % ��ĭ�������                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                else
                    V_s(t,x)=V_s_ass2; % �����ٶȣ�m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % ���������٣�m/s
                    Vsr(t,x)=Vsr_ass2; % ��м����ĩ�٣�m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % ��ĭ�ܶȣ�kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % ��ĭ�������
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    V_f_a(t,x)=V_f_ass2; % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % ���ջ�����ٶȣ�m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % ����ѹ������ֵ��Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % ����ѹ��������
                P_a(t,x)=P_new; % �µĻ���ѹ������ֵ��Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % ����ѹ��������        
        if (P_a(t,1)-OutPressure)>0 % ���ݳ���ѹ�������������Ի��չܵ�ѹ������ֵ���е���
            P_a(t,Nx_Lp+t-1)=P_a(t,Nx_Lp+t-1)-(P_a(t,1)-OutPressure)/2; % �µĻ��չܵ�ѹ������ֵ��Pa
        else
            P_a(t,Nx_Lp+t-1)=P_a(t,Nx_Lp+t-1)-(P_a(t,1)-OutPressure)/2*0.3; % �µĻ��չܵ�ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Lp+t��Nx_Dsp_pen(1)���ռ�ڵ㴦��ز�������
    for x=Nx_Lp+t:1:Nx_Dsp_pen(1)
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=0; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(1)+1��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)+1:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1-PHI; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=PHI; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
end

%% ��1��Ntʱ��ڵ���ͷѹ������
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Lp+t-1); % ������ĭ���������m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Lp+t-1)); % �����ɰ����ѹ����Pa����������������V_nozzle��m/s��
end

%% ��1��Ntʱ��ڵ���ز������㣨�������ڣ�
for t=1:1:Nt
    % �ܵף���Nx_Lp+t-1���ռ�ڵ㣩����ز�������
    P_ct(t,Nx_Lp+t-1)=P_a(t,Nx_Lp+t-1)+delta_P_SWT(t); % ����ѹ����Pa
    rho_g_ct(t,Nx_Lp+t-1)=DensityG(T_ct(t,Nx_Lp+t-1),P_ct(t,Nx_Lp+t-1)); % ���������ܶȣ�kg/m^3
    rho_l_ct(t,Nx_Lp+t-1)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Lp+t-1),P_ct(t,Nx_Lp+t-1)); % ���ڻ�Һ�ܶȣ�kg/m^3
    mu_g_ct(t,Nx_Lp+t-1)=RheologyG(T_ct(t,Nx_Lp+t-1),P_ct(t,Nx_Lp+t-1)); % ��������ճ�ȣ�Pa*s
    mu_l_ct(t,Nx_Lp+t-1)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Lp+t-1),P_ct(t,Nx_Lp+t-1)); % ���ڻ�Һճ�ȣ�Pa*s
    alpha_f_ct(t,Nx_Lp+t-1)=1; % ������ĭ����
    gamma_g_ct(t,Nx_Lp+t-1)=(Qm_g_0/rho_g_ct(t,Nx_Lp+t-1))/(Qm_g_0/rho_g_ct(t,Nx_Lp+t-1)+Qm_l_0/rho_l_ct(t,Nx_Lp+t-1)); % ������ĭ����
    gamma_l_ct(t,Nx_Lp+t-1)=1-gamma_g_ct(t,Nx_Lp+t-1); % ����Һ��������
    alpha_g_ct(t,Nx_Lp+t-1)=alpha_f_ct(t,Nx_Lp+t-1)*gamma_g_ct(t,Nx_Lp+t-1); % �������庬��
    alpha_l_ct(t,Nx_Lp+t-1)=alpha_f_ct(t,Nx_Lp+t-1)*gamma_l_ct(t,Nx_Lp+t-1); % ���ڻ�Һ����
    rho_f_ct(t,Nx_Lp+t-1)=rho_g_ct(t,Nx_Lp+t-1)*gamma_g_ct(t,Nx_Lp+t-1)+rho_l_ct(t,Nx_Lp+t-1)*gamma_l_ct(t,Nx_Lp+t-1); % ������ĭ�ܶȣ�kg/m^3
    mu_f_ct(t,Nx_Lp+t-1)=mu_g_ct(t,Nx_Lp+t-1)*gamma_g_ct(t,Nx_Lp+t-1)+mu_l_ct(t,Nx_Lp+t-1)*gamma_l_ct(t,Nx_Lp+t-1); % ������ĭճ�ȣ�Pa*s
    V_f_ct(t,Nx_Lp+t-1)=Qm_f_0/(rho_f_ct(t,Nx_Lp+t-1)*A_ct(t,Nx_Lp+t-1)); % �����������٣�m/s
    [Ff_ct(t,Nx_Lp+t-1),f_ct(t,Nx_Lp+t-1),Re_ct(t,Nx_Lp+t-1),flow_pattern_ct(t,Nx_Lp+t-1)]=Friction_coiledtubing(rho_f_ct(t,Nx_Lp+t-1),V_f_ct(t,Nx_Lp+t-1),mu_f_ct(t,Nx_Lp+t-1),D_ct_i(t,Nx_Lp+t-1),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    % ��Nx_Lp+t-2��1���ռ�ڵ㴦��ز�������
    for x=Nx_Lp+t-2:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % ����ѹ������ֵ��Pa
        
        err_DriPipePressure=1; % ����ѹ��������
        COUNT_DriPipePressure=0; % ����ѹ������������ֵ
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % ����ѹ����������
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һ�ܶȣ�kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һճ�ȣ�Pa*s
            alpha_f_ct(t,x)=1; % ������ĭ����
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % ������ĭ����
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % �������庬��
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % ���ڻ�Һ����
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭճ�ȣ�Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % ������ĭ���٣�m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % ����ѹ����Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % �������ѹ������ֵ�����ֵ֮���������
            P_ct_ass(t,x)=P_ct(t,x); % �µĹ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Lp+t��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Lp+t:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % �������庬��
        alpha_l_ct(t,x)=alpha_l_a(t,x); % ���ڻ�Һ����
        gamma_g_ct(t,x)=gamma_g_a(t,x); % ������ĭ����
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
        rho_f_ct(t,x)=rho_m(t,x); % ���������ܶȣ�kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % ��������ճ�ȣ�Pa*s
        V_f_ct(t,t)=V_m(t,x); % �����������٣�sm/s
        Ff_ct(t,x)=Ff_a(t,x); % �������嵥λ����Ħ��ѹ����Pa/m��
        f_ct(t,x)=f_a(t,x); % �������巶��Ħ������
        Re_ct(t,x)=Re_a(t,x); % ����������ŵ����������̬
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % ����������̬
        P_ct(t,x)=P_a(t,x); % ����ѹ����Pa
    end
end

%% ��1��Ntʱ��ڵ��̹ܶγ���ѹ������ѹ���㣨�������ڣ�
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % ��������������٣�m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % ����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % �̹ܶ����嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % �����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % ���ڵ�ע��ͷ������Ħ��ѹ����Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % ��������Ħ��ѹ����Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % ����������Ͳ��Ħ��ѹ����Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % �̹ܶ�Ħ��ѹ����Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % �̹ܶγ���ѹ����Pa
    P_pumb(t)=P_coil(t)+delta_P_coil(t); % ��ѹ��Pa
end

%% ���㾮���ۻ���ɰ��
M_w_tem(1)=0; % ����˲ʱ��ɰ����ֵ��kg
M_w_tot(1)=0; % �����ۻ���ɰ����ֵ��kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮��˲ʱ��ɰ����kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮���ۻ���ɰ����kg
end

%% ����ĩ̬��Ͳ���ɰŨ��
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(2)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% ���㻷����ĭ����ƽ��ֵV_l_mre
V_f_mre=0; % ���շ���ƽ��ֵ��m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % ���շ���ƽ��ֵ��m/s
    end
end

%% ����������ĩ��ƽ��ֵVsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % ��Vsr��������ȡֵ����ĳ���ĩ�٣�m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % �������ĩ��ƽ��ֵ��m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % �������ĩ��ƽ��ֵ��m/s
end

%% �ж��Ƿ���Ч��ɰ�����շ���ƽ��ֵ����2��ɰ������ĩ��ƽ��ֵ�����Ƿ���ɳ�ɰ
fprintf("�������������1��\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % �����д���˵����ǰ�����¡���ɳ�ɰ��
    else
        fprintf("Sand Cleanout UnFinished!\n"); % �����д���˵����ǰ�����¡�δ��ɳ�ɰ��
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
end

%% ��ͲECD����
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % ����ECD��kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % ����ECD��kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
end

%% �������������ĩ״̬����
alpha_f_a_4=alpha_f_a(Nt,:); % �������������ĩ״̬��ĭ����
alpha_g_a_4=alpha_g_a(Nt,:); % �������������ĩ״̬���ຬ��
alpha_l_a_4=alpha_l_a(Nt,:); % �������������ĩ״̬Һ�ຬ��
alpha_s_4=alpha_s(Nt,:); % �������������ĩ״̬���ຬ��
f_a_4=f_a(Nt,:); % �������������ĩ״̬����Ħ������
Ff_a_4=Ff_a(Nt,:); % �������������ĩ״̬���յ�λ����Ħ��ѹ����Pa/m
flow_pattern_a_4=flow_pattern_a(Nt,:); % �������������ĩ״̬����������̬
gamma_g_a_4=gamma_g_a(Nt,:); % �������������ĩ״̬��ĭ����
gamma_l_a_4=gamma_l_a(Nt,:); % �������������ĩ״̬Һ��������
mu_f_a_4=mu_f_a(Nt,:); % �������������ĩ״̬������ĭճ�ȣ�Pa*s
mu_g_a_4=mu_g_a(Nt,:); % �������������ĩ״̬����ճ�ȣ�Pa*s
mu_l_a_4=mu_l_a(Nt,:); % �������������ĩ״̬����Һ��ճ�ȣ�Pa*s
mu_m_4=mu_m(Nt,:); % �������������ĩ״̬���ջ����ճ�ȣ�Pa*s
mu_s_4=mu_s(Nt,:); % �������������ĩ״̬����ճ�ȣ�Pa*s
P_a_4=P_a(Nt,:); % �������������ĩ״̬����ѹ����Pa
Re_a_4=Re_a(Nt,:); % �������������ĩ״̬������ŵ��
rho_f_a_4=rho_f_a(Nt,:); % �������������ĩ״̬������ĭ�ܶȣ�kg/m^3
rho_g_a_4=rho_g_a(Nt,:); % �������������ĩ״̬���������ܶȣ�kg/m^3
rho_l_a_4=rho_l_a(Nt,:); % �������������ĩ״̬����Һ���ܶȣ�kg/m^3
rho_m_4=rho_m(Nt,:); % �������������ĩ״̬���ջ�����ܶȣ�kg/m^3
V_f_a_4=V_f_a(Nt,:); % �������������ĩ״̬������ĭ���٣�m/s
V_g_a_4=V_g_a(Nt,:); % �������������ĩ״̬�������٣�m/s
V_l_a_4=V_l_a(Nt,:); % �������������ĩ״̬����Һ�����٣�m/s
V_m_4=V_m(Nt,:); % �������������ĩ״̬���ջ�������٣�m/s
V_s_4=V_s(Nt,:); % �������������ĩ״̬�������٣�m/s
Va_f_a_4=Va_f_a(Nt,:); % �������������ĩ״̬������ĭ������٣�m/s
Va_g_a_4=Va_g_a(Nt,:); % �������������ĩ״̬���������٣�m/s
Va_l_a_4=Va_l_a(Nt,:); % �������������ĩ״̬����Һ�������٣�m/s
Va_s_4=Va_s(Nt,:); % �������������ĩ״̬���������٣�m/s
Vsr_4=Vsr(Nt,:); % �������������ĩ״̬���໬���ٶȣ�m/s

%% ���ݴ洢�����������룩
ANS_Nt_4=Nt; % ʱ��ڵ���
ANS_Time_4=Time(1:Nt); % ����ѭ��ʱ�䣬s
ANS_alpha_g_a_4=alpha_g_a(1:Nt,:); % ���������������
ANS_alpha_l_a_4=alpha_l_a(1:Nt,:); % ����Һ���������
ANS_alpha_s_4=alpha_s(1:Nt,:); % ���ຬ��
ANS_P_ct_4=P_ct(1:Nt,:); % ����ѹ����Pa��
ANS_delta_P_SWT_4=delta_P_SWT(1:Nt); % ��ɰ����ѹ����Pa��
ANS_P_a_4=P_a(1:Nt,:); % ����ѹ����Pa��
ANS_T_a_4=T_a(1:Nt,:); % �����¶ȣ��棩
ANS_T_ct_4=T_ct(1:Nt,:); % �����¶ȣ��棩
ANS_P_coil_4=P_coil(1:Nt); % �̹ܶγ���ѹ����Pa��
ANS_P_pump_4=P_pump(1:Nt); % ��ѹ��Pa��
ANS_M_w_tem_4=M_w_tem(1:Nt); % ����˲ʱ��ɰ����kg��
ANS_M_w_tot_4=M_w_tot(1:Nt); % �����ۻ���ɰ����kg��
ANS_Va_s_4=Va_s(1:Nt,:); % ��м����ٶȣ�m/s��
ANS_Va_f_a_4=Va_f_a(1:Nt,:); % ������ĭ������٣�m/s��
ANS_V_s_4=V_s(1:Nt,:); % ��м�����ٶȣ�m/s��
ANS_V_f_a_4=V_f_a(1:Nt,:); % ������ĭ���٣�m/s��
ANS_alpha_g_ct_4=alpha_g_ct(1:Nt,:); % ���������������
ANS_alpha_l_ct_4=alpha_l_ct(1:Nt,:); % ����Һ���������
ANS_gamma_g_ct_4=gamma_g_ct(1:Nt,:); % ������ĭ����
ANS_gamma_g_a_4=gamma_g_a(1:Nt,:); % ������ĭ����
ANS_ECD_a_4=ECD_a(1:Nt,:); % ����ECD��kg/m^3��
ANS_L_coil_4=L_coil(1:Nt); % ���������m��



%% ��ϴ�������2�����ʵ�λ�ƣ�
%% ������ֵ����
dt=zeros(); % ʱ�䲽����s
dx=zeros(); % �ռ䲽����m
Time=zeros(); % ������������ʱ����s
L_coil=zeros(); % ���������m
L_reel=zeros(); % �̹ܶγ��ȣ�m
D_t_i=zeros(); % �͹��ھ���m
D_ct_o=zeros(); % �����͹��⾶��m
D_ct_i=zeros(); % �����͹��ھ���m
M_s=zeros(); % ���׽�ɰ����kg/s
P_a=zeros(); % ����ѹ����Pa
rho_g_a=zeros(); % �����ܶȣ�kg/m^3
rho_l_a=zeros(); % ����Һ���ܶȣ�kg/m^3
rho_f_a=zeros(); % ������ĭ�ܶȣ�kg/m^3
mu_g_a=zeros(); % ����ճ�ȣ�Pa*s
mu_l_a=zeros(); % ����Һ��ճ�ȣ�Pa*s
mu_f_a=zeros(); % ������ĭճ�ȣ�Pa*s
mu_s=zeros(); % ����ճ�ȣ�Pa*s
alpha_g_a=zeros(); % ���ຬ��
alpha_l_a=zeros(); % Һ�ຬ��
alpha_f_a=zeros(); % ��ĭ����
alpha_s=zeros(); % ���ຬ��
Va_g_a=zeros(); % ���������٣�m/s
Va_l_a=zeros(); % ����Һ�������٣�m/s
Va_f_a=zeros(); % ������ĭ������٣�m/s
Va_s=zeros(); % ���������٣�m/s
V_g_a=zeros(); % �������٣�m/s
V_l_a=zeros(); % ����Һ�����٣�m/s
V_f_a=zeros(); % ������ĭ���٣�m/s
V_s=zeros(); % �������٣�m/s
Vsr=zeros(); % ���໬���ٶȣ�m/s
gamma_g_a=zeros(); % ��ĭ����
gamma_l_a=zeros(); % Һ��������
V_m=zeros(); % ���ջ�������٣�m/s
rho_m=zeros(); % ���ջ�����ܶȣ�kg/m^3
mu_m=zeros(); % ���ջ����ճ�ȣ�Pa*s
Ff_a=zeros(); % ���յ�λ����Ħ��ѹ����Pa/m
f_a=zeros(); % ����Ħ������
Re_a=zeros(); % ������ŵ��
flow_pattern_a=zeros(); % ����������̬

%% ����ռ䲽��dx����Ӧʱ�䲽��dt
V5=(Dsp_pen(2)-Dsp_pen(1))/(1*3600); % ��ϴ����ٶȣ�m/s�����룩
t_5=(Dsp_pen(2)-Dsp_pen(1))/V5; % ��ϴ�����ʱ����s

Nt=Nx_Dsp_pen(2)-Nx_Dsp_pen(1)+1; % ʱ��ڵ���
nt=Nt-1; % ʱ��������

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % �ռ䲽����m
end

Time(1)=0; % ��ϴ�����ʱ����ֵ��s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Dsp_pen(1)+t-1)/V5; % ÿ��ϴ���һ���ռ䲽������ʱ�䣬s
    Time(t+1)=Time(t)+dt(t); % ��ϴ�������(t+1)���ռ�ڵ�����������ʱ����s
end

%% ���㲻ͬʱ�������ܳ�ϴ������L_coil���̹ܶ������ܳ���L_reel
L=10000; % �����͹��ܳ���m�����룩
L_wg=8; % ���ڵ�ע��ͷ�����������ܳ��ȣ�m�����룩
L_goose=3; % �������������ܳ��ȣ�m�����룩
D_goose=2; % �������ΰ뾶��m�����룩
H_goose=10; % �����������߶ȣ�m�����룩
L_gr=20; % ����������Ͳ�������ܳ��ȣ�m�����룩
theta_gr=acosd(H_goose/L_gr); % ����������Ͳ����������Ǧ���߼нǣ���
D_r_i=3; % ��Ͳ�ھ���m�����룩
D_r_o=5; % ��Ͳ�⾶��m�����룩
W_r=5; % ��Ͳ��ȣ�m�����룩
D_cable=0.005; % �����⾶��m

L_coil(1)=Dsp_pen(1); % ��ʼʱ�������ܵײ���ȣ�m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % ��ʼʱ���̹ܶ������ܳ��ȣ�m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)+dx(Nx_Dsp_pen(1)+t-2); % ��ϴ�����ȣ�m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % �̹ܶ������ܳ��ȣ�m
end
L_cable=L_coil; % �������볤�ȣ�m

%% �������뼰Ԥ����
D_ct_o_0=0.04445; % �����͹��⾶��m�����룩
L1=2000; % �����ڵ�һ�������ܳ��ȣ�m�����룩
D_ct_i_1=0.03709; % �����ڵ�һ���������ھ���m�����룩
L2=2000; % �����ڵڶ��������ܳ��ȣ�m�����룩
D_ct_i_2=0.03653; % �����ڵڶ����������ھ���m�����룩
L3=2000; % �����ڵ����������ܳ��ȣ�m�����룩
D_ct_i_3=0.03555; % �����ڵ������������ھ���m�����룩
L4=2000; % �����ڵ��Ķ������ܳ��ȣ�m�����룩
D_ct_i_4=0.03489; % �����ڵ��Ķ��������ھ���m�����룩
L5=2000; % �����ڵ���������ܳ��ȣ�m�����룩
D_ct_i_5=0.03409; % �����ڵ�����������ھ���m�����룩

L_t_1=4000; % �ϲ��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_1=0.09718; %0.068;% �ϲ��͹ܣ����׹ܻ����ۣ��ھ���m�����룩
L_t_2=2200; % �²��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_2=0.09718; %0.13970; % �²��͹ܣ����׹ܻ����ۣ��ھ���m�����룩

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % �͹ܣ����׹ܻ����ۣ��ھ���m
        else
            D_t_i(t,x)=D_t_i_2; % �͹ܣ����׹ܻ����ۣ��ھ���m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(2)
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % �����͹��⾶��m
            else
                D_ct_o(t,x)=0; % �����͹��⾶��m
            end
        else
            D_ct_o(t,x)=0; % �����͹��⾶��m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(2)
            if Depth(x)<=L_coil(t) % �����ļ�����������ھ�
                if (L_coil(t)-Depth(x))<=L1 %&& (L_coil(t)-Depth(x))>=0
                    D_ct_i(t,x)=D_ct_i_1; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % �����͹��ھ���m
                end
            else
                D_ct_i(t,x)=0; % �����͹��ھ���m
            end
        else
            D_ct_i(t,x)=0; % �����͹��ھ���m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % �����ļ������������ֵ
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % �����͹��ڽ������m^2
        else
            A_ct(t,x)=0; % �����͹��ڽ������m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % ����ˮ��ֱ����m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % ���ս������m^2
    end
end

%% �������뼰Ԥ����
h_t=25.4*10^(-6); % �͹ܣ����׹ܻ����ۣ����Դֲڶȣ�m�����룩
h_ct=25.4*10^(-6); % �����͹ܾ��Դֲڶȣ�m�����룩
h_a=(h_t+h_ct)/2; % ����ƽ�����Դֲڶȣ�m
epsilon_e=1*10^(-3); % �����������ޣ����룩
epsilon_t=1*10^3; % ���������������룩
g=9.81; % �������ٶȣ�m/s^2��Ĭ�ϣ�

T_0=20; % ��Һ�����¶ȣ��棨���룩
P_0=0.1*10^6; % ��Һ����ѹ����Pa�����룩
rho_l_0=1150; % T_0��P_0�»�Һ�ܶȣ�kg/m^3�����룩
mu_l_0=0.03; % T_0��P_0�»�Һճ�ȣ�Pa*s�����룩
Qv_l_0=0.1/60; % ��Һ���������m^3/s�����룩
Qm_l_0=Qv_l_0*rho_l_0; % ��Һ����������kg/s
rho_g_0=0.655; % ע�����ܶȣ�kg/m^3
mu_g_0=RheologyG(T_0,P_0); % ����ճ�ȣ�Pa*s
Qv_g_0=4/60; % ע�������������m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % ע��������������kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % ��ĭ����������kg/s

D_s=1*10^(-3); % ɰ��ֱ����m�����룩
rho_s=2000; % ɰ���ܶȣ�kg/m^3�����룩
H_s=L_s_b-L_s_t; % �ײ�ɰ���߶ȣ�m
PHI=0.6; % ɰ����ӯ�ȣ����룩
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % ����ɰ����������kg

D_nozzle=4/1000; % ����ֱ����m�����룩
N_nozzle=3; % ������������룩
C=0.95; % ��������ϵ����ȡ0.95�����룩

C0=1.2; % Ư��������ϵ����Ĭ�ϣ�

M_s(1)=0; % ���׽�ɰ����kg/s
for t=2:1:Nt
    M_s(t)=PHI*rho_s*V5*0.25*pi*D_t_i(t,Nx_Dsp_pen(1)+t-1)^2; % ���׽�ɰ����kg/s
end

OutPressure=1*10^6; % ����ѹ����Pa�����룩

%% �¶�����
T_i=20; % ��ĭע���¶ȣ���
T_g=0.02; % �����ݶȣ���/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % ����������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������
    end
end
T_a=T_ct; % ������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������

%% ��1��ʱ��ڵ㣨��ʼʱ�̣���ز�����ֵ���㣨���գ�
for x=1:1:Nx
    P_a(1,x)=P_a_4(x); % ����ѹ����Pa
    rho_g_a(1,x)=rho_g_a_4(x); % �����ܶȣ�kg/m^3
    rho_l_a(1,x)=rho_l_a_4(x); % ����Һ���ܶȣ�kg/m^3
    rho_f_a(1,x)=rho_f_a_4(x); % ������ĭ�ܶȣ�kg/m^3
    mu_g_a(1,x)=mu_g_a_4(x); % ����ճ�ȣ�Pa*s
    mu_l_a(1,x)=mu_l_a_4(x); % ����Һ��ճ�ȣ�Pa*s
    mu_f_a(1,x)=mu_f_a_4(x); % ������ĭճ�ȣ�Pa*s
    mu_s(1,x)=mu_s_4(x); % ����ճ�ȣ�Pa*s
    alpha_g_a(1,x)=alpha_g_a_4(x); % ���ຬ��
    alpha_l_a(1,x)=alpha_l_a_4(x); % Һ�ຬ��
    alpha_f_a(1,x)=alpha_f_a_4(x); % ��ĭ����
    alpha_s(1,x)=alpha_s_4(1,x); % ���ຬ��
    Va_g_a(1,x)=Va_g_a_4(x); % ���������٣�m/s
    Va_l_a(1,x)=Va_l_a_4(x); % ����Һ�������٣�m/s
    Va_f_a(1,x)=Va_f_a_4(x); % ������ĭ������٣�m/s
    Va_s(1,x)=Va_s_4(x); % ���������٣�m/s
    V_g_a(1,x)=V_g_a_4(x); % �������٣�m/s
    V_l_a(1,x)=V_l_a_4(x); % ����Һ�����٣�m/s
    V_f_a(1,x)=V_f_a_4(x); % ������ĭ���٣�m/s
    V_s(1,x)=V_s_4(x); % �������٣�m/s
    Vsr(1,x)=Vsr_4(x); % ���໬���ٶȣ�m/s
    gamma_g_a(1,x)=gamma_g_a_4(x); % ��ĭ����
    gamma_l_a(1,x)=gamma_l_a_4(1,x); % Һ��������
    V_m(1,x)=V_m_4(x); % ���ջ�������٣�m/s
    rho_m(1,x)=rho_m_4(x); % ���ջ�����ܶȣ�kg/m^3
    mu_m(1,x)=mu_m_4(x); % ���ջ����ճ�ȣ�Pa*s
    Ff_a(1,x)=Ff_a_4(x); % ���յ�λ����Ħ��ѹ����Pa/m
    f_a(1,x)=f_a_4(x); % ����Ħ������
    Re_a(1,x)=Re_a_4(x); % ������ŵ��
    flow_pattern_a(1,x)=flow_pattern_a_4(x); % ����������̬
end

%% ��2��Nt��ʱ��ڵ���ز������㣨���գ�
for t=2:1:Nt
    P_a(t,Nx_Dsp_pen(1)+t-1)=P_a(t-1,Nx_Dsp_pen(1)+t-1);  % ���չܵ�ѹ������ֵ��Pa
    
    err_OutPressure=1; % ����ѹ��������
    COUNT_OutPressure=0; % ����ѹ������������ֵ
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % ����ѹ����������
        
        % �ܵף���Nx_Dsp_pen(1)+t-1���ռ�ڵ㣩����ز�������
        rho_g_a(t,Nx_Dsp_pen(1)+t-1)=DensityG(T_a(t,Nx_Dsp_pen(1)+t-1),P_a(t,Nx_Dsp_pen(1)+t-1)); % ���������ܶȣ�kg/m^3
        rho_l_a(t,Nx_Dsp_pen(1)+t-1)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)+t-1),P_a(t,Nx_Dsp_pen(1)+t-1)); % ���ջ�Һ�ܶȣ�kg/m^3
        mu_g_a(t,Nx_Dsp_pen(1)+t-1)=RheologyG(T_a(t,Nx_Dsp_pen(1)+t-1),P_a(t,Nx_Dsp_pen(1)+t-1)); % ��������ճ�ȣ�Pa*s
        mu_l_a(t,Nx_Dsp_pen(1)+t-1)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(1)+t-1),P_a(t,Nx_Dsp_pen(1)+t-1)); % ���ջ�Һճ�ȣ�Pa*s
        gamma_g_a(t,Nx_Dsp_pen(1)+t-1)=(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)+t-1))/(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(1)+t-1)+Qm_l_0/rho_l_a(t,Nx_Dsp_pen(1)+t-1)); % ��ĭ����
        gamma_l_a(t,Nx_Dsp_pen(1)+t-1)=1-gamma_g_a(t,Nx_Dsp_pen(1)+t-1); % Һ��������
        alpha_g_a(t,Nx_Dsp_pen(1)+t-1)=alpha_f_a(t-1,Nx_Dsp_pen(1)+t-1)*gamma_g_a(t,Nx_Dsp_pen(1)+t-1); % �������庬��
        alpha_l_a(t,Nx_Dsp_pen(1)+t-1)=alpha_f_a(t-1,Nx_Dsp_pen(1)+t-1)*gamma_l_a(t,Nx_Dsp_pen(1)+t-1); % ���ջ�Һ����
        rho_f_a(t,Nx_Dsp_pen(1)+t-1)=rho_g_a(t,Nx_Dsp_pen(1)+t-1)*gamma_g_a(t,Nx_Dsp_pen(1)+t-1)+rho_l_a(t,Nx_Dsp_pen(1)+t-1)*gamma_l_a(t,Nx_Dsp_pen(1)+t-1); % ������ĭ�ܶȣ�kg/m^3
        mu_f_a(t,Nx_Dsp_pen(1)+t-1)=mu_g_a(t,Nx_Dsp_pen(1)+t-1)*gamma_g_a(t,Nx_Dsp_pen(1)+t-1)+mu_l_a(t,Nx_Dsp_pen(1)+t-1)*gamma_l_a(t,Nx_Dsp_pen(1)+t-1); % ������ĭճ�ȣ�Pa*s
        mu_s(t,Nx_Dsp_pen(1)+t-1)=mu_f_a(t,Nx_Dsp_pen(1)+t-1); % ����ճ�ȣ�Pa*s
        Va_f_a(t,Nx_Dsp_pen(1)+t-1)=Qm_f_0/(A_a(t,Nx_Dsp_pen(1)+t-1)*rho_f_a(t,Nx_Dsp_pen(1)+t-1)); % ������ĭ������٣�m/s
        Va_g_a(t,Nx_Dsp_pen(1)+t-1)=Va_f_a(t,Nx_Dsp_pen(1)+t-1); % �������������٣�m/s
        Va_l_a(t,Nx_Dsp_pen(1)+t-1)=Va_f_a(t,Nx_Dsp_pen(1)+t-1); % ���ջ�Һ������٣�m/s
        Va_s(t,Nx_Dsp_pen(1)+t-1)=M_s(t)/(A_a(t,Nx_Dsp_pen(1)+t-1)*rho_s); % ��м����ٶȣ�m/s
        Vsr(t,Nx_Dsp_pen(1)+t-1)=12*(mu_f_a(t,Nx_Dsp_pen(1)+t-1)/(rho_f_a(t,Nx_Dsp_pen(1)+t-1)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Dsp_pen(1)+t-1))/rho_f_a(t,Nx_Dsp_pen(1)+t-1))*((rho_f_a(t,Nx_Dsp_pen(1)+t-1)*D_s/mu_f_a(t,Nx_Dsp_pen(1)+t-1))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
        alpha_s(t,Nx_Dsp_pen(1)+t-1)=Va_s(t,Nx_Dsp_pen(1)+t-1)/(C0*(Va_s(t,Nx_Dsp_pen(1)+t-1)+Va_f_a(t,Nx_Dsp_pen(1)+t-1))-Vsr(t,Nx_Dsp_pen(1)+t-1));  % �����������
        V_s(t,Nx_Dsp_pen(1)+t-1)=Va_s(t,Nx_Dsp_pen(1)+t-1)/alpha_s(t,Nx_Dsp_pen(1)+t-1); % ��м�ٶȣ�m/s
        
        alpha_f_a(t,Nx_Dsp_pen(1)+t-1)=1-alpha_s(t,Nx_Dsp_pen(1)+t-1); % ������ĭ����
        alpha_g_a(t,Nx_Dsp_pen(1)+t-1)=alpha_f_a(t,Nx_Dsp_pen(1)+t-1)*gamma_g_a(t,Nx_Dsp_pen(1)+t-1); % �������庬��
        alpha_l_a(t,Nx_Dsp_pen(1)+t-1)=alpha_f_a(t,Nx_Dsp_pen(1)+t-1)*gamma_l_a(t,Nx_Dsp_pen(1)+t-1); % ���ջ�Һ����
        V_f_a(t,Nx_Dsp_pen(1)+t-1)=Va_f_a(t,Nx_Dsp_pen(1)+t-1)/alpha_f_a(t,Nx_Dsp_pen(1)+t-1); % ������ĭ���٣�m/s
        V_g_a(t,Nx_Dsp_pen(1)+t-1)=V_f_a(t,Nx_Dsp_pen(1)+t-1); % �����������٣�m/s
        V_l_a(t,Nx_Dsp_pen(1)+t-1)=V_f_a(t,Nx_Dsp_pen(1)+t-1); % ���ջ�Һ���٣�m/s
        
        V_m(t,Nx_Dsp_pen(1)+t-1)=Va_s(t,Nx_Dsp_pen(1)+t-1)+Va_f_a(t,Nx_Dsp_pen(1)+t-1); % ���ջ�����ٶȣ�m/s
        rho_m(t,Nx_Dsp_pen(1)+t-1)=alpha_s(t,Nx_Dsp_pen(1)+t-1)*rho_s+alpha_f_a(t,Nx_Dsp_pen(1)+t-1)*rho_f_a(t,Nx_Dsp_pen(1)+t-1); % ���ջ�����ܶȣ�kg/m^3
        mu_m(t,Nx_Dsp_pen(1)+t-1)=alpha_s(t,Nx_Dsp_pen(1)+t-1)*mu_s(t,Nx_Dsp_pen(1)+t-1)+alpha_f_a(t,Nx_Dsp_pen(1)+t-1)*mu_f_a(t,Nx_Dsp_pen(1)+t-1); % ���ջ����ճ�ȣ�Pa*s
        [Ff_a(t,Nx_Dsp_pen(1)+t-1),f_a(t,Nx_Dsp_pen(1)+t-1),Re_a(t,Nx_Dsp_pen(1)+t-1),flow_pattern_a(t,Nx_Dsp_pen(1)+t-1)]=Friction_annulus(rho_m(t,Nx_Dsp_pen(1)+t-1),V_m(t,Nx_Dsp_pen(1)+t-1),mu_m(t,Nx_Dsp_pen(1)+t-1),D_h(t,Nx_Dsp_pen(1)+t-1),h_a,rho_f_a(t,Nx_Dsp_pen(1)+t-1),V_f_a(t,Nx_Dsp_pen(1)+t-1)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
        
        % ��Nx_Dsp_pen(1)+t-2��1���ռ�ڵ㴦��ز�������
        for x=Nx_Dsp_pen(1)+t-2:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % ����ѹ������ֵ��Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % ���ջ�Һ����
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            
            err_NodePressure=1; % ����ѹ��������
            COUNT_NodePressure=0; % ����ѹ������������ֵ
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % ����ѹ����������
                
                % ���߷��������������
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % ���������������ֵ1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % ���������������ֵ2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % ������������������
                COUNT_NodeEg=0; % ��������
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % �����������Ϊalpha_s_ass1ʱ������ֵ��
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % �����ٶȣ�m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % ���������٣�m/s
                    alpha_f_ass1=1-alpha_s_ass1; % ��ĭ�������
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % ��ĭ�ٶȣ�m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % ��ĭ������٣�m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % ���ջ�����ٶȣ�m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % ���������������ֵ
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % ����ĺ��������Ľ������ʵ�����������
                    
                    % �����������Ϊalpha_s_ass2ʱ������ֵ��
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % ���������غ㷽����ɢ��ʽ�м�ֵ����
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % �����ٶȣ�m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % ���������٣�m/s
                    alpha_f_ass2=1-alpha_s_ass2; % ��ĭ�������
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % ��ĭ�ٶȣ�m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % ��ĭ������٣�m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % ���ջ�����ٶȣ�m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % ���������������ֵ
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % ����ĺ��������Ľ������ʵ�����������
                    
                    % ���߷���������������
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % �µĹ��������������ֵalpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % ������������������
                    alpha_s_ass1=alpha_s_ass2; % �µĹ��������������ֵ1
                    alpha_s_ass2=alpha_s_ass3; % �µĹ��������������ֵ2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % ���������õ�����ʵ�����������ֵ����alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % �������������С��һ��ֵʱ����Ϊ���������Ϊ0�����ڷ�ֹ���ֺ����ļ������
                    alpha_s(t,x)=0; % �����������
                    V_s(t,x)=0; % �����ٶȣ�m/s
                    Va_s(t,x)=0; % �������ٶȣ�m/s
                    mu_s(t,x)=0; % ����ճ�ȣ�Pa*s
                    Vsr(t,x)=0; % ��м����ĩ�٣�m/s
                    
                    alpha_f_a(t,x)=1; % ��ĭ�������                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                else
                    V_s(t,x)=V_s_ass2; % �����ٶȣ�m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % ���������٣�m/s
                    Vsr(t,x)=Vsr_ass2; % ��м����ĩ�٣�m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % ��ĭ�ܶȣ�kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % ��ĭ�������
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    V_f_a(t,x)=V_f_ass2; % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % ���ջ�����ٶȣ�m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % ����ѹ������ֵ��Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % ����ѹ��������
                P_a(t,x)=P_new; % �µĻ���ѹ������ֵ��Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % ����ѹ��������        
        if (P_a(t,1)-OutPressure)>0 % ���ݳ���ѹ�������������Ի��չܵ�ѹ������ֵ���е���
            P_a(t,Nx_Dsp_pen(1)+t-1)=P_a(t,Nx_Dsp_pen(1)+t-1)-(P_a(t,1)-OutPressure)/2; % �µĻ��չܵ�ѹ������ֵ��Pa
        else
            P_a(t,Nx_Dsp_pen(1)+t-1)=P_a(t,Nx_Dsp_pen(1)+t-1)-(P_a(t,1)-OutPressure)/2*0.3; % �µĻ��չܵ�ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(1)+t��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)+t:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1-PHI; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=PHI; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
end

%% ��1��Ntʱ��ڵ���ͷѹ������
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Dsp_pen(1)+t-1); % ������ĭ���������m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Dsp_pen(1)+t-1)); % �����ɰ����ѹ����Pa����������������V_nozzle��m/s��
end

%% ��1��Ntʱ��ڵ���ز������㣨�������ڣ�
for t=1:1:Nt
    % �ܵף���Nx_Dsp_pen(1)+t-1���ռ�ڵ㣩����ز�������
    P_ct(t,Nx_Dsp_pen(1)+t-1)=P_a(t,Nx_Dsp_pen(1)+t-1)+delta_P_SWT(t); % ����ѹ����Pa
    rho_g_ct(t,Nx_Dsp_pen(1)+t-1)=DensityG(T_ct(t,Nx_Dsp_pen(1)+t-1),P_ct(t,Nx_Dsp_pen(1)+t-1)); % ���������ܶȣ�kg/m^3
    rho_l_ct(t,Nx_Dsp_pen(1)+t-1)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)+t-1),P_ct(t,Nx_Dsp_pen(1)+t-1)); % ���ڻ�Һ�ܶȣ�kg/m^3
    mu_g_ct(t,Nx_Dsp_pen(1)+t-1)=RheologyG(T_ct(t,Nx_Dsp_pen(1)+t-1),P_ct(t,Nx_Dsp_pen(1)+t-1)); % ��������ճ�ȣ�Pa*s
    mu_l_ct(t,Nx_Dsp_pen(1)+t-1)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(1)+t-1),P_ct(t,Nx_Dsp_pen(1)+t-1)); % ���ڻ�Һճ�ȣ�Pa*s
    alpha_f_ct(t,Nx_Dsp_pen(1)+t-1)=1; % ������ĭ����
    gamma_g_ct(t,Nx_Dsp_pen(1)+t-1)=(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)+t-1))/(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(1)+t-1)+Qm_l_0/rho_l_ct(t,Nx_Dsp_pen(1)+t-1)); % ������ĭ����
    gamma_l_ct(t,Nx_Dsp_pen(1)+t-1)=1-gamma_g_ct(t,Nx_Dsp_pen(1)+t-1); % ����Һ��������
    alpha_g_ct(t,Nx_Dsp_pen(1)+t-1)=alpha_f_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_g_ct(t,Nx_Dsp_pen(1)+t-1); % �������庬��
    alpha_l_ct(t,Nx_Dsp_pen(1)+t-1)=alpha_f_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_l_ct(t,Nx_Dsp_pen(1)+t-1); % ���ڻ�Һ����
    rho_f_ct(t,Nx_Dsp_pen(1)+t-1)=rho_g_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_g_ct(t,Nx_Dsp_pen(1)+t-1)+rho_l_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_l_ct(t,Nx_Dsp_pen(1)+t-1); % ������ĭ�ܶȣ�kg/m^3
    mu_f_ct(t,Nx_Dsp_pen(1)+t-1)=mu_g_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_g_ct(t,Nx_Dsp_pen(1)+t-1)+mu_l_ct(t,Nx_Dsp_pen(1)+t-1)*gamma_l_ct(t,Nx_Dsp_pen(1)+t-1); % ������ĭճ�ȣ�Pa*s
    V_f_ct(t,Nx_Dsp_pen(1)+t-1)=Qm_f_0/(rho_f_ct(t,Nx_Dsp_pen(1)+t-1)*A_ct(t,Nx_Dsp_pen(1)+t-1)); % �����������٣�m/s
    [Ff_ct(t,Nx_Dsp_pen(1)+t-1),f_ct(t,Nx_Dsp_pen(1)+t-1),Re_ct(t,Nx_Dsp_pen(1)+t-1),flow_pattern_ct(t,Nx_Dsp_pen(1)+t-1)]=Friction_coiledtubing(rho_f_ct(t,Nx_Dsp_pen(1)+t-1),V_f_ct(t,Nx_Dsp_pen(1)+t-1),mu_f_ct(t,Nx_Dsp_pen(1)+t-1),D_ct_i(t,Nx_Dsp_pen(1)+t-1),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    % ��Nx_Dsp_pen(1)+t-2��1���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)+t-2:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % ����ѹ������ֵ��Pa
        
        err_DriPipePressure=1; % ����ѹ��������
        COUNT_DriPipePressure=0; % ����ѹ������������ֵ
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % ����ѹ����������
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һ�ܶȣ�kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һճ�ȣ�Pa*s
            alpha_f_ct(t,x)=1; % ������ĭ����
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % ������ĭ����
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % �������庬��
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % ���ڻ�Һ����
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭճ�ȣ�Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % ������ĭ���٣�m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % ����ѹ����Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % �������ѹ������ֵ�����ֵ֮���������
            P_ct_ass(t,x)=P_ct(t,x); % �µĹ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(1)+t��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)+t:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % �������庬��
        alpha_l_ct(t,x)=alpha_l_a(t,x); % ���ڻ�Һ����
        gamma_g_ct(t,x)=gamma_g_a(t,x); % ������ĭ����
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
        rho_f_ct(t,x)=rho_m(t,x); % ���������ܶȣ�kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % ��������ճ�ȣ�Pa*s
        V_f_ct(t,t)=V_m(t,x); % �����������٣�sm/s
        Ff_ct(t,x)=Ff_a(t,x); % �������嵥λ����Ħ��ѹ����Pa/m��
        f_ct(t,x)=f_a(t,x); % �������巶��Ħ������
        Re_ct(t,x)=Re_a(t,x); % ����������ŵ����������̬
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % ����������̬
        P_ct(t,x)=P_a(t,x); % ����ѹ����Pa
    end
end

%% ��1��Ntʱ��ڵ��̹ܶγ���ѹ������ѹ���㣨�������ڣ�
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % ��������������٣�m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % ����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % �̹ܶ����嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % �����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % ���ڵ�ע��ͷ������Ħ��ѹ����Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % ��������Ħ��ѹ����Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % ����������Ͳ��Ħ��ѹ����Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % �̹ܶ�Ħ��ѹ����Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % �̹ܶγ���ѹ����Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % ��ѹ��Pa
end

%% ���㾮���ۻ���ɰ��
M_w_tem(1)=0; % ����˲ʱ��ɰ����ֵ��kg
M_w_tot(1)=0; % �����ۻ���ɰ����ֵ��kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮��˲ʱ��ɰ����kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮���ۻ���ɰ����kg
end

%% ����ĩ̬��Ͳ���ɰŨ��
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(2)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% ���㻷����ĭ����ƽ��ֵV_l_mre
V_f_mre=0; % ���շ���ƽ��ֵ��m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % ���շ���ƽ��ֵ��m/s
    end
end

%% ����������ĩ��ƽ��ֵVsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % ��Vsr��������ȡֵ����ĳ���ĩ�٣�m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % �������ĩ��ƽ��ֵ��m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % �������ĩ��ƽ��ֵ��m/s
end

%% �ж��Ƿ���Ч��ɰ�����շ���ƽ��ֵ����2��ɰ������ĩ��ƽ��ֵ�����Ƿ���ɳ�ɰ
fprintf("��ϴ�������2��\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % �����д���˵����ǰ�����¡���ɳ�ɰ��
    else
        fprintf("Sand Cleanout UnFinished!\n"); % �����д���˵����ǰ�����¡�δ��ɳ�ɰ��
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
end

%% ��ͲECD����
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % ����ECD��kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % ����ECD��kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
end

%% ��ϴ�������ĩ״̬����
alpha_f_a_5=alpha_f_a(Nt,:); % ��ϴ�������ĩ״̬��ĭ����
alpha_g_a_5=alpha_g_a(Nt,:); % ��ϴ�������ĩ״̬���ຬ��
alpha_l_a_5=alpha_l_a(Nt,:); % ��ϴ�������ĩ״̬Һ�ຬ��
alpha_s_5=alpha_s(Nt,:); % ��ϴ�������ĩ״̬���ຬ��
f_a_5=f_a(Nt,:); % ��ϴ�������ĩ״̬����Ħ������
Ff_a_5=Ff_a(Nt,:); % ��ϴ�������ĩ״̬���յ�λ����Ħ��ѹ����Pa/m
flow_pattern_a_5=flow_pattern_a(Nt,:); % ��ϴ�������ĩ״̬����������̬
gamma_g_a_5=gamma_g_a(Nt,:); % ��ϴ�������ĩ״̬��ĭ����
gamma_l_a_5=gamma_l_a(Nt,:); % ��ϴ�������ĩ״̬Һ��������
mu_f_a_5=mu_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭճ�ȣ�Pa*s
mu_g_a_5=mu_g_a(Nt,:); % ��ϴ�������ĩ״̬����ճ�ȣ�Pa*s
mu_l_a_5=mu_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ��ճ�ȣ�Pa*s
mu_m_5=mu_m(Nt,:); % ��ϴ�������ĩ״̬���ջ����ճ�ȣ�Pa*s
mu_s_5=mu_s(Nt,:); % ��ϴ�������ĩ״̬����ճ�ȣ�Pa*s
P_a_5=P_a(Nt,:); % ��ϴ�������ĩ״̬����ѹ����Pa
Re_a_5=Re_a(Nt,:); % ��ϴ�������ĩ״̬������ŵ��
rho_f_a_5=rho_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ�ܶȣ�kg/m^3
rho_g_a_5=rho_g_a(Nt,:); % ��ϴ�������ĩ״̬���������ܶȣ�kg/m^3
rho_l_a_5=rho_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ���ܶȣ�kg/m^3
rho_m_5=rho_m(Nt,:); % ��ϴ�������ĩ״̬���ջ�����ܶȣ�kg/m^3
V_f_a_5=V_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ���٣�m/s
V_g_a_5=V_g_a(Nt,:); % ��ϴ�������ĩ״̬�������٣�m/s
V_l_a_5=V_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ�����٣�m/s
V_m_5=V_m(Nt,:); % ��ϴ�������ĩ״̬���ջ�������٣�m/s
V_s_5=V_s(Nt,:); % ��ϴ�������ĩ״̬�������٣�m/s
Va_f_a_5=Va_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ������٣�m/s
Va_g_a_5=Va_g_a(Nt,:); % ��ϴ�������ĩ״̬���������٣�m/s
Va_l_a_5=Va_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ�������٣�m/s
Va_s_5=Va_s(Nt,:); % ��ϴ�������ĩ״̬���������٣�m/s
Vsr_5=Vsr(Nt,:); % ��ϴ�������ĩ״̬���໬���ٶȣ�m/s

%% ���ݴ洢����ϴ������̣�
ANS_Nt_5=Nt; % ʱ��ڵ���
ANS_Time_5=Time(1:Nt); % ��ϴʱ�䣬s
ANS_alpha_g_a_5=alpha_g_a(1:Nt,:); % ���������������
ANS_alpha_l_a_5=alpha_l_a(1:Nt,:); % ����Һ���������
ANS_alpha_s_5=alpha_s(1:Nt,:); % ���ຬ��
ANS_P_ct_5=P_ct(1:Nt,:); % ����ѹ����Pa��
ANS_delta_P_SWT_5=delta_P_SWT(1:Nt); % ��ɰ����ѹ����Pa��
ANS_P_a_5=P_a(1:Nt,:); % ����ѹ����Pa��
ANS_T_a_5=T_a(1:Nt,:); % �����¶ȣ��棩
ANS_T_ct_5=T_ct(1:Nt,:); % �����¶ȣ��棩
ANS_P_coil_5=P_coil(1:Nt); % �̹ܶγ���ѹ����Pa��
ANS_P_pump_5=P_pump(1:Nt); % ��ѹ��Pa��
ANS_M_w_tem_5=M_w_tem(1:Nt); % ����˲ʱ��ɰ����kg��
ANS_M_w_tot_5=M_w_tot(1:Nt); % �����ۻ���ɰ����kg��
ANS_Va_s_5=Va_s(1:Nt,:); % ��м����ٶȣ�m/s��
ANS_Va_f_a_5=Va_f_a(1:Nt,:); % ������ĭ������٣�m/s��
ANS_V_s_5=V_s(1:Nt,:); % ��м�����ٶȣ�m/s��
ANS_V_f_a_5=V_f_a(1:Nt,:); % ������ĭ���٣�m/s��
ANS_alpha_g_ct_5=alpha_g_ct(1:Nt,:); % ���������������
ANS_alpha_l_ct_5=alpha_l_ct(1:Nt,:); % ����Һ���������
ANS_gamma_g_ct_5=gamma_g_ct(1:Nt,:); % ������ĭ����
ANS_gamma_g_a_5=gamma_g_a(1:Nt,:); % ������ĭ����
ANS_ECD_a_5=ECD_a(1:Nt,:); % ����ECD��kg/m^3��
ANS_L_coil_5=L_coil(1:Nt); % ���������m��



%% ���϶������2�����ʵ�λ�ƣ�
%% ������ֵ����
dt=zeros(); % ʱ�䲽����s
dx=zeros(); % �ռ䲽����m
Time=zeros(); % ������������ʱ����s
L_coil=zeros(); % ���������m
L_reel=zeros(); % �̹ܶγ��ȣ�m
D_t_i=zeros(); % �͹��ھ���m
D_ct_o=zeros(); % �����͹��⾶��m
D_ct_i=zeros(); % �����͹��ھ���m
M_s=zeros(); % ���׽�ɰ����kg/s
P_a=zeros(); % ����ѹ����Pa
rho_g_a=zeros(); % �����ܶȣ�kg/m^3
rho_l_a=zeros(); % ����Һ���ܶȣ�kg/m^3
rho_f_a=zeros(); % ������ĭ�ܶȣ�kg/m^3
mu_g_a=zeros(); % ����ճ�ȣ�Pa*s
mu_l_a=zeros(); % ����Һ��ճ�ȣ�Pa*s
mu_f_a=zeros(); % ������ĭճ�ȣ�Pa*s
mu_s=zeros(); % ����ճ�ȣ�Pa*s
alpha_g_a=zeros(); % ���ຬ��
alpha_l_a=zeros(); % Һ�ຬ��
alpha_f_a=zeros(); % ��ĭ����
alpha_s=zeros(); % ���ຬ��
Va_g_a=zeros(); % ���������٣�m/s
Va_l_a=zeros(); % ����Һ�������٣�m/s
Va_f_a=zeros(); % ������ĭ������٣�m/s
Va_s=zeros(); % ���������٣�m/s
V_g_a=zeros(); % �������٣�m/s
V_l_a=zeros(); % ����Һ�����٣�m/s
V_f_a=zeros(); % ������ĭ���٣�m/s
V_s=zeros(); % �������٣�m/s
Vsr=zeros(); % ���໬���ٶȣ�m/s
gamma_g_a=zeros(); % ��ĭ����
gamma_l_a=zeros(); % Һ��������
V_m=zeros(); % ���ջ�������٣�m/s
rho_m=zeros(); % ���ջ�����ܶȣ�kg/m^3
mu_m=zeros(); % ���ջ����ճ�ȣ�Pa*s
Ff_a=zeros(); % ���յ�λ����Ħ��ѹ����Pa/m
f_a=zeros(); % ����Ħ������
Re_a=zeros(); % ������ŵ��
flow_pattern_a=zeros(); % ����������̬

%% ����ռ䲽��dx����Ӧʱ�䲽��dt
V6=(Dsp_pen(2)-Lp)/(3*3600); % ���϶����ٶȣ�m/s�����룩
t_6=(Dsp_pen(2)-Lp)/V6; % ���϶�����ʱ����s

Nt=Nx_Dsp_pen(2)-Nx_Lp+1; % ʱ��ڵ���
nt=Nt-1; % ʱ��������

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % �ռ䲽����m
end

Time(1)=0; % �����������ʱ����ֵ��s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Dsp_pen(2)-t)/V6; % ÿ���һ���ռ䲽������ʱ�䣬s
    Time(t+1)=Time(t)+dt(t); % �������������(t+1)���ռ�ڵ�����������ʱ����s
end

%% ���㲻ͬʱ��������������L_coil���̹ܶ������ܳ���L_reel
L=10000; % �����͹��ܳ���m�����룩
L_wg=8; % ���ڵ�ע��ͷ�����������ܳ��ȣ�m�����룩
L_goose=3; % �������������ܳ��ȣ�m�����룩
D_goose=2; % �������ΰ뾶��m�����룩
H_goose=10; % �����������߶ȣ�m�����룩
L_gr=20; % ����������Ͳ�������ܳ��ȣ�m�����룩
theta_gr=acosd(H_goose/L_gr); % ����������Ͳ����������Ǧ���߼нǣ���
D_r_i=3; % ��Ͳ�ھ���m�����룩
D_r_o=5; % ��Ͳ�⾶��m�����룩
W_r=5; % ��Ͳ��ȣ�m�����룩
D_cable=0.005; % �����⾶��m

L_coil(1)=Dsp_pen(2); % ��ʼʱ�������ܵײ���ȣ�m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % ��ʼʱ���̹ܶ������ܳ��ȣ�m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)-dx(Nx_Dsp_pen(2)-t+1); % �����������ȣ�m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % �̹ܶ������ܳ��ȣ�m
end
L_cable=L_coil; % �������볤�ȣ�m

%% �������뼰Ԥ����
D_ct_o_0=0.04445; % �����͹��⾶��m�����룩
L1=2000; % �����ڵ�һ�������ܳ��ȣ�m�����룩
D_ct_i_1=0.03709; % �����ڵ�һ���������ھ���m�����룩
L2=2000; % �����ڵڶ��������ܳ��ȣ�m�����룩
D_ct_i_2=0.03653; % �����ڵڶ����������ھ���m�����룩
L3=2000; % �����ڵ����������ܳ��ȣ�m�����룩
D_ct_i_3=0.03555; % �����ڵ������������ھ���m�����룩
L4=2000; % �����ڵ��Ķ������ܳ��ȣ�m�����룩
D_ct_i_4=0.03489; % �����ڵ��Ķ��������ھ���m�����룩
L5=2000; % �����ڵ���������ܳ��ȣ�m�����룩
D_ct_i_5=0.03409; % �����ڵ�����������ھ���m�����룩

L_t_1=4000; % �ϲ��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_1=0.09718; %0.068;%0.09718; % �ϲ��͹ܣ����׹ܻ����ۣ��ھ���m�����룩
L_t_2=2200; % �²��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_2=0.09718; %0.13970; % �²��͹ܣ����׹ܻ����ۣ��ھ���m�����룩

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % �͹ܣ����׹ܻ����ۣ��ھ���m
        else
            D_t_i(t,x)=D_t_i_2; % �͹ܣ����׹ܻ����ۣ��ھ���m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(2)
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % �����͹��⾶��m
            else
                D_ct_o(t,x)=0; % �����͹��⾶��m
            end
        else
            D_ct_o(t,x)=0; % �����͹��⾶��m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Dsp_pen(2)
            if Depth(x)<=L_coil(t) % �����ļ�����������ھ�
                if (L_coil(t)-Depth(x))<=L1 && (L_coil(t)-Depth(x))>=0
                    D_ct_i(t,x)=D_ct_i_1; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % �����͹��ھ���m
                end
            else
                D_ct_i(t,x)=0; % �����͹��ھ���m
            end
        else
            D_ct_i(t,x)=0; % �����͹��ھ���m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % �����ļ������������ֵ
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % �����͹��ڽ������m^2
        else
            A_ct(t,x)=0; % �����͹��ڽ������m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % ����ˮ��ֱ����m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % ���ս������m^2
    end
end

%% �������뼰Ԥ����
h_t=25.4*10^(-6); % �͹ܣ����׹ܻ����ۣ����Դֲڶȣ�m�����룩
h_ct=25.4*10^(-6); % �����͹ܾ��Դֲڶȣ�m�����룩
h_a=(h_t+h_ct)/2; % ����ƽ�����Դֲڶȣ�m
epsilon_e=1*10^(-3); % �����������ޣ����룩
epsilon_t=1*10^3; % ���������������룩
g=9.81; % �������ٶȣ�m/s^2��Ĭ�ϣ�

T_0=20; % ��Һ�����¶ȣ��棨���룩
P_0=0.1*10^6; % ��Һ����ѹ����Pa�����룩
rho_l_0=1150; % T_0��P_0�»�Һ�ܶȣ�kg/m^3�����룩
mu_l_0=0.03; % T_0��P_0�»�Һճ�ȣ�Pa*s�����룩
Qv_l_0=0.1/60; % ��Һ���������m^3/s�����룩
Qm_l_0=Qv_l_0*rho_l_0; % ��Һ����������kg/s
rho_g_0=0.655; % ע�����ܶȣ�kg/m^3
mu_g_0=RheologyG(T_0,P_0); % ����ճ�ȣ�Pa*s
Qv_g_0=4/60; % ע�������������m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % ע��������������kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % ��ĭ����������kg/s

D_s=1*10^(-3); % ɰ��ֱ����m�����룩
rho_s=2000; % ɰ���ܶȣ�kg/m^3�����룩
H_s=L_s_b-Nx_Dsp_pen(1); % �ײ�ɰ���߶ȣ�m
PHI=0.6; % ɰ����ӯ�ȣ����룩
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % ����ɰ����������kg

D_nozzle=4/1000; % ����ֱ����m�����룩
N_nozzle=3; % ������������룩
C=0.95; % ��������ϵ����ȡ0.95�����룩

C0=1.2; % Ư��������ϵ����Ĭ�ϣ�

M_s(1)=0; % ���׽�ɰ����kg/s
for t=2:1:Nt
    M_s(t)=0; % ���׽�ɰ����kg/s
end

OutPressure=1*10^6; % ����ѹ����Pa�����룩

%% �¶�����
T_i=20; % ��ĭע���¶ȣ���
T_g=0.02; % �����ݶȣ���/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % ����������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������
    end
end
T_a=T_ct; % ������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������

%% ��1��ʱ��ڵ㣨��ʼʱ�̣���ز�����ֵ���㣨���գ�
for x=1:1:Nx
    P_a(1,x)=P_a_5(x); % ����ѹ����Pa
    rho_g_a(1,x)=rho_g_a_5(x); % �����ܶȣ�kg/m^3
    rho_l_a(1,x)=rho_l_a_5(x); % ����Һ���ܶȣ�kg/m^3
    rho_f_a(1,x)=rho_f_a_5(x); % ������ĭ�ܶȣ�kg/m^3
    mu_g_a(1,x)=mu_g_a_5(x); % ����ճ�ȣ�Pa*s
    mu_l_a(1,x)=mu_l_a_5(x); % ����Һ��ճ�ȣ�Pa*s
    mu_f_a(1,x)=mu_f_a_5(x); % ������ĭճ�ȣ�Pa*s
    mu_s(1,x)=mu_s_5(x); % ����ճ�ȣ�Pa*s
    alpha_g_a(1,x)=alpha_g_a_5(x); % ���ຬ��
    alpha_l_a(1,x)=alpha_l_a_5(x); % Һ�ຬ��
    alpha_f_a(1,x)=alpha_f_a_5(x); % ��ĭ����
    alpha_s(1,x)=alpha_s_5(1,x); % ���ຬ��
    Va_g_a(1,x)=Va_g_a_5(x); % ���������٣�m/s
    Va_l_a(1,x)=Va_l_a_5(x); % ����Һ�������٣�m/s
    Va_f_a(1,x)=Va_f_a_5(x); % ������ĭ������٣�m/s
    Va_s(1,x)=Va_s_5(x); % ���������٣�m/s
    V_g_a(1,x)=V_g_a_5(x); % �������٣�m/s
    V_l_a(1,x)=V_l_a_5(x); % ����Һ�����٣�m/s
    V_f_a(1,x)=V_f_a_5(x); % ������ĭ���٣�m/s
    V_s(1,x)=V_s_5(x); % �������٣�m/s
    Vsr(1,x)=Vsr_5(x); % ���໬���ٶȣ�m/s
    gamma_g_a(1,x)=gamma_g_a_5(x); % ��ĭ����
    gamma_l_a(1,x)=gamma_l_a_5(1,x); % Һ��������
    V_m(1,x)=V_m_5(x); % ���ջ�������٣�m/s
    rho_m(1,x)=rho_m_5(x); % ���ջ�����ܶȣ�kg/m^3
    mu_m(1,x)=mu_m_5(x); % ���ջ����ճ�ȣ�Pa*s
    Ff_a(1,x)=Ff_a_5(x); % ���յ�λ����Ħ��ѹ����Pa/m
    f_a(1,x)=f_a_5(x); % ����Ħ������
    Re_a(1,x)=Re_a_5(x); % ������ŵ��
    flow_pattern_a(1,x)=flow_pattern_a_5(x); % ����������̬
end

%% ��2��Nt��ʱ��ڵ���ز������㣨���գ�
for t=2:1:Nt
    P_a(t,Nx_Dsp_pen(2)-t+1)=P_a(t-1,Nx_Dsp_pen(2)-t+1);  % ���չܵ�ѹ������ֵ��Pa
    
    err_OutPressure=1; % ����ѹ��������
    COUNT_OutPressure=0; % ����ѹ������������ֵ
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % ����ѹ����������
        
        % �ܵף���Nx_Dsp_pen(2)-t+1���ռ�ڵ㣩����ز�������
        rho_g_a(t,Nx_Dsp_pen(2)-t+1)=DensityG(T_a(t,Nx_Dsp_pen(2)-t+1),P_a(t,Nx_Dsp_pen(2)-t+1)); % ���������ܶȣ�kg/m^3
        rho_l_a(t,Nx_Dsp_pen(2)-t+1)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(2)-t+1),P_a(t,Nx_Dsp_pen(2)-t+1)); % ���ջ�Һ�ܶȣ�kg/m^3
        mu_g_a(t,Nx_Dsp_pen(2)-t+1)=RheologyG(T_a(t,Nx_Dsp_pen(2)-t+1),P_a(t,Nx_Dsp_pen(2)-t+1)); % ��������ճ�ȣ�Pa*s
        mu_l_a(t,Nx_Dsp_pen(2)-t+1)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Dsp_pen(2)-t+1),P_a(t,Nx_Dsp_pen(2)-t+1)); % ���ջ�Һճ�ȣ�Pa*s
        gamma_g_a(t,Nx_Dsp_pen(2)-t+1)=(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(2)-t+1))/(Qm_g_0/rho_g_a(t,Nx_Dsp_pen(2)-t+1)+Qm_l_0/rho_l_a(t,Nx_Dsp_pen(2)-t+1)); % ��ĭ����
        gamma_l_a(t,Nx_Dsp_pen(2)-t+1)=1-gamma_g_a(t,Nx_Dsp_pen(2)-t+1); % Һ��������
        alpha_g_a(t,Nx_Dsp_pen(2)-t+1)=alpha_f_a(t-1,Nx_Dsp_pen(2)-t+1)*gamma_g_a(t,Nx_Dsp_pen(2)-t+1); % �������庬��
        alpha_l_a(t,Nx_Dsp_pen(2)-t+1)=alpha_f_a(t-1,Nx_Dsp_pen(2)-t+1)*gamma_l_a(t,Nx_Dsp_pen(2)-t+1); % ���ջ�Һ����
        rho_f_a(t,Nx_Dsp_pen(2)-t+1)=rho_g_a(t,Nx_Dsp_pen(2)-t+1)*gamma_g_a(t,Nx_Dsp_pen(2)-t+1)+rho_l_a(t,Nx_Dsp_pen(2)-t+1)*gamma_l_a(t,Nx_Dsp_pen(2)-t+1); % ������ĭ�ܶȣ�kg/m^3
        mu_f_a(t,Nx_Dsp_pen(2)-t+1)=mu_g_a(t,Nx_Dsp_pen(2)-t+1)*gamma_g_a(t,Nx_Dsp_pen(2)-t+1)+mu_l_a(t,Nx_Dsp_pen(2)-t+1)*gamma_l_a(t,Nx_Dsp_pen(2)-t+1); % ������ĭճ�ȣ�Pa*s
        mu_s(t,Nx_Dsp_pen(2)-t+1)=mu_f_a(t,Nx_Dsp_pen(2)-t+1); % ����ճ�ȣ�Pa*s
        Va_f_a(t,Nx_Dsp_pen(2)-t+1)=Qm_f_0/(A_a(t,Nx_Dsp_pen(2)-t+1)*rho_f_a(t,Nx_Dsp_pen(2)-t+1)); % ������ĭ������٣�m/s
        Va_g_a(t,Nx_Dsp_pen(2)-t+1)=Va_f_a(t,Nx_Dsp_pen(2)-t+1); % �������������٣�m/s
        Va_l_a(t,Nx_Dsp_pen(2)-t+1)=Va_f_a(t,Nx_Dsp_pen(2)-t+1); % ���ջ�Һ������٣�m/s
        Va_s(t,Nx_Dsp_pen(2)-t+1)=M_s(t)/(A_a(t,Nx_Dsp_pen(2)-t+1)*rho_s); % ��м����ٶȣ�m/s
        Vsr(t,Nx_Dsp_pen(2)-t+1)=12*(mu_f_a(t,Nx_Dsp_pen(2)-t+1)/(rho_f_a(t,Nx_Dsp_pen(2)-t+1)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Dsp_pen(2)-t+1))/rho_f_a(t,Nx_Dsp_pen(2)-t+1))*((rho_f_a(t,Nx_Dsp_pen(2)-t+1)*D_s/mu_f_a(t,Nx_Dsp_pen(2)-t+1))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
        alpha_s(t,Nx_Dsp_pen(2)-t+1)=Va_s(t,Nx_Dsp_pen(2)-t+1)/(C0*(Va_s(t,Nx_Dsp_pen(2)-t+1)+Va_f_a(t,Nx_Dsp_pen(2)-t+1))-Vsr(t,Nx_Dsp_pen(2)-t+1));  % �����������
        V_s(t,Nx_Dsp_pen(2)-t+1)=0;%Va_s(t,Nx_Dsp_pen(2)-t+1)/alpha_s(t,Nx_Dsp_pen(2)-t+1); % ��м�ٶȣ�m/s
        
        alpha_f_a(t,Nx_Dsp_pen(2)-t+1)=1-alpha_s(t,Nx_Dsp_pen(2)-t+1); % ������ĭ����
        alpha_g_a(t,Nx_Dsp_pen(2)-t+1)=alpha_f_a(t,Nx_Dsp_pen(2)-t+1)*gamma_g_a(t,Nx_Dsp_pen(2)-t+1); % �������庬��
        alpha_l_a(t,Nx_Dsp_pen(2)-t+1)=alpha_f_a(t,Nx_Dsp_pen(2)-t+1)*gamma_l_a(t,Nx_Dsp_pen(2)-t+1); % ���ջ�Һ����
        V_f_a(t,Nx_Dsp_pen(2)-t+1)=Va_f_a(t,Nx_Dsp_pen(2)-t+1)/alpha_f_a(t,Nx_Dsp_pen(2)-t+1); % ������ĭ���٣�m/s
        V_g_a(t,Nx_Dsp_pen(2)-t+1)=V_f_a(t,Nx_Dsp_pen(2)-t+1); % �����������٣�m/s
        V_l_a(t,Nx_Dsp_pen(2)-t+1)=V_f_a(t,Nx_Dsp_pen(2)-t+1); % ���ջ�Һ���٣�m/s
        
        V_m(t,Nx_Dsp_pen(2)-t+1)=Va_s(t,Nx_Dsp_pen(2)-t+1)+Va_f_a(t,Nx_Dsp_pen(2)-t+1); % ���ջ�����ٶȣ�m/s
        rho_m(t,Nx_Dsp_pen(2)-t+1)=alpha_s(t,Nx_Dsp_pen(2)-t+1)*rho_s+alpha_f_a(t,Nx_Dsp_pen(2)-t+1)*rho_f_a(t,Nx_Dsp_pen(2)-t+1); % ���ջ�����ܶȣ�kg/m^3
        mu_m(t,Nx_Dsp_pen(2)-t+1)=alpha_s(t,Nx_Dsp_pen(2)-t+1)*mu_s(t,Nx_Dsp_pen(2)-t+1)+alpha_f_a(t,Nx_Dsp_pen(2)-t+1)*mu_f_a(t,Nx_Dsp_pen(2)-t+1); % ���ջ����ճ�ȣ�Pa*s
        [Ff_a(t,Nx_Dsp_pen(2)-t+1),f_a(t,Nx_Dsp_pen(2)-t+1),Re_a(t,Nx_Dsp_pen(2)-t+1),flow_pattern_a(t,Nx_Dsp_pen(2)-t+1)]=Friction_annulus(rho_m(t,Nx_Dsp_pen(2)-t+1),V_m(t,Nx_Dsp_pen(2)-t+1),mu_m(t,Nx_Dsp_pen(2)-t+1),D_h(t,Nx_Dsp_pen(2)-t+1),h_a,rho_f_a(t,Nx_Dsp_pen(2)-t+1),V_f_a(t,Nx_Dsp_pen(2)-t+1)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
        
        % ��Nx_Dsp_pen(2)-t��1���ռ�ڵ㴦��ز�������
        for x=Nx_Dsp_pen(2)-t:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % ����ѹ������ֵ��Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % ���ջ�Һ����
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            
            err_NodePressure=1; % ����ѹ��������
            COUNT_NodePressure=0; % ����ѹ������������ֵ
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % ����ѹ����������
                
                % ���߷��������������
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % ���������������ֵ1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % ���������������ֵ2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % ������������������
                COUNT_NodeEg=0; % ��������
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % �����������Ϊalpha_s_ass1ʱ������ֵ��
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % �����ٶȣ�m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % ���������٣�m/s
                    alpha_f_ass1=1-alpha_s_ass1; % ��ĭ�������
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % ��ĭ�ٶȣ�m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % ��ĭ������٣�m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % ���ջ�����ٶȣ�m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % ���������������ֵ
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % ����ĺ��������Ľ������ʵ�����������
                    
                    % �����������Ϊalpha_s_ass2ʱ������ֵ��
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % ���������غ㷽����ɢ��ʽ�м�ֵ����
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % �����ٶȣ�m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % ���������٣�m/s
                    alpha_f_ass2=1-alpha_s_ass2; % ��ĭ�������
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % ��ĭ�ٶȣ�m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % ��ĭ������٣�m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % ���ջ�����ٶȣ�m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % ���������������ֵ
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % ����ĺ��������Ľ������ʵ�����������
                    
                    % ���߷���������������
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % �µĹ��������������ֵalpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % ������������������
                    alpha_s_ass1=alpha_s_ass2; % �µĹ��������������ֵ1
                    alpha_s_ass2=alpha_s_ass3; % �µĹ��������������ֵ2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % ���������õ�����ʵ�����������ֵ����alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % �������������С��һ��ֵʱ����Ϊ���������Ϊ0�����ڷ�ֹ���ֺ����ļ������
                    alpha_s(t,x)=0; % �����������
                    V_s(t,x)=0; % �����ٶȣ�m/s
                    Va_s(t,x)=0; % �������ٶȣ�m/s
                    mu_s(t,x)=0; % ����ճ�ȣ�Pa*s
                    Vsr(t,x)=0; % ��м����ĩ�٣�m/s
                    
                    alpha_f_a(t,x)=1; % ��ĭ�������                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                else
                    V_s(t,x)=V_s_ass2; % �����ٶȣ�m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % ���������٣�m/s
                    Vsr(t,x)=Vsr_ass2; % ��м����ĩ�٣�m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % ��ĭ�ܶȣ�kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % ��ĭ�������
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    V_f_a(t,x)=V_f_ass2; % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % ���ջ�����ٶȣ�m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % ����ѹ������ֵ��Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % ����ѹ��������
                P_a(t,x)=P_new; % �µĻ���ѹ������ֵ��Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % ����ѹ��������
        if (P_a(t,1)-OutPressure)>0 % ���ݳ���ѹ�������������Ի��չܵ�ѹ������ֵ���е���
            P_a(t,Nx_Dsp_pen(2)-t+1)=P_a(t,Nx_Dsp_pen(2)-t+1)-(P_a(t,1)-OutPressure)/2; % �µĻ��չܵ�ѹ������ֵ��Pa
        else
            P_a(t,Nx_Dsp_pen(2)-t+1)=P_a(t,Nx_Dsp_pen(2)-t+1)-(P_a(t,1)-OutPressure)/2*0.3; % �µĻ��չܵ�ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(2)-t+2��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(2)-t+2:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=0; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
end

%% ��1��Ntʱ��ڵ���ͷѹ������
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Dsp_pen(2)-t+1); % ������ĭ���������m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Dsp_pen(2)-t+1)); % �����ɰ����ѹ����Pa����������������V_nozzle��m/s��
end

%% ��1��Ntʱ��ڵ���ز������㣨�������ڣ�
for t=1:1:Nt
    % �ܵף���Nx_Dsp_pen(2)-t+1���ռ�ڵ㣩����ز�������
    P_ct(t,Nx_Dsp_pen(2)-t+1)=P_a(t,Nx_Dsp_pen(2)-t+1)+delta_P_SWT(t); % ����ѹ����Pa
    rho_g_ct(t,Nx_Dsp_pen(2)-t+1)=DensityG(T_ct(t,Nx_Dsp_pen(2)-t+1),P_ct(t,Nx_Dsp_pen(2)-t+1)); % ���������ܶȣ�kg/m^3
    rho_l_ct(t,Nx_Dsp_pen(2)-t+1)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(2)-t+1),P_ct(t,Nx_Dsp_pen(2)-t+1)); % ���ڻ�Һ�ܶȣ�kg/m^3
    mu_g_ct(t,Nx_Dsp_pen(2)-t+1)=RheologyG(T_ct(t,Nx_Dsp_pen(2)-t+1),P_ct(t,Nx_Dsp_pen(2)-t+1)); % ��������ճ�ȣ�Pa*s
    mu_l_ct(t,Nx_Dsp_pen(2)-t+1)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Dsp_pen(2)-t+1),P_ct(t,Nx_Dsp_pen(2)-t+1)); % ���ڻ�Һճ�ȣ�Pa*s
    alpha_f_ct(t,Nx_Dsp_pen(2)-t+1)=1; % ������ĭ����
    gamma_g_ct(t,Nx_Dsp_pen(2)-t+1)=(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(2)-t+1))/(Qm_g_0/rho_g_ct(t,Nx_Dsp_pen(2)-t+1)+Qm_l_0/rho_l_ct(t,Nx_Dsp_pen(2)-t+1)); % ������ĭ����
    gamma_l_ct(t,Nx_Dsp_pen(2)-t+1)=1-gamma_g_ct(t,Nx_Dsp_pen(2)-t+1); % ����Һ��������
    alpha_g_ct(t,Nx_Dsp_pen(2)-t+1)=alpha_f_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(2)-t+1); % �������庬��
    alpha_l_ct(t,Nx_Dsp_pen(2)-t+1)=alpha_f_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(2)-t+1); % ���ڻ�Һ����
    rho_f_ct(t,Nx_Dsp_pen(2)-t+1)=rho_g_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(2)-t+1)+rho_l_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(2)-t+1); % ������ĭ�ܶȣ�kg/m^3
    mu_f_ct(t,Nx_Dsp_pen(2)-t+1)=mu_g_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_g_ct(t,Nx_Dsp_pen(2)-t+1)+mu_l_ct(t,Nx_Dsp_pen(2)-t+1)*gamma_l_ct(t,Nx_Dsp_pen(2)-t+1); % ������ĭճ�ȣ�Pa*s
    V_f_ct(t,Nx_Dsp_pen(2)-t+1)=Qm_f_0/(rho_f_ct(t,Nx_Dsp_pen(2)-t+1)*A_ct(t,Nx_Dsp_pen(2)-t+1)); % �����������٣�m/s
    [Ff_ct(t,Nx_Dsp_pen(2)-t+1),f_ct(t,Nx_Dsp_pen(2)-t+1),Re_ct(t,Nx_Dsp_pen(2)-t+1),flow_pattern_ct(t,Nx_Dsp_pen(2)-t+1)]=Friction_coiledtubing(rho_f_ct(t,Nx_Dsp_pen(2)-t+1),V_f_ct(t,Nx_Dsp_pen(2)-t+1),mu_f_ct(t,Nx_Dsp_pen(2)-t+1),D_ct_i(t,Nx_Dsp_pen(2)-t+1),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    % ��Nx_Dsp_pen(2)-t��1���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(2)-t:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % ����ѹ������ֵ��Pa
        
        err_DriPipePressure=1; % ����ѹ��������
        COUNT_DriPipePressure=0; % ����ѹ������������ֵ
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % ����ѹ����������
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һ�ܶȣ�kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һճ�ȣ�Pa*s
            alpha_f_ct(t,x)=1; % ������ĭ����
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % ������ĭ����
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % �������庬��
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % ���ڻ�Һ����
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭճ�ȣ�Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % ������ĭ���٣�m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % ����ѹ����Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % �������ѹ������ֵ�����ֵ֮���������
            P_ct_ass(t,x)=P_ct(t,x); % �µĹ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(2)-t+2��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(2)-t+2:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % �������庬��
        alpha_l_ct(t,x)=alpha_l_a(t,x); % ���ڻ�Һ����
        gamma_g_ct(t,x)=gamma_g_a(t,x); % ������ĭ����
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
        rho_f_ct(t,x)=rho_m(t,x); % ���������ܶȣ�kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % ��������ճ�ȣ�Pa*s
        V_f_ct(t,t)=V_m(t,x); % �����������٣�sm/s
        Ff_ct(t,x)=Ff_a(t,x); % �������嵥λ����Ħ��ѹ����Pa/m��
        f_ct(t,x)=f_a(t,x); % �������巶��Ħ������
        Re_ct(t,x)=Re_a(t,x); % ����������ŵ����������̬
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % ����������̬
        P_ct(t,x)=P_a(t,x); % ����ѹ����Pa
    end
end

%% ��1��Ntʱ��ڵ��̹ܶγ���ѹ������ѹ���㣨�������ڣ�
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % ��������������٣�m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % ����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % �̹ܶ����嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % �����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % ���ڵ�ע��ͷ������Ħ��ѹ����Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % ��������Ħ��ѹ����Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % ����������Ͳ��Ħ��ѹ����Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % �̹ܶ�Ħ��ѹ����Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % �̹ܶγ���ѹ����Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % ��ѹ��Pa
end

%% ���㾮���ۻ���ɰ��
M_w_tem(1)=0; % ����˲ʱ��ɰ����ֵ��kg
M_w_tot(1)=0; % �����ۻ���ɰ����ֵ��kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮��˲ʱ��ɰ����kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮���ۻ���ɰ����kg
end

%% ����ĩ̬��Ͳ���ɰŨ��
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(1)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% ���㻷����ĭ����ƽ��ֵV_l_mre
V_f_mre=0; % ���շ���ƽ��ֵ��m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % ���շ���ƽ��ֵ��m/s
    end
end

%% ����������ĩ��ƽ��ֵVsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % ��Vsr��������ȡֵ����ĳ���ĩ�٣�m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % �������ĩ��ƽ��ֵ��m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % �������ĩ��ƽ��ֵ��m/s
end

%% �ж��Ƿ���Ч��ɰ�����շ���ƽ��ֵ����2��ɰ������ĩ��ƽ��ֵ�����Ƿ���ɳ�ɰ
fprintf("���϶������2��\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % �����д���˵����ǰ�����¡���ɳ�ɰ��
    else
        fprintf("Sand Cleanout UnFinished!\n"); % �����д���˵����ǰ�����¡�δ��ɳ�ɰ��
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
end

%% ��ͲECD����
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % ����ECD��kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % ����ECD��kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
end

%% ���϶������ĩ״̬����
alpha_f_a_6=alpha_f_a(Nt,:); % ��ϴ�������ĩ״̬��ĭ����
alpha_g_a_6=alpha_g_a(Nt,:); % ��ϴ�������ĩ״̬���ຬ��
alpha_l_a_6=alpha_l_a(Nt,:); % ��ϴ�������ĩ״̬Һ�ຬ��
alpha_s_6=alpha_s(Nt,:); % ��ϴ�������ĩ״̬���ຬ��
f_a_6=f_a(Nt,:); % ��ϴ�������ĩ״̬����Ħ������
Ff_a_6=Ff_a(Nt,:); % ��ϴ�������ĩ״̬���յ�λ����Ħ��ѹ����Pa/m
flow_pattern_a_6=flow_pattern_a(Nt,:); % ��ϴ�������ĩ״̬����������̬
gamma_g_a_6=gamma_g_a(Nt,:); % ��ϴ�������ĩ״̬��ĭ����
gamma_l_a_6=gamma_l_a(Nt,:); % ��ϴ�������ĩ״̬Һ��������
mu_f_a_6=mu_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭճ�ȣ�Pa*s
mu_g_a_6=mu_g_a(Nt,:); % ��ϴ�������ĩ״̬����ճ�ȣ�Pa*s
mu_l_a_6=mu_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ��ճ�ȣ�Pa*s
mu_m_6=mu_m(Nt,:); % ��ϴ�������ĩ״̬���ջ����ճ�ȣ�Pa*s
mu_s_6=mu_s(Nt,:); % ��ϴ�������ĩ״̬����ճ�ȣ�Pa*s
P_a_6=P_a(Nt,:); % ��ϴ�������ĩ״̬����ѹ����Pa
Re_a_6=Re_a(Nt,:); % ��ϴ�������ĩ״̬������ŵ��
rho_f_a_6=rho_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ�ܶȣ�kg/m^3
rho_g_a_6=rho_g_a(Nt,:); % ��ϴ�������ĩ״̬���������ܶȣ�kg/m^3
rho_l_a_6=rho_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ���ܶȣ�kg/m^3
rho_m_6=rho_m(Nt,:); % ��ϴ�������ĩ״̬���ջ�����ܶȣ�kg/m^3
V_f_a_6=V_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ���٣�m/s
V_g_a_6=V_g_a(Nt,:); % ��ϴ�������ĩ״̬�������٣�m/s
V_l_a_6=V_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ�����٣�m/s
V_m_6=V_m(Nt,:); % ��ϴ�������ĩ״̬���ջ�������٣�m/s
V_s_6=V_s(Nt,:); % ��ϴ�������ĩ״̬�������٣�m/s
Va_f_a_6=Va_f_a(Nt,:); % ��ϴ�������ĩ״̬������ĭ������٣�m/s
Va_g_a_6=Va_g_a(Nt,:); % ��ϴ�������ĩ״̬���������٣�m/s
Va_l_a_6=Va_l_a(Nt,:); % ��ϴ�������ĩ״̬����Һ�������٣�m/s
Va_s_6=Va_s(Nt,:); % ��ϴ�������ĩ״̬���������٣�m/s
Vsr_6=Vsr(Nt,:); % ��ϴ�������ĩ״̬���໬���ٶȣ�m/s

%% ���ݴ洢�����϶�����̣�
ANS_Nt_6=Nt; % ʱ��ڵ���
ANS_Time_6=Time(1:Nt); % ���������ʱ�䣬s
ANS_alpha_g_a_6=alpha_g_a(1:Nt,:); % ���������������
ANS_alpha_l_a_6=alpha_l_a(1:Nt,:); % ����Һ���������
ANS_alpha_s_6=alpha_s(1:Nt,:); % ���ຬ��
ANS_P_ct_6=P_ct(1:Nt,:); % ����ѹ����Pa��
ANS_delta_P_SWT_6=delta_P_SWT(1:Nt); % ��ɰ����ѹ����Pa��
ANS_P_a_6=P_a(1:Nt,:); % ����ѹ����Pa��
ANS_T_a_6=T_a(1:Nt,:); % �����¶ȣ��棩
ANS_T_ct_6=T_ct(1:Nt,:); % �����¶ȣ��棩
ANS_P_coil_6=P_coil(1:Nt); % �̹ܶγ���ѹ����Pa��
ANS_P_pump_6=P_pump(1:Nt); % ��ѹ��Pa��
ANS_M_w_tem_6=M_w_tem(1:Nt); % ����˲ʱ��ɰ����kg��
ANS_M_w_tot_6=M_w_tot(1:Nt); % �����ۻ���ɰ����kg��
ANS_Va_s_6=Va_s(1:Nt,:); % ��м����ٶȣ�m/s��
ANS_Va_f_a_6=Va_f_a(1:Nt,:); % ������ĭ������٣�m/s��
ANS_V_s_6=V_s(1:Nt,:); % ��м�����ٶȣ�m/s��
ANS_V_f_a_6=V_f_a(1:Nt,:); % ������ĭ���٣�m/s��
ANS_alpha_g_ct_6=alpha_g_ct(1:Nt,:); % ���������������
ANS_alpha_l_ct_6=alpha_l_ct(1:Nt,:); % ����Һ���������
ANS_gamma_g_ct_6=gamma_g_ct(1:Nt,:); % ������ĭ����
ANS_gamma_g_a_6=gamma_g_a(1:Nt,:); % ������ĭ����
ANS_ECD_a_6=ECD_a(1:Nt,:); % ����ECD��kg/m^3��
ANS_L_coil_6=L_coil(1:Nt); % ���������m��



%% ����ѭ������2�����ʵ�λ�ƣ�
%% ������ֵ����
dt=zeros(); % ʱ�䲽����s
dx=zeros(); % �ռ䲽����m
Time=zeros(); % ������������ʱ����s
L_coil=zeros(); % ���������m
L_reel=zeros(); % �̹ܶγ��ȣ�m
D_t_i=zeros(); % �͹��ھ���m
D_ct_o=zeros(); % �����͹��⾶��m
D_ct_i=zeros(); % �����͹��ھ���m
M_s=zeros(); % ���׽�ɰ����kg/s
P_a=zeros(); % ����ѹ����Pa
rho_g_a=zeros(); % �����ܶȣ�kg/m^3
rho_l_a=zeros(); % ����Һ���ܶȣ�kg/m^3
rho_f_a=zeros(); % ������ĭ�ܶȣ�kg/m^3
mu_g_a=zeros(); % ����ճ�ȣ�Pa*s
mu_l_a=zeros(); % ����Һ��ճ�ȣ�Pa*s
mu_f_a=zeros(); % ������ĭճ�ȣ�Pa*s
mu_s=zeros(); % ����ճ�ȣ�Pa*s
alpha_g_a=zeros(); % ���ຬ��
alpha_l_a=zeros(); % Һ�ຬ��
alpha_f_a=zeros(); % ��ĭ����
alpha_s=zeros(); % ���ຬ��
Va_g_a=zeros(); % ���������٣�m/s
Va_l_a=zeros(); % ����Һ�������٣�m/s
Va_f_a=zeros(); % ������ĭ������٣�m/s
Va_s=zeros(); % ���������٣�m/s
V_g_a=zeros(); % �������٣�m/s
V_l_a=zeros(); % ����Һ�����٣�m/s
V_f_a=zeros(); % ������ĭ���٣�m/s
V_s=zeros(); % �������٣�m/s
Vsr=zeros(); % ���໬���ٶȣ�m/s
gamma_g_a=zeros(); % ��ĭ����
gamma_l_a=zeros(); % Һ��������
V_m=zeros(); % ���ջ�������٣�m/s
rho_m=zeros(); % ���ջ�����ܶȣ�kg/m^3
mu_m=zeros(); % ���ջ����ճ�ȣ�Pa*s
Ff_a=zeros(); % ���յ�λ����Ħ��ѹ����Pa/m
f_a=zeros(); % ����Ħ������
Re_a=zeros(); % ������ŵ��
flow_pattern_a=zeros(); % ����������̬

%% ����ռ䲽��dx����Ӧʱ�䲽��dt
t_7=3600*1; % ����ѭ����ʱ����s

Nt=101; % ʱ��ڵ��������룩
nt=Nt-1; % ʱ��������

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % �ռ䲽����m
end

for t=1:1:Nt-1
    dt(t)=t_7/nt; % ʱ�䲽��
end

Time(1)=0; % ����ѭ����ʱ����s
for t=1:1:Nt-1
    Time(t+1)=Time(t)+dt(t); % ����ѭ������(t+1)��ʱ������������ʱ����s
end

%% ���㲻ͬʱ���������������L_coil���̹ܶ������ܳ���L_reel
L=10000; % �����͹��ܳ���m�����룩
L_wg=8; % ���ڵ�ע��ͷ�����������ܳ��ȣ�m�����룩
L_goose=3; % �������������ܳ��ȣ�m�����룩
D_goose=2; % �������ΰ뾶��m�����룩
H_goose=10; % �����������߶ȣ�m�����룩
L_gr=20; % ����������Ͳ�������ܳ��ȣ�m�����룩
theta_gr=acosd(H_goose/L_gr); % ����������Ͳ����������Ǧ���߼нǣ���
D_r_i=3; % ��Ͳ�ھ���m�����룩
D_r_o=5; % ��Ͳ�⾶��m�����룩
W_r=5; % ��Ͳ��ȣ�m�����룩
D_cable=0.005; % �����⾶��m

L_coil(1)=Lp; % ������������ȣ�m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % ��ʼʱ���̹ܶ������ܳ��ȣ�m
for t=2:1:Nt
    L_coil(t)=Lp; % ������������ȣ�m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % �̹ܶ������ܳ��ȣ�m
end
L_cable=L_coil; % �������볤�ȣ�m

%% �������뼰Ԥ����
D_ct_o_0=0.04445; % �����͹��⾶��m�����룩
L1=2000; % �����ڵ�һ�������ܳ��ȣ�m�����룩
D_ct_i_1=0.03709; % �����ڵ�һ���������ھ���m�����룩
L2=2000; % �����ڵڶ��������ܳ��ȣ�m�����룩
D_ct_i_2=0.03653; % �����ڵڶ����������ھ���m�����룩
L3=2000; % �����ڵ����������ܳ��ȣ�m�����룩
D_ct_i_3=0.03555; % �����ڵ������������ھ���m�����룩
L4=2000; % �����ڵ��Ķ������ܳ��ȣ�m�����룩
D_ct_i_4=0.03489; % �����ڵ��Ķ��������ھ���m�����룩
L5=2000; % �����ڵ���������ܳ��ȣ�m�����룩
D_ct_i_5=0.03409; % �����ڵ�����������ھ���m�����룩

L_t_1=4000; % �ϲ��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_1=0.09718; %0.068;% �ϲ��͹ܣ����׹ܻ����ۣ��ھ���m�����룩
L_t_2=2200; % �²��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_2=0.09718; %0.13970; % �²��͹ܣ����׹ܻ����ۣ��ھ���m�����룩

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % �͹ܣ����׹ܻ����ۣ��ھ���m
        else
            D_t_i(t,x)=D_t_i_2; % �͹ܣ����׹ܻ����ۣ��ھ���m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            D_ct_o(t,x)=D_ct_o_0; % �����͹��⾶��m
        else
            D_ct_o(t,x)=0; % �����͹��⾶��m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            if (L_coil(t)-Depth(x))<=L1 && (L_coil(t)-Depth(x))>=0
                D_ct_i(t,x)=D_ct_i_1; % �����͹��ھ���m
            elseif (L_coil(t)-Depth(x))<=L1+L2
                D_ct_i(t,x)=D_ct_i_2; % �����͹��ھ���m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                D_ct_i(t,x)=D_ct_i_3; % �����͹��ھ���m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                D_ct_i(t,x)=D_ct_i_4; % �����͹��ھ���m
            elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                D_ct_i(t,x)=D_ct_i_5; % �����͹��ھ���m
            end
        else
            D_ct_i(t,x)=0; % �����͹��ھ���m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % �����ļ������������ֵ
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % �����͹��ڽ������m^2
        else
            A_ct(t,x)=0; % �����͹��ڽ������m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % ����ˮ��ֱ����m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % ���ս������m^2
    end
end

%% �������뼰Ԥ����
h_t=25.4*10^(-6); % �͹ܣ����׹ܻ����ۣ����Դֲڶȣ�m�����룩
h_ct=25.4*10^(-6); % �����͹ܾ��Դֲڶȣ�m�����룩
h_a=(h_t+h_ct)/2; % ����ƽ�����Դֲڶȣ�m
epsilon_e=1*10^(-3); % �����������ޣ����룩
epsilon_t=1*10^3; % ���������������룩
g=9.81; % �������ٶȣ�m/s^2��Ĭ�ϣ�

T_0=20; % ��Һ�����¶ȣ��棨���룩
P_0=0.1*10^6; % ��Һ����ѹ����Pa�����룩
rho_l_0=1150; % T_0��P_0�»�Һ�ܶȣ�kg/m^3�����룩
mu_l_0=0.03; % T_0��P_0�»�Һճ�ȣ�Pa*s�����룩
Qv_l_0=0.1/60; % ��Һ���������m^3/s�����룩
Qm_l_0=Qv_l_0*rho_l_0; % ��Һ����������kg/s
rho_g_0=0.655; % ע�����ܶȣ�kg/m^3
mu_g_0=RheologyG(T_0,P_0); % ����ճ�ȣ�Pa*s
Qv_g_0=4/60; % ע�������������m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % ע��������������kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % ��ĭ����������kg/s

D_s=1*10^(-3); % ɰ��ֱ����m�����룩
rho_s=2000; % ɰ���ܶȣ�kg/m^3�����룩
H_s=L_s_b-Dsp_pen(1); % �ײ�ɰ���߶ȣ�m
PHI=0.6; % ɰ����ӯ�ȣ����룩
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % ����ɰ����������kg

D_nozzle=4/1000; % ����ֱ����m�����룩
N_nozzle=3; % ������������룩
C=0.95; % ��������ϵ����ȡ0.95�����룩

C0=1.2; % Ư��������ϵ����Ĭ�ϣ�

M_s(1)=0; % ���׽�ɰ����kg/s
for t=2:1:Nt
    M_s(t)=0; % ���׽�ɰ����kg/s
end

OutPressure=1*10^6; % ����ѹ����Pa�����룩

%% �¶�����
T_i=20; % ��ĭע���¶ȣ���
T_g=0.02; % �����ݶȣ���/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % ����������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������
    end
end
T_a=T_ct; % ������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������

%% ��1��ʱ��ڵ㣨��ʼʱ�̣���ز�����ֵ���㣨���գ�
for x=1:1:Nx
    P_a(1,x)=P_a_6(x); % ����ѹ����Pa
    rho_g_a(1,x)=rho_g_a_6(x); % �����ܶȣ�kg/m^3
    rho_l_a(1,x)=rho_l_a_6(x); % ����Һ���ܶȣ�kg/m^3
    rho_f_a(1,x)=rho_f_a_6(x); % ������ĭ�ܶȣ�kg/m^3
    mu_g_a(1,x)=mu_g_a_6(x); % ����ճ�ȣ�Pa*s
    mu_l_a(1,x)=mu_l_a_6(x); % ����Һ��ճ�ȣ�Pa*s
    mu_f_a(1,x)=mu_f_a_6(x); % ������ĭճ�ȣ�Pa*s
    mu_s(1,x)=mu_s_6(x); % ����ճ�ȣ�Pa*s
    alpha_g_a(1,x)=alpha_g_a_6(x); % ���ຬ��
    alpha_l_a(1,x)=alpha_l_a_6(x); % Һ�ຬ��
    alpha_f_a(1,x)=alpha_f_a_6(x); % ��ĭ����
    alpha_s(1,x)=alpha_s_6(1,x); % ���ຬ��
    Va_g_a(1,x)=Va_g_a_6(x); % ���������٣�m/s
    Va_l_a(1,x)=Va_l_a_6(x); % ����Һ�������٣�m/s
    Va_f_a(1,x)=Va_f_a_6(x); % ������ĭ������٣�m/s
    Va_s(1,x)=Va_s_6(x); % ���������٣�m/s
    V_g_a(1,x)=V_g_a_6(x); % �������٣�m/s
    V_l_a(1,x)=V_l_a_6(x); % ����Һ�����٣�m/s
    V_f_a(1,x)=V_f_a_6(x); % ������ĭ���٣�m/s
    V_s(1,x)=V_s_6(x); % �������٣�m/s
    Vsr(1,x)=Vsr_6(x); % ���໬���ٶȣ�m/s
    gamma_g_a(1,x)=gamma_g_a_6(x); % ��ĭ����
    gamma_l_a(1,x)=gamma_l_a_6(1,x); % Һ��������
    V_m(1,x)=V_m_6(x); % ���ջ�������٣�m/s
    rho_m(1,x)=rho_m_6(x); % ���ջ�����ܶȣ�kg/m^3
    mu_m(1,x)=mu_m_6(x); % ���ջ����ճ�ȣ�Pa*s
    Ff_a(1,x)=Ff_a_6(x); % ���յ�λ����Ħ��ѹ����Pa/m
    f_a(1,x)=f_a_6(x); % ����Ħ������
    Re_a(1,x)=Re_a_6(x); % ������ŵ��
    flow_pattern_a(1,x)=flow_pattern_a_6(x); % ����������̬
end

%% ��2��Nt��ʱ��ڵ���ز������㣨���գ�
for t=2:1:Nt
    P_a(t,Nx_Lp)=P_a(t-1,Nx_Lp);  % ���չܵ�ѹ������ֵ��Pa
    
    err_OutPressure=1; % ����ѹ��������
    COUNT_OutPressure=0; % ����ѹ������������ֵ
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % ����ѹ����������
        
        % �ܵף���Nx_Lp���ռ�ڵ㣩����ز�������
        rho_g_a(t,Nx_Lp)=DensityG(T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % ���������ܶȣ�kg/m^3
        rho_l_a(t,Nx_Lp)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % ���ջ�Һ�ܶȣ�kg/m^3
        mu_g_a(t,Nx_Lp)=RheologyG(T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % ��������ճ�ȣ�Pa*s
        mu_l_a(t,Nx_Lp)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Lp),P_a(t,Nx_Lp)); % ���ջ�Һճ�ȣ�Pa*s
        gamma_g_a(t,Nx_Lp)=(Qm_g_0/rho_g_a(t,Nx_Lp))/(Qm_g_0/rho_g_a(t,Nx_Lp)+Qm_l_0/rho_l_a(t,Nx_Lp)); % ��ĭ����
        gamma_l_a(t,Nx_Lp)=1-gamma_g_a(t,Nx_Lp); % Һ��������
        alpha_g_a(t,Nx_Lp)=alpha_f_a(t-1,Nx_Lp)*gamma_g_a(t,Nx_Lp); % �������庬��
        alpha_l_a(t,Nx_Lp)=alpha_f_a(t-1,Nx_Lp)*gamma_l_a(t,Nx_Lp); % ���ջ�Һ����
        rho_f_a(t,Nx_Lp)=rho_g_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp)+rho_l_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % ������ĭ�ܶȣ�kg/m^3
        mu_f_a(t,Nx_Lp)=mu_g_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp)+mu_l_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % ������ĭճ�ȣ�Pa*s
        mu_s(t,Nx_Lp)=mu_f_a(t,Nx_Lp); % ����ճ�ȣ�Pa*s
        Va_f_a(t,Nx_Lp)=Qm_f_0/(A_a(t,Nx_Lp)*rho_f_a(t,Nx_Lp)); % ������ĭ������٣�m/s
        Va_g_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp); % �������������٣�m/s
        Va_l_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp); % ���ջ�Һ������٣�m/s
        Va_s(t,Nx_Lp)=M_s(t)/(A_a(t,Nx_Lp)*rho_s); % ��м����ٶȣ�m/s
        Vsr(t,Nx_Lp)=12*(mu_f_a(t,Nx_Lp)/(rho_f_a(t,Nx_Lp)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Lp))/rho_f_a(t,Nx_Lp))*((rho_f_a(t,Nx_Lp)*D_s/mu_f_a(t,Nx_Lp))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
        alpha_s(t,Nx_Lp)=Va_s(t,Nx_Lp)/(C0*(Va_s(t,Nx_Lp)+Va_f_a(t,Nx_Lp))-Vsr(t,Nx_Lp));  % �����������
        V_s(t,Nx_Lp)=0;%Va_s(t,Nx_Lp)/alpha_s(t,Nx_Lp); % ��м�ٶȣ�m/s
        
        alpha_f_a(t,Nx_Lp)=1-alpha_s(t,Nx_Lp); % ������ĭ����
        alpha_g_a(t,Nx_Lp)=alpha_f_a(t,Nx_Lp)*gamma_g_a(t,Nx_Lp); % �������庬��
        alpha_l_a(t,Nx_Lp)=alpha_f_a(t,Nx_Lp)*gamma_l_a(t,Nx_Lp); % ���ջ�Һ����
        V_f_a(t,Nx_Lp)=Va_f_a(t,Nx_Lp)/alpha_f_a(t,Nx_Lp); % ������ĭ���٣�m/s
        V_g_a(t,Nx_Lp)=V_f_a(t,Nx_Lp); % �����������٣�m/s
        V_l_a(t,Nx_Lp)=V_f_a(t,Nx_Lp); % ���ջ�Һ���٣�m/s
        
        V_m(t,Nx_Lp)=Va_s(t,Nx_Lp)+Va_f_a(t,Nx_Lp); % ���ջ�����ٶȣ�m/s
        rho_m(t,Nx_Lp)=alpha_s(t,Nx_Lp)*rho_s+alpha_f_a(t,Nx_Lp)*rho_f_a(t,Nx_Lp); % ���ջ�����ܶȣ�kg/m^3
        mu_m(t,Nx_Lp)=alpha_s(t,Nx_Lp)*mu_s(t,Nx_Lp)+alpha_f_a(t,Nx_Lp)*mu_f_a(t,Nx_Lp); % ���ջ����ճ�ȣ�Pa*s
        [Ff_a(t,Nx_Lp),f_a(t,Nx_Lp),Re_a(t,Nx_Lp),flow_pattern_a(t,Nx_Lp)]=Friction_annulus(rho_m(t,Nx_Lp),V_m(t,Nx_Lp),mu_m(t,Nx_Lp),D_h(t,Nx_Lp),h_a,rho_f_a(t,Nx_Lp),V_f_a(t,Nx_Lp)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
        
        % ��Nx_Lp-1��1���ռ�ڵ㴦��ز�������
        for x=Nx_Lp-1:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % ����ѹ������ֵ��Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % ���ջ�Һ����
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            
            err_NodePressure=1; % ����ѹ��������
            COUNT_NodePressure=0; % ����ѹ������������ֵ
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % ����ѹ����������
                
                % ���߷��������������
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % ���������������ֵ1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % ���������������ֵ2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % ������������������
                COUNT_NodeEg=0; % ��������
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % �����������Ϊalpha_s_ass1ʱ������ֵ��
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % �����ٶȣ�m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % ���������٣�m/s
                    alpha_f_ass1=1-alpha_s_ass1; % ��ĭ�������
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % ��ĭ�ٶȣ�m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % ��ĭ������٣�m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % ���ջ�����ٶȣ�m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % ���������������ֵ
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % ����ĺ��������Ľ������ʵ�����������
                    
                    % �����������Ϊalpha_s_ass2ʱ������ֵ��
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % ���������غ㷽����ɢ��ʽ�м�ֵ����
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % �����ٶȣ�m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % ���������٣�m/s
                    alpha_f_ass2=1-alpha_s_ass2; % ��ĭ�������
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % ��ĭ�ٶȣ�m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % ��ĭ������٣�m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % ���ջ�����ٶȣ�m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % ���������������ֵ
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % ����ĺ��������Ľ������ʵ�����������
                    
                    % ���߷���������������
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % �µĹ��������������ֵalpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % ������������������
                    alpha_s_ass1=alpha_s_ass2; % �µĹ��������������ֵ1
                    alpha_s_ass2=alpha_s_ass3; % �µĹ��������������ֵ2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % ���������õ�����ʵ�����������ֵ����alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % �������������С��һ��ֵʱ����Ϊ���������Ϊ0�����ڷ�ֹ���ֺ����ļ������
                    alpha_s(t,x)=0; % �����������
                    V_s(t,x)=0; % �����ٶȣ�m/s
                    Va_s(t,x)=0; % �������ٶȣ�m/s
                    mu_s(t,x)=0; % ����ճ�ȣ�Pa*s
                    Vsr(t,x)=0; % ��м����ĩ�٣�m/s
                    
                    alpha_f_a(t,x)=1; % ��ĭ�������                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                else
                    V_s(t,x)=V_s_ass2; % �����ٶȣ�m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % ���������٣�m/s
                    Vsr(t,x)=Vsr_ass2; % ��м����ĩ�٣�m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % ��ĭ�ܶȣ�kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % ��ĭ�������
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    V_f_a(t,x)=V_f_ass2; % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % ���ջ�����ٶȣ�m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % ����ѹ������ֵ��Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % ����ѹ��������
                P_a(t,x)=P_new; % �µĻ���ѹ������ֵ��Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % ����ѹ��������        
        if (P_a(t,1)-OutPressure)>0 % ���ݳ���ѹ�������������Ի��չܵ�ѹ������ֵ���е���
            P_a(t,Nx_Lp)=P_a(t,Nx_Lp)-(P_a(t,1)-OutPressure)/2; % �µĻ��չܵ�ѹ������ֵ��Pa
        else
            P_a(t,Nx_Lp)=P_a(t,Nx_Lp)-(P_a(t,1)-OutPressure)/2*0.3; % �µĻ��չܵ�ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Lp+1��Nx_Dsp_pen(1)���ռ�ڵ㴦��ز�������
    for x=Nx_Lp+1:1:Nx_Dsp_pen(1)
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=0; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Dsp_pen(1)+1��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Dsp_pen(1)+1:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=0; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
end

%% ��1��Ntʱ��ڵ���ͷѹ������
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Lp); % ������ĭ���������m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Lp)); % �����ɰ����ѹ����Pa����������������V_nozzle��m/s��
end

%% ��1��Ntʱ��ڵ���ز������㣨�������ڣ�
for t=1:1:Nt
    % �ܵף���Nx_Lp���ռ�ڵ㣩����ز�������
    P_ct(t,Nx_Lp)=P_a(t,Nx_Lp)+delta_P_SWT(t); % ����ѹ����Pa
    rho_g_ct(t,Nx_Lp)=DensityG(T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % ���������ܶȣ�kg/m^3
    rho_l_ct(t,Nx_Lp)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % ���ڻ�Һ�ܶȣ�kg/m^3
    mu_g_ct(t,Nx_Lp)=RheologyG(T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % ��������ճ�ȣ�Pa*s
    mu_l_ct(t,Nx_Lp)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Lp),P_ct(t,Nx_Lp)); % ���ڻ�Һճ�ȣ�Pa*s
    alpha_f_ct(t,Nx_Lp)=1; % ������ĭ����
    gamma_g_ct(t,Nx_Lp)=(Qm_g_0/rho_g_ct(t,Nx_Lp))/(Qm_g_0/rho_g_ct(t,Nx_Lp)+Qm_l_0/rho_l_ct(t,Nx_Lp)); % ������ĭ����
    gamma_l_ct(t,Nx_Lp)=1-gamma_g_ct(t,Nx_Lp); % ����Һ��������
    alpha_g_ct(t,Nx_Lp)=alpha_f_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp); % �������庬��
    alpha_l_ct(t,Nx_Lp)=alpha_f_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % ���ڻ�Һ����
    rho_f_ct(t,Nx_Lp)=rho_g_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp)+rho_l_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % ������ĭ�ܶȣ�kg/m^3
    mu_f_ct(t,Nx_Lp)=mu_g_ct(t,Nx_Lp)*gamma_g_ct(t,Nx_Lp)+mu_l_ct(t,Nx_Lp)*gamma_l_ct(t,Nx_Lp); % ������ĭճ�ȣ�Pa*s
    V_f_ct(t,Nx_Lp)=Qm_f_0/(rho_f_ct(t,Nx_Lp)*A_ct(t,Nx_Lp)); % �����������٣�m/s
    [Ff_ct(t,Nx_Lp),f_ct(t,Nx_Lp),Re_ct(t,Nx_Lp),flow_pattern_ct(t,Nx_Lp)]=Friction_coiledtubing(rho_f_ct(t,Nx_Lp),V_f_ct(t,Nx_Lp),mu_f_ct(t,Nx_Lp),D_ct_i(t,Nx_Lp),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    % ��Nx_Lp-1��1���ռ�ڵ㴦��ز�������
    for x=Nx_Lp-1:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % ����ѹ������ֵ��Pa
        
        err_DriPipePressure=1; % ����ѹ��������
        COUNT_DriPipePressure=0; % ����ѹ������������ֵ
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % ����ѹ����������
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һ�ܶȣ�kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һճ�ȣ�Pa*s
            alpha_f_ct(t,x)=1; % ������ĭ����
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % ������ĭ����
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % �������庬��
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % ���ڻ�Һ����
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭճ�ȣ�Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % ������ĭ���٣�m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % ����ѹ����Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % �������ѹ������ֵ�����ֵ֮���������
            P_ct_ass(t,x)=P_ct(t,x); % �µĹ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Lp+1��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Lp+1:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % �������庬��
        alpha_l_ct(t,x)=alpha_l_a(t,x); % ���ڻ�Һ����
        gamma_g_ct(t,x)=gamma_g_a(t,x); % ������ĭ����
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
        rho_f_ct(t,x)=rho_m(t,x); % ���������ܶȣ�kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % ��������ճ�ȣ�Pa*s
        V_f_ct(t,t)=V_m(t,x); % �����������٣�sm/s
        Ff_ct(t,x)=Ff_a(t,x); % �������嵥λ����Ħ��ѹ����Pa/m��
        f_ct(t,x)=f_a(t,x); % �������巶��Ħ������
        Re_ct(t,x)=Re_a(t,x); % ����������ŵ����������̬
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % ����������̬
        P_ct(t,x)=P_a(t,x); % ����ѹ����Pa
    end
end

%% ��1��Ntʱ��ڵ��̹ܶγ���ѹ������ѹ���㣨�������ڣ�
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % ��������������٣�m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % ����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % �̹ܶ����嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % �����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % ���ڵ�ע��ͷ������Ħ��ѹ����Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % ��������Ħ��ѹ����Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % ����������Ͳ��Ħ��ѹ����Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % �̹ܶ�Ħ��ѹ����Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % �̹ܶγ���ѹ����Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % ��ѹ��Pa
end

%% ���㾮���ۻ���ɰ��
M_w_tem(1)=0; % ����˲ʱ��ɰ����ֵ��kg
M_w_tot(1)=0; % �����ۻ���ɰ����ֵ��kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮��˲ʱ��ɰ����kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮���ۻ���ɰ����kg
end

%% ����ĩ̬��Ͳ���ɰŨ��
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(1)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% ���㻷����ĭ����ƽ��ֵV_l_mre
V_f_mre=0; % ���շ���ƽ��ֵ��m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % ���շ���ƽ��ֵ��m/s
    end
end

%% ����������ĩ��ƽ��ֵVsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % ��Vsr��������ȡֵ����ĳ���ĩ�٣�m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % �������ĩ��ƽ��ֵ��m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % �������ĩ��ƽ��ֵ��m/s
end

%% �ж��Ƿ���Ч��ɰ�����շ���ƽ��ֵ����2��ɰ������ĩ��ƽ��ֵ�����Ƿ���ɳ�ɰ
fprintf("����ѭ������2��\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % �����д���˵����ǰ�����¡���ɳ�ɰ��
    else
        fprintf("Sand Cleanout UnFinished!\n"); % �����д���˵����ǰ�����¡�δ��ɳ�ɰ��
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
end

%% ��ͲECD����
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % ����ECD��kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % ����ECD��kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
end

%% ����ѭ������ĩ״̬����
alpha_f_a_7=alpha_f_a(Nt,:); % ����ѭ������ĩ״̬��ĭ����
alpha_g_a_7=alpha_g_a(Nt,:); % ����ѭ������ĩ״̬���ຬ��
alpha_l_a_7=alpha_l_a(Nt,:); % ����ѭ������ĩ״̬Һ�ຬ��
alpha_s_7=alpha_s(Nt,:); % ����ѭ������ĩ״̬���ຬ��
f_a_7=f_a(Nt,:); % ����ѭ������ĩ״̬����Ħ������
Ff_a_7=Ff_a(Nt,:); % ����ѭ������ĩ״̬���յ�λ����Ħ��ѹ����Pa/m
flow_pattern_a_7=flow_pattern_a(Nt,:); % ����ѭ������ĩ״̬����������̬
gamma_g_a_7=gamma_g_a(Nt,:); % ����ѭ������ĩ״̬��ĭ����
gamma_l_a_7=gamma_l_a(Nt,:); % ����ѭ������ĩ״̬Һ��������
mu_f_a_7=mu_f_a(Nt,:); % ����ѭ������ĩ״̬������ĭճ�ȣ�Pa*s
mu_g_a_7=mu_g_a(Nt,:); % ����ѭ������ĩ״̬����ճ�ȣ�Pa*s
mu_l_a_7=mu_l_a(Nt,:); % ����ѭ������ĩ״̬����Һ��ճ�ȣ�Pa*s
mu_m_7=mu_m(Nt,:); % ����ѭ������ĩ״̬���ջ����ճ�ȣ�Pa*s
mu_s_7=mu_s(Nt,:); % ����ѭ������ĩ״̬����ճ�ȣ�Pa*s
P_a_7=P_a(Nt,:); % ����ѭ������ĩ״̬����ѹ����Pa
Re_a_7=Re_a(Nt,:); % ����ѭ������ĩ״̬������ŵ��
rho_f_a_7=rho_f_a(Nt,:); % ����ѭ������ĩ״̬������ĭ�ܶȣ�kg/m^3
rho_g_a_7=rho_g_a(Nt,:); % ����ѭ������ĩ״̬���������ܶȣ�kg/m^3
rho_l_a_7=rho_l_a(Nt,:); % ����ѭ������ĩ״̬����Һ���ܶȣ�kg/m^3
rho_m_7=rho_m(Nt,:); % ����ѭ������ĩ״̬���ջ�����ܶȣ�kg/m^3
V_f_a_7=V_f_a(Nt,:); % ����ѭ������ĩ״̬������ĭ���٣�m/s
V_g_a_7=V_g_a(Nt,:); % ����ѭ������ĩ״̬�������٣�m/s
V_l_a_7=V_l_a(Nt,:); % ����ѭ������ĩ״̬����Һ�����٣�m/s
V_m_7=V_m(Nt,:); % ����ѭ������ĩ״̬���ջ�������٣�m/s
V_s_7=V_s(Nt,:); % ����ѭ������ĩ״̬�������٣�m/s
Va_f_a_7=Va_f_a(Nt,:); % ����ѭ������ĩ״̬������ĭ������٣�m/s
Va_g_a_7=Va_g_a(Nt,:); % ����ѭ������ĩ״̬���������٣�m/s
Va_l_a_7=Va_l_a(Nt,:); % ����ѭ������ĩ״̬����Һ�������٣�m/s
Va_s_7=Va_s(Nt,:); % ����ѭ������ĩ״̬���������٣�m/s
Vsr_7=Vsr(Nt,:); % ����ѭ������ĩ״̬���໬���ٶȣ�m/s

%% ���ݴ洢������ѭ����
ANS_Nt_7=Nt; % ʱ��ڵ���
ANS_Time_7=Time(1:Nt); % ����ѭ��ʱ�䣬s
ANS_alpha_g_a_7=alpha_g_a(1:Nt,:); % ���������������
ANS_alpha_l_a_7=alpha_l_a(1:Nt,:); % ����Һ���������
ANS_alpha_s_7=alpha_s(1:Nt,:); % ���ຬ��
ANS_P_ct_7=P_ct(1:Nt,:); % ����ѹ����Pa��
ANS_delta_P_SWT_7=delta_P_SWT(1:Nt); % ��ɰ����ѹ����Pa��
ANS_P_a_7=P_a(1:Nt,:); % ����ѹ����Pa��
ANS_T_a_7=T_a(1:Nt,:); % �����¶ȣ��棩
ANS_T_ct_7=T_ct(1:Nt,:); % �����¶ȣ��棩
ANS_P_coil_7=P_coil(1:Nt); % �̹ܶγ���ѹ����Pa��
ANS_P_pump_7=P_pump(1:Nt); % ��ѹ��Pa��
ANS_M_w_tem_7=M_w_tem(1:Nt); % ����˲ʱ��ɰ����kg��
ANS_M_w_tot_7=M_w_tot(1:Nt); % �����ۻ���ɰ����kg��
ANS_Va_s_7=Va_s(1:Nt,:); % ��м����ٶȣ�m/s��
ANS_Va_f_a_7=Va_f_a(1:Nt,:); % ������ĭ������٣�m/s��
ANS_V_s_7=V_s(1:Nt,:); % ��м�����ٶȣ�m/s��
ANS_V_f_a_7=V_f_a(1:Nt,:); % ������ĭ���٣�m/s��
ANS_alpha_g_ct_7=alpha_g_ct(1:Nt,:); % ���������������
ANS_alpha_l_ct_7=alpha_l_ct(1:Nt,:); % ����Һ���������
ANS_gamma_g_ct_7=gamma_g_ct(1:Nt,:); % ������ĭ����
ANS_gamma_g_a_7=gamma_g_a(1:Nt,:); % ������ĭ����
ANS_ECD_a_7=ECD_a(1:Nt,:); % ����ECD��kg/m^3��
ANS_L_coil_7=L_coil(1:Nt); % ���������m��



%% ������������̣����ʵ�λ�ƣ�
%% ������ֵ����
dt=zeros(); % ʱ�䲽����s
dx=zeros(); % �ռ䲽����m
Time=zeros(); % ������������ʱ����s
L_coil=zeros(); % ���������m
L_reel=zeros(); % �̹ܶγ��ȣ�m
D_t_i=zeros(); % �͹��ھ���m
D_ct_o=zeros(); % �����͹��⾶��m
D_ct_i=zeros(); % �����͹��ھ���m
M_s=zeros(); % ���׽�ɰ����kg/s
P_a=zeros(); % ����ѹ����Pa
rho_g_a=zeros(); % �����ܶȣ�kg/m^3
rho_l_a=zeros(); % ����Һ���ܶȣ�kg/m^3
rho_f_a=zeros(); % ������ĭ�ܶȣ�kg/m^3
mu_g_a=zeros(); % ����ճ�ȣ�Pa*s
mu_l_a=zeros(); % ����Һ��ճ�ȣ�Pa*s
mu_f_a=zeros(); % ������ĭճ�ȣ�Pa*s
mu_s=zeros(); % ����ճ�ȣ�Pa*s
alpha_g_a=zeros(); % ���ຬ��
alpha_l_a=zeros(); % Һ�ຬ��
alpha_f_a=zeros(); % ��ĭ����
alpha_s=zeros(); % ���ຬ��
Va_g_a=zeros(); % ���������٣�m/s
Va_l_a=zeros(); % ����Һ�������٣�m/s
Va_f_a=zeros(); % ������ĭ������٣�m/s
Va_s=zeros(); % ���������٣�m/s
V_g_a=zeros(); % �������٣�m/s
V_l_a=zeros(); % ����Һ�����٣�m/s
V_f_a=zeros(); % ������ĭ���٣�m/s
V_s=zeros(); % �������٣�m/s
Vsr=zeros(); % ���໬���ٶȣ�m/s
gamma_g_a=zeros(); % ��ĭ����
gamma_l_a=zeros(); % Һ��������
V_m=zeros(); % ���ջ�������٣�m/s
rho_m=zeros(); % ���ջ�����ܶȣ�kg/m^3
mu_m=zeros(); % ���ջ����ճ�ȣ�Pa*s
Ff_a=zeros(); % ���յ�λ����Ħ��ѹ����Pa/m
f_a=zeros(); % ����Ħ������
Re_a=zeros(); % ������ŵ��
flow_pattern_a=zeros(); % ����������̬

%% ����ռ䲽��dx����Ӧʱ�䲽��dt
V8=Lp/(1.5*3600); % ����������ٶȣ�m/s�����룩
t_8=Lp/V8; % �����������ʱ����s

Nt=Nx_Lp; % ʱ��ڵ���
nt=Nt-1; % ʱ��������

for x=1:1:Nx-1
    dx(x)=Depth(x+1)-Depth(x); % �ռ䲽����m
end

Time(1)=0; % �����������ʱ����ֵ��s
for t=1:1:Nt-1
    dt(t)=dx(Nx_Lp-t)/V8; % ÿ���һ���ռ䲽������ʱ�䣬s
    Time(t+1)=Time(t)+dt(t); % �������������(t+1)���ռ�ڵ�����������ʱ����s
end

%% ���㲻ͬʱ��������������L_coil���̹ܶ������ܳ���L_reel
L=10000; % �����͹��ܳ���m�����룩
L_wg=8; % ���ڵ�ע��ͷ�����������ܳ��ȣ�m�����룩
L_goose=3; % �������������ܳ��ȣ�m�����룩
D_goose=2; % �������ΰ뾶��m�����룩
H_goose=10; % �����������߶ȣ�m�����룩
L_gr=20; % ����������Ͳ�������ܳ��ȣ�m�����룩
theta_gr=acosd(H_goose/L_gr); % ����������Ͳ����������Ǧ���߼нǣ���
D_r_i=3; % ��Ͳ�ھ���m�����룩
D_r_o=5; % ��Ͳ�⾶��m�����룩
W_r=5; % ��Ͳ��ȣ�m�����룩
D_cable=0.005; % �����⾶��m

L_coil(1)=Lp; % ��ʼʱ�������ܵײ���ȣ�m
L_reel(1)=L-L_coil(1)-L_wg-L_goose-L_gr; % ��ʼʱ���̹ܶ������ܳ��ȣ�m
for t=2:1:Nt
    L_coil(t)=L_coil(t-1)-dx(Nx_Lp-t+1); % �����������ȣ�m
    L_reel(t)=L-L_coil(t)-L_wg-L_goose-L_gr; % �̹ܶ������ܳ��ȣ�m
end
L_cable=L_coil; % �������볤�ȣ�m

%% �������뼰Ԥ����
D_ct_o_0=0.04445; % �����͹��⾶��m�����룩
L1=2000; % �����ڵ�һ�������ܳ��ȣ�m�����룩
D_ct_i_1=0.03709; % �����ڵ�һ���������ھ���m�����룩
L2=2000; % �����ڵڶ��������ܳ��ȣ�m�����룩
D_ct_i_2=0.03653; % �����ڵڶ����������ھ���m�����룩
L3=2000; % �����ڵ����������ܳ��ȣ�m�����룩
D_ct_i_3=0.03555; % �����ڵ������������ھ���m�����룩
L4=2000; % �����ڵ��Ķ������ܳ��ȣ�m�����룩
D_ct_i_4=0.03489; % �����ڵ��Ķ��������ھ���m�����룩
L5=2000; % �����ڵ���������ܳ��ȣ�m�����룩
D_ct_i_5=0.03409; % �����ڵ�����������ھ���m�����룩

L_t_1=4000; % �ϲ��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_1=0.09718; %0.068;%0.09718; % �ϲ��͹ܣ����׹ܻ����ۣ��ھ���m�����룩
L_t_2=2200; % �²��͹ܣ����׹ܻ����ۣ����ȣ�m�����룩
D_t_i_2=0.09718; %0.13970; % �²��͹ܣ����׹ܻ����ۣ��ھ���m�����룩

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_t_1
            D_t_i(t,x)=D_t_i_1; % �͹ܣ����׹ܻ����ۣ��ھ���m
        else
            D_t_i(t,x)=D_t_i_2; % �͹ܣ����׹ܻ����ۣ��ھ���m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            if Depth(x)<=L_coil(t)
                D_ct_o(t,x)=D_ct_o_0; % �����͹��⾶��m
            else
                D_ct_o(t,x)=0; % �����͹��⾶��m
            end
        else
            D_ct_o(t,x)=0; % �����͹��⾶��m
        end
    end
end
for t=1:1:Nt
    for x=1:1:Nx
        if x<=Nx_Lp
            if Depth(x)<=L_coil(t) % �����ļ�����������ھ�
                if (L_coil(t)-Depth(x))<=L1 && (L_coil(t)-Depth(x))>=0
                    D_ct_i(t,x)=D_ct_i_1; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2
                    D_ct_i(t,x)=D_ct_i_2; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3
                    D_ct_i(t,x)=D_ct_i_3; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4
                    D_ct_i(t,x)=D_ct_i_4; % �����͹��ھ���m
                elseif (L_coil(t)-Depth(x))<=L1+L2+L3+L4+L5
                    D_ct_i(t,x)=D_ct_i_5; % �����͹��ھ���m
                end
            else
                D_ct_i(t,x)=0; % �����͹��ھ���m
            end
        else
            D_ct_i(t,x)=0; % �����͹��ھ���m
        end
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if Depth(x)<=L_coil(t) % �����ļ������������ֵ
            A_ct(t,x)=1/4*pi*(D_ct_i(t,x)^2-D_cable^2); % �����͹��ڽ������m^2
        else
            A_ct(t,x)=0; % �����͹��ڽ������m^2
        end
        D_h(t,x)=D_t_i(t,x)-D_ct_o(t,x); % ����ˮ��ֱ����m
        A_a(t,x)=1/4*pi*(D_t_i(t,x)^2-D_ct_o(t,x)^2); % ���ս������m^2
    end
end

%% �������뼰Ԥ����
h_t=25.4*10^(-6); % �͹ܣ����׹ܻ����ۣ����Դֲڶȣ�m�����룩
h_ct=25.4*10^(-6); % �����͹ܾ��Դֲڶȣ�m�����룩
h_a=(h_t+h_ct)/2; % ����ƽ�����Դֲڶȣ�m
epsilon_e=1*10^(-3); % �����������ޣ����룩
epsilon_t=1*10^3; % ���������������룩
g=9.81; % �������ٶȣ�m/s^2��Ĭ�ϣ�

T_0=20; % ��Һ�����¶ȣ��棨���룩
P_0=0.1*10^6; % ��Һ����ѹ����Pa�����룩
rho_l_0=1150; % T_0��P_0�»�Һ�ܶȣ�kg/m^3�����룩
mu_l_0=0.03; % T_0��P_0�»�Һճ�ȣ�Pa*s�����룩
Qv_l_0=0.1/60; % ��Һ���������m^3/s�����룩
Qm_l_0=Qv_l_0*rho_l_0; % ��Һ����������kg/s
rho_g_0=0.655; % ע�����ܶȣ�kg/m^3
mu_g_0=RheologyG(T_0,P_0); % ����ճ�ȣ�Pa*s
Qv_g_0=4/60; % ע�������������m^3/s
Qm_g_0=Qv_g_0*rho_g_0;  % ע��������������kg/s
Qm_f_0=Qm_g_0+Qm_l_0; % ��ĭ����������kg/s

D_s=1*10^(-3); % ɰ��ֱ����m�����룩
rho_s=2000; % ɰ���ܶȣ�kg/m^3�����룩
H_s=L_s_b-L_s_t; % �ײ�ɰ���߶ȣ�m
PHI=0.6; % ɰ����ӯ�ȣ����룩
M_s_total=PHI*rho_s*1/4*pi*D_t_i(1,Nx)^2*H_s; % ����ɰ����������kg

D_nozzle=4/1000; % ����ֱ����m�����룩
N_nozzle=3; % ������������룩
C=0.95; % ��������ϵ����ȡ0.95�����룩

C0=1.2; % Ư��������ϵ����Ĭ�ϣ�

M_s(1)=0; % ���׽�ɰ����kg/s
for t=2:1:Nt
    M_s(t)=0; % ���׽�ɰ����kg/s
end

OutPressure=1*10^6; % ����ѹ����Pa�����룩

%% �¶�����
T_i=20; % ��ĭע���¶ȣ���
T_g=0.02; % �����ݶȣ���/m

for t=1:1:Nt
    T_ct(t,1)=T_i;
    for x=2:1:Nx
        T_ct(t,x)=T_ct(t,x-1)+T_g*dx(x-1); % ����������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������
    end
end
T_a=T_ct; % ������ĭ�¶ȣ�����ÿ��ʱ�̶�һ��������

%% ��1��ʱ��ڵ㣨��ʼʱ�̣���ز�����ֵ���㣨���գ�
for x=1:1:Nx
    P_a(1,x)=P_a_7(x); % ����ѹ����Pa
    rho_g_a(1,x)=rho_g_a_7(x); % �����ܶȣ�kg/m^3
    rho_l_a(1,x)=rho_l_a_7(x); % ����Һ���ܶȣ�kg/m^3
    rho_f_a(1,x)=rho_f_a_7(x); % ������ĭ�ܶȣ�kg/m^3
    mu_g_a(1,x)=mu_g_a_7(x); % ����ճ�ȣ�Pa*s
    mu_l_a(1,x)=mu_l_a_7(x); % ����Һ��ճ�ȣ�Pa*s
    mu_f_a(1,x)=mu_f_a_7(x); % ������ĭճ�ȣ�Pa*s
    mu_s(1,x)=mu_s_7(x); % ����ճ�ȣ�Pa*s
    alpha_g_a(1,x)=alpha_g_a_7(x); % ���ຬ��
    alpha_l_a(1,x)=alpha_l_a_7(x); % Һ�ຬ��
    alpha_f_a(1,x)=alpha_f_a_7(x); % ��ĭ����
    alpha_s(1,x)=alpha_s_7(1,x); % ���ຬ��
    Va_g_a(1,x)=Va_g_a_7(x); % ���������٣�m/s
    Va_l_a(1,x)=Va_l_a_7(x); % ����Һ�������٣�m/s
    Va_f_a(1,x)=Va_f_a_7(x); % ������ĭ������٣�m/s
    Va_s(1,x)=Va_s_7(x); % ���������٣�m/s
    V_g_a(1,x)=V_g_a_7(x); % �������٣�m/s
    V_l_a(1,x)=V_l_a_7(x); % ����Һ�����٣�m/s
    V_f_a(1,x)=V_f_a_7(x); % ������ĭ���٣�m/s
    V_s(1,x)=V_s_7(x); % �������٣�m/s
    Vsr(1,x)=Vsr_7(x); % ���໬���ٶȣ�m/s
    gamma_g_a(1,x)=gamma_g_a_7(x); % ��ĭ����
    gamma_l_a(1,x)=gamma_l_a_7(1,x); % Һ��������
    V_m(1,x)=V_m_7(x); % ���ջ�������٣�m/s
    rho_m(1,x)=rho_m_7(x); % ���ջ�����ܶȣ�kg/m^3
    mu_m(1,x)=mu_m_7(x); % ���ջ����ճ�ȣ�Pa*s
    Ff_a(1,x)=Ff_a_7(x); % ���յ�λ����Ħ��ѹ����Pa/m
    f_a(1,x)=f_a_7(x); % ����Ħ������
    Re_a(1,x)=Re_a_7(x); % ������ŵ��
    flow_pattern_a(1,x)=flow_pattern_a_7(x); % ����������̬
end

%% ��2��Nt��ʱ��ڵ���ز������㣨���գ�
for t=2:1:Nt
    P_a(t,Nx_Lp-t+1)=P_a(t-1,Nx_Lp-t+1);  % ���չܵ�ѹ������ֵ��Pa
    
    err_OutPressure=1; % ����ѹ��������
    COUNT_OutPressure=0; % ����ѹ������������ֵ
    while abs(err_OutPressure)>epsilon_e && COUNT_OutPressure<epsilon_t
        COUNT_OutPressure=COUNT_OutPressure+1;  % ����ѹ����������
        
        % �ܵף���Nx_Lp-t+1���ռ�ڵ㣩����ز�������
        rho_g_a(t,Nx_Lp-t+1)=DensityG(T_a(t,Nx_Lp-t+1),P_a(t,Nx_Lp-t+1)); % ���������ܶȣ�kg/m^3
        rho_l_a(t,Nx_Lp-t+1)=DensityL(rho_l_0,T_0,P_0,T_a(t,Nx_Lp-t+1),P_a(t,Nx_Lp-t+1)); % ���ջ�Һ�ܶȣ�kg/m^3
        mu_g_a(t,Nx_Lp-t+1)=RheologyG(T_a(t,Nx_Lp-t+1),P_a(t,Nx_Lp-t+1)); % ��������ճ�ȣ�Pa*s
        mu_l_a(t,Nx_Lp-t+1)=RheologyL(mu_l_0,T_0,P_0,T_a(t,Nx_Lp-t+1),P_a(t,Nx_Lp-t+1)); % ���ջ�Һճ�ȣ�Pa*s
        gamma_g_a(t,Nx_Lp-t+1)=(Qm_g_0/rho_g_a(t,Nx_Lp-t+1))/(Qm_g_0/rho_g_a(t,Nx_Lp-t+1)+Qm_l_0/rho_l_a(t,Nx_Lp-t+1)); % ��ĭ����
        gamma_l_a(t,Nx_Lp-t+1)=1-gamma_g_a(t,Nx_Lp-t+1); % Һ��������
        alpha_g_a(t,Nx_Lp-t+1)=alpha_f_a(t-1,Nx_Lp-t+1)*gamma_g_a(t,Nx_Lp-t+1); % �������庬��
        alpha_l_a(t,Nx_Lp-t+1)=alpha_f_a(t-1,Nx_Lp-t+1)*gamma_l_a(t,Nx_Lp-t+1); % ���ջ�Һ����
        rho_f_a(t,Nx_Lp-t+1)=rho_g_a(t,Nx_Lp-t+1)*gamma_g_a(t,Nx_Lp-t+1)+rho_l_a(t,Nx_Lp-t+1)*gamma_l_a(t,Nx_Lp-t+1); % ������ĭ�ܶȣ�kg/m^3
        mu_f_a(t,Nx_Lp-t+1)=mu_g_a(t,Nx_Lp-t+1)*gamma_g_a(t,Nx_Lp-t+1)+mu_l_a(t,Nx_Lp-t+1)*gamma_l_a(t,Nx_Lp-t+1); % ������ĭճ�ȣ�Pa*s
        mu_s(t,Nx_Lp-t+1)=mu_f_a(t,Nx_Lp-t+1); % ����ճ�ȣ�Pa*s
        Va_f_a(t,Nx_Lp-t+1)=Qm_f_0/(A_a(t,Nx_Lp-t+1)*rho_f_a(t,Nx_Lp-t+1)); % ������ĭ������٣�m/s
        Va_g_a(t,Nx_Lp-t+1)=Va_f_a(t,Nx_Lp-t+1); % �������������٣�m/s
        Va_l_a(t,Nx_Lp-t+1)=Va_f_a(t,Nx_Lp-t+1); % ���ջ�Һ������٣�m/s
        Va_s(t,Nx_Lp-t+1)=M_s(t)/(A_a(t,Nx_Lp-t+1)*rho_s); % ��м����ٶȣ�m/s
        Vsr(t,Nx_Lp-t+1)=12*(mu_f_a(t,Nx_Lp-t+1)/(rho_f_a(t,Nx_Lp-t+1)*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_a(t,Nx_Lp-t+1))/rho_f_a(t,Nx_Lp-t+1))*((rho_f_a(t,Nx_Lp-t+1)*D_s/mu_f_a(t,Nx_Lp-t+1))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
        alpha_s(t,Nx_Lp-t+1)=Va_s(t,Nx_Lp-t+1)/(C0*(Va_s(t,Nx_Lp-t+1)+Va_f_a(t,Nx_Lp-t+1))-Vsr(t,Nx_Lp-t+1));  % �����������
        V_s(t,Nx_Lp-t+1)=0;%Va_s(t,Nx_Lp-t+1)/alpha_s(t,Nx_Lp-t+1); % ��м�ٶȣ�m/s
        
        alpha_f_a(t,Nx_Lp-t+1)=1-alpha_s(t,Nx_Lp-t+1); % ������ĭ����
        alpha_g_a(t,Nx_Lp-t+1)=alpha_f_a(t,Nx_Lp-t+1)*gamma_g_a(t,Nx_Lp-t+1); % �������庬��
        alpha_l_a(t,Nx_Lp-t+1)=alpha_f_a(t,Nx_Lp-t+1)*gamma_l_a(t,Nx_Lp-t+1); % ���ջ�Һ����
        V_f_a(t,Nx_Lp-t+1)=Va_f_a(t,Nx_Lp-t+1)/alpha_f_a(t,Nx_Lp-t+1); % ������ĭ���٣�m/s
        V_g_a(t,Nx_Lp-t+1)=V_f_a(t,Nx_Lp-t+1); % �����������٣�m/s
        V_l_a(t,Nx_Lp-t+1)=V_f_a(t,Nx_Lp-t+1); % ���ջ�Һ���٣�m/s
        
        V_m(t,Nx_Lp-t+1)=Va_s(t,Nx_Lp-t+1)+Va_f_a(t,Nx_Lp-t+1); % ���ջ�����ٶȣ�m/s
        rho_m(t,Nx_Lp-t+1)=alpha_s(t,Nx_Lp-t+1)*rho_s+alpha_f_a(t,Nx_Lp-t+1)*rho_f_a(t,Nx_Lp-t+1); % ���ջ�����ܶȣ�kg/m^3
        mu_m(t,Nx_Lp-t+1)=alpha_s(t,Nx_Lp-t+1)*mu_s(t,Nx_Lp-t+1)+alpha_f_a(t,Nx_Lp-t+1)*mu_f_a(t,Nx_Lp-t+1); % ���ջ����ճ�ȣ�Pa*s
        [Ff_a(t,Nx_Lp-t+1),f_a(t,Nx_Lp-t+1),Re_a(t,Nx_Lp-t+1),flow_pattern_a(t,Nx_Lp-t+1)]=Friction_annulus(rho_m(t,Nx_Lp-t+1),V_m(t,Nx_Lp-t+1),mu_m(t,Nx_Lp-t+1),D_h(t,Nx_Lp-t+1),h_a,rho_f_a(t,Nx_Lp-t+1),V_f_a(t,Nx_Lp-t+1)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
        
        % ��Nx_Lp-t��1���ռ�ڵ㴦��ز�������
        for x=Nx_Lp-t:-1:1
            P_a(t,x)=P_a(t,x+1)-rho_m(t,x+1)*g*dx(x)*cosd(theta(x)); % ����ѹ������ֵ��Pa
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x+1)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x+1)*gamma_l_a(t,x); % ���ջ�Һ����
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            
            err_NodePressure=1; % ����ѹ��������
            COUNT_NodePressure=0; % ����ѹ������������ֵ
            while abs(err_NodePressure)>epsilon_e && COUNT_NodePressure<epsilon_t
                COUNT_NodePressure=COUNT_NodePressure+1; % ����ѹ����������
                
                % ���߷��������������
                alpha_s_ass1=alpha_s(t,x+1)+0.001; % ���������������ֵ1
                alpha_s_ass2=alpha_s(t,x+1)+10000; % ���������������ֵ2
                err_NodeEg=abs(alpha_s_ass1-alpha_s_ass2); % ������������������
                COUNT_NodeEg=0; % ��������
                while abs(err_NodeEg)>epsilon_e && COUNT_NodeEg<epsilon_t
                    COUNT_NodeEg=COUNT_NodeEg+1;
                    
                    % �����������Ϊalpha_s_ass1ʱ������ֵ��
                    int1=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass1)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x)));
                    V_s_ass1=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int1)/(rho_s*alpha_s_ass1); % �����ٶȣ�m/s
                    Va_s_ass1=V_s_ass1*alpha_s_ass1; % ���������٣�m/s
                    alpha_f_ass1=1-alpha_s_ass1; % ��ĭ�������
                    rho_f_ass1=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int1=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass1*alpha_f_ass1)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass1=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int1)/(rho_f_ass1*alpha_f_ass1); % ��ĭ�ٶȣ�m/s
                    Va_f_ass1=V_f_ass1*alpha_f_ass1; % ��ĭ������٣�m/s
                    V_m_ass1=Va_s_ass1+Va_f_ass1; % ���ջ�����ٶȣ�m/s
                    Vsr_ass1=12*(mu_f_a(t,x)/(rho_f_ass1*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass1)/rho_f_ass1)*((rho_f_ass1*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass1_new=Va_s_ass1/(C0*V_m_ass1-Vsr_ass1); % ���������������ֵ
                    
                    Y1=alpha_s_ass1_new-alpha_s_ass1; % ����ĺ��������Ľ������ʵ�����������
                    
                    % �����������Ϊalpha_s_ass2ʱ������ֵ��
                    int2=-dx(x)/(2*dt(t-1))*((rho_s*alpha_s(t,x+1))+(rho_s*alpha_s_ass2)-(rho_s*alpha_s(t-1,x+1))-(rho_s*alpha_s(t-1,x))); % ���������غ㷽����ɢ��ʽ�м�ֵ����
                    V_s_ass2=(rho_s*alpha_s(t,x+1)*V_s(t,x+1)+int2)/(rho_s*alpha_s_ass2); % �����ٶȣ�m/s
                    Va_s_ass2=V_s_ass2*alpha_s_ass2; % ���������٣�m/s
                    alpha_f_ass2=1-alpha_s_ass2; % ��ĭ�������
                    rho_f_ass2=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ��ĭ�ܶȣ�kg/m^3
                    int2=-dx(x)/(2*dt(t-1))*((rho_f_a(t,x+1)*alpha_f_a(t,x+1))+(rho_f_ass2*alpha_f_ass2)-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)));
                    V_f_ass2=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int2)/(rho_f_ass2*alpha_f_ass2); % ��ĭ�ٶȣ�m/s
                    Va_f_ass2=V_f_ass2*alpha_f_ass2; % ��ĭ������٣�m/s
                    V_m_ass2=Va_s_ass2+Va_f_ass2; % ���ջ�����ٶȣ�m/s
                    Vsr_ass2=12*(mu_f_a(t,x)/(rho_f_ass2*D_s))*((1+0.0073*g*D_s*((rho_s-rho_f_ass2)/rho_f_ass2)*((rho_f_ass2*D_s/mu_f_a(t,x))^2))^0.5-1); % ��м����ĩ�٣�m/s��Chien��
                    alpha_s_ass2_new=Va_s_ass2/(C0*V_m_ass2-Vsr_ass2); % ���������������ֵ
                    Y2=alpha_s_ass2_new-alpha_s_ass2; % ����ĺ��������Ľ������ʵ�����������
                    
                    % ���߷���������������
                    alpha_s_ass3=alpha_s_ass2-Y2*(alpha_s_ass2-alpha_s_ass1)/(Y2-Y1); % �µĹ��������������ֵalpha_s_ass3
                    err_NodeEg=abs(alpha_s_ass3-alpha_s_ass2); % ������������������
                    alpha_s_ass1=alpha_s_ass2; % �µĹ��������������ֵ1
                    alpha_s_ass2=alpha_s_ass3; % �µĹ��������������ֵ2
                end
                
                alpha_s(t,x)=alpha_s_ass1; % ���������õ�����ʵ�����������ֵ����alpha_s(t,x)
                
                if alpha_s(t,x)<1e-4 % �������������С��һ��ֵʱ����Ϊ���������Ϊ0�����ڷ�ֹ���ֺ����ļ������
                    alpha_s(t,x)=0; % �����������
                    V_s(t,x)=0; % �����ٶȣ�m/s
                    Va_s(t,x)=0; % �������ٶȣ�m/s
                    mu_s(t,x)=0; % ����ճ�ȣ�Pa*s
                    Vsr(t,x)=0; % ��м����ĩ�٣�m/s
                    
                    alpha_f_a(t,x)=1; % ��ĭ�������                   
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
                    mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
                    int=-dx(x)/(2*dt(t-1))*(rho_f_a(t,x)*alpha_f_a(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)-rho_f_a(t-1,x)*alpha_f_a(t-1,x)-rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1));
                    V_f_a(t,x)=(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_f_a(t,x+1)+int)/(rho_f_a(t,x)*alpha_f_a(t,x)); % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                else
                    V_s(t,x)=V_s_ass2; % �����ٶȣ�m/s
                    Va_s(t,x)=V_s(t,x)*alpha_s(t,x); % ���������٣�m/s
                    Vsr(t,x)=Vsr_ass2; % ��м����ĩ�٣�m/s
                    
                    rho_f_a(t,x)=rho_f_ass2; % ��ĭ�ܶȣ�kg/m^3
                    alpha_f_a(t,x)=1-alpha_s(t,x); % ��ĭ�������
                    alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
                    alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
                    V_f_a(t,x)=V_f_ass2; % ��ĭ�ٶȣ�m/s
                    V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
                    V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
                    Va_f_a(t,x)=V_f_a(t,x)*alpha_f_a(t,x); % ��ĭ������٣�m/s
                    Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
                    Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
                end
                
                mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
                V_m(t,x)=alpha_s(t,x)*V_s(t,x)+alpha_f_a(t,x)*V_f_a(t,x); % ���ջ�����ٶȣ�m/s
                rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
                mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa`s
                [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
                
                M1=-(((rho_f_a(t,x)*alpha_f_a(t,x)*V_f_a(t,x)+rho_s*alpha_s(t,x)*V_s(t,x)+rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)+rho_s*alpha_s(t,x+1)*V_s(t,x+1))-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)+rho_s*alpha_s(t-1,x)*V_s(t-1,x)+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=-((rho_f_a(t,x)*alpha_f_a(t,x)*V_g_a(t,x)^2+rho_s*alpha_s(t,x)*V_s(t,x)^2+rho_f_a(t-1,x)*alpha_f_a(t-1,x)*V_g_a(t-1,x)^2+rho_s*alpha_s(t-1,x)*V_s(t-1,x)^2)-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)*V_g_a(t,x+1)^2+rho_s*alpha_s(t,x+1)*V_s(t,x+1)^2+rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)*V_g_a(t-1,x+1)^2+rho_s*alpha_s(t-1,x+1)*V_s(t-1,x+1)^2))/2;
                M3=dx(x)*((-(rho_f_a(t,x)*alpha_f_a(t,x)+rho_s*alpha_s(t,x))*g*cosd(theta(x))-Ff_a(t,x))+(-(rho_f_a(t-1,x)*alpha_f_a(t-1,x)+rho_s*alpha_s(t-1,x))*g*cosd(theta(x))-Ff_a(t-1,x))+(-(rho_f_a(t,x+1)*alpha_f_a(t,x+1)+rho_s*alpha_s(t,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1))+(-(rho_f_a(t-1,x+1)*alpha_f_a(t-1,x+1)+rho_s*alpha_s(t-1,x+1))*g*cosd(theta(x+1))-Ff_a(t,x+1)))/4;
                P_new=P_a(t,x+1)+M1+M2+M3; % ����ѹ������ֵ��Pa

                err_NodePressure=abs(P_new-P_a(t,x))/P_a(t,x); % ����ѹ��������
                P_a(t,x)=P_new; % �µĻ���ѹ������ֵ��Pa
            end
        end
        
        err_OutPressure=abs(P_a(t,1)-OutPressure)/OutPressure; % ����ѹ��������
        if (P_a(t,1)-OutPressure)>0 % ���ݳ���ѹ�������������Ի��չܵ�ѹ������ֵ���е���
            P_a(t,Nx_Lp-t+1)=P_a(t,Nx_Lp-t+1)-(P_a(t,1)-OutPressure)/2; % �µĻ��չܵ�ѹ������ֵ��Pa
        else
            P_a(t,Nx_Lp-t+1)=P_a(t,Nx_Lp-t+1)-(P_a(t,1)-OutPressure)/2*0.3; % �µĻ��չܵ�ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Lp-t+2��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Lp-t+2:1:Nx
        P_a_ass(t,x)=P_a(t,x-1)+rho_f_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
        
        err_AnnPressure=1; % ����ѹ��������
        COUNT_AnnPressure=0; % ����ѹ������������ֵ
        while abs(err_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
            COUNT_AnnPressure=COUNT_AnnPressure+1; % ����ѹ����������
            
            rho_g_a(t,x)=DensityG(T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_a(t,x)=DensityL(rho_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һ�ܶȣ�kg/m^3
            mu_g_a(t,x)=RheologyG(T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_a(t,x)=RheologyL(mu_l_0,T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���ջ�Һճ�ȣ�Pa*s
            alpha_f_a(t,x)=1; % ������ĭ����
            gamma_g_a(t,x)=(Qm_g_0/rho_g_a(t,x))/(Qm_g_0/rho_g_a(t,x)+Qm_l_0/rho_l_a(t,x)); % ��ĭ����
            gamma_l_a(t,x)=1-gamma_g_a(t,x); % Һ��������
            alpha_g_a(t,x)=alpha_f_a(t,x)*gamma_g_a(t,x); % �������庬��
            alpha_l_a(t,x)=alpha_f_a(t,x)*gamma_l_a(t,x); % ���ջ�Һ����
            alpha_s(t,x)=0; % �����������
            rho_f_a(t,x)=rho_g_a(t,x)*gamma_g_a(t,x)+rho_l_a(t,x)*gamma_l_a(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_a(t,x)=mu_g_a(t,x)*gamma_g_a(t,x)+mu_l_a(t,x)*gamma_l_a(t,x); % ������ĭճ�ȣ�Pa*s
            V_s(t,x)=0; % �����ٶȣ�m/s
            Va_s(t,x)=0; % ���������٣�m/s
            mu_s(t,x)=mu_f_a(t,x); % ����ճ�ȣ�Pa*s
            Vsr(t,x)=0; % ɰ������ĩ�٣�m/s
            Va_f_a(t,x)=0;%Qm_f_0/(rho_f_a(t,x)*A_a(t,x)); % ������ĭ������٣�m/s
            Va_g_a(t,x)=Va_f_a(t,x); % �������������٣�m/s
            Va_l_a(t,x)=Va_f_a(t,x); % ���ջ�Һ������٣�m/s
            V_f_a(t,x)=Va_f_a(t,x)/alpha_f_a(t,x); % ������ĭ���٣�m/s
            V_g_a(t,x)=V_f_a(t,x); % �����������٣�m/s
            V_l_a(t,x)=V_f_a(t,x); % ���ջ�Һ���٣�m/s
            
            V_m(t,x)=Va_s(t,x)+Va_f_a(t,x); % ���ջ�����ٶȣ�m/s
            rho_m(t,x)=alpha_s(t,x)*rho_s+alpha_f_a(t,x)*rho_f_a(t,x); % ���ջ�����ܶȣ�kg/m^3
            mu_m(t,x)=alpha_s(t,x)*mu_s(t,x)+alpha_f_a(t,x)*mu_f_a(t,x); % ���ջ����ճ�ȣ�Pa*s
            [Ff_a(t,x),f_a(t,x),Re_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_m(t,x),V_m(t,x),mu_m(t,x),D_h(t,x),h_a,rho_f_a(t,x),V_f_a(t,x)); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_a(t,x)=-rho_m(t,x)*V_m(t,x)^2+P_a(t,x-1)+rho_m(t,x-1)*V_m(t,x-1)^2+((rho_m(t,x)*g*cosd(theta(x))+Ff_a(t,x)+rho_m(t,x-1)*g*cosd(theta(x-1))+Ff_a(t,x-1))*dx(x-1))/2; % ����ѹ����Pa
            
            err_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ���㻷��ѹ������ֵ�����ֵ֮���������
            P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
        end
    end
end

%% ��1��Ntʱ��ڵ���ͷѹ������
for t=1:1:Nt
    Qv_nozzle(t)=Qm_f_0/rho_f_a(t,Nx_Lp-t+1); % ������ĭ���������m^3/s
    [delta_P_SWT(t),V_nozzle(t)]=PressureDrop_SandWashingTool(C,D_nozzle,N_nozzle,Qv_nozzle(t),rho_f_a(t,Nx_Lp-t+1)); % �����ɰ����ѹ����Pa����������������V_nozzle��m/s��
end

%% ��1��Ntʱ��ڵ���ز������㣨�������ڣ�
for t=1:1:Nt
    % �ܵף���Nx_Lp-t+1���ռ�ڵ㣩����ز�������
    P_ct(t,Nx_Lp-t+1)=P_a(t,Nx_Lp-t+1)+delta_P_SWT(t); % ����ѹ����Pa
    rho_g_ct(t,Nx_Lp-t+1)=DensityG(T_ct(t,Nx_Lp-t+1),P_ct(t,Nx_Lp-t+1)); % ���������ܶȣ�kg/m^3
    rho_l_ct(t,Nx_Lp-t+1)=DensityL(rho_l_0,T_0,P_0,T_ct(t,Nx_Lp-t+1),P_ct(t,Nx_Lp-t+1)); % ���ڻ�Һ�ܶȣ�kg/m^3
    mu_g_ct(t,Nx_Lp-t+1)=RheologyG(T_ct(t,Nx_Lp-t+1),P_ct(t,Nx_Lp-t+1)); % ��������ճ�ȣ�Pa*s
    mu_l_ct(t,Nx_Lp-t+1)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,Nx_Lp-t+1),P_ct(t,Nx_Lp-t+1)); % ���ڻ�Һճ�ȣ�Pa*s
    alpha_f_ct(t,Nx_Lp-t+1)=1; % ������ĭ����
    gamma_g_ct(t,Nx_Lp-t+1)=(Qm_g_0/rho_g_ct(t,Nx_Lp-t+1))/(Qm_g_0/rho_g_ct(t,Nx_Lp-t+1)+Qm_l_0/rho_l_ct(t,Nx_Lp-t+1)); % ������ĭ����
    gamma_l_ct(t,Nx_Lp-t+1)=1-gamma_g_ct(t,Nx_Lp-t+1); % ����Һ��������
    alpha_g_ct(t,Nx_Lp-t+1)=alpha_f_ct(t,Nx_Lp-t+1)*gamma_g_ct(t,Nx_Lp-t+1); % �������庬��
    alpha_l_ct(t,Nx_Lp-t+1)=alpha_f_ct(t,Nx_Lp-t+1)*gamma_l_ct(t,Nx_Lp-t+1); % ���ڻ�Һ����
    rho_f_ct(t,Nx_Lp-t+1)=rho_g_ct(t,Nx_Lp-t+1)*gamma_g_ct(t,Nx_Lp-t+1)+rho_l_ct(t,Nx_Lp-t+1)*gamma_l_ct(t,Nx_Lp-t+1); % ������ĭ�ܶȣ�kg/m^3
    mu_f_ct(t,Nx_Lp-t+1)=mu_g_ct(t,Nx_Lp-t+1)*gamma_g_ct(t,Nx_Lp-t+1)+mu_l_ct(t,Nx_Lp-t+1)*gamma_l_ct(t,Nx_Lp-t+1); % ������ĭճ�ȣ�Pa*s
    V_f_ct(t,Nx_Lp-t+1)=Qm_f_0/(rho_f_ct(t,Nx_Lp-t+1)*A_ct(t,Nx_Lp-t+1)); % �����������٣�m/s
    [Ff_ct(t,Nx_Lp-t+1),f_ct(t,Nx_Lp-t+1),Re_ct(t,Nx_Lp-t+1),flow_pattern_ct(t,Nx_Lp-t+1)]=Friction_coiledtubing(rho_f_ct(t,Nx_Lp-t+1),V_f_ct(t,Nx_Lp-t+1),mu_f_ct(t,Nx_Lp-t+1),D_ct_i(t,Nx_Lp-t+1),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    % ��Nx_Lp-t��1���ռ�ڵ㴦��ز�������
    for x=Nx_Lp-t:-1:1
        P_ct_ass(t,x)=P_ct(t,x+1)-rho_f_ct(t,x+1)*g*cosd(theta(x+1))*dx(x); % ����ѹ������ֵ��Pa
        
        err_DriPipePressure=1; % ����ѹ��������
        COUNT_DriPipePressure=0; % ����ѹ������������ֵ
        while abs(err_DriPipePressure)>epsilon_e && COUNT_DriPipePressure<epsilon_t
            COUNT_DriPipePressure=COUNT_DriPipePressure+1; % ����ѹ����������
            
            rho_g_ct(t,x)=DensityG(T_ct(t,x),P_ct_ass(t,x)); % ���������ܶȣ�kg/m^3
            rho_l_ct(t,x)=DensityL(rho_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һ�ܶȣ�kg/m^3
            mu_g_ct(t,x)=RheologyG(T_ct(t,x),P_ct_ass(t,x)); % ��������ճ�ȣ�Pa*s
            mu_l_ct(t,x)=RheologyL(mu_l_0,T_0,P_0,T_ct(t,x),P_ct_ass(t,x)); % ���ڻ�Һճ�ȣ�Pa*s
            alpha_f_ct(t,x)=1; % ������ĭ����
            gamma_g_ct(t,x)=(Qm_g_0/rho_g_ct(t,x))/(Qm_g_0/rho_g_ct(t,x)+Qm_l_0/rho_l_ct(t,x)); % ������ĭ����
            gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
            alpha_g_ct(t,x)=alpha_f_ct(t,x)*gamma_g_ct(t,x); % �������庬��
            alpha_l_ct(t,x)=alpha_f_ct(t,x)*gamma_l_ct(t,x); % ���ڻ�Һ����
            rho_f_ct(t,x)=rho_g_ct(t,x)*gamma_g_ct(t,x)+rho_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭ�ܶȣ�kg/m^3
            mu_f_ct(t,x)=mu_g_ct(t,x)*gamma_g_ct(t,x)+mu_l_ct(t,x)*gamma_l_ct(t,x); % ������ĭճ�ȣ�Pa*s
            V_f_ct(t,x)=Qm_f_0/(rho_f_ct(t,x)*A_ct(t,x)); % ������ĭ���٣�m/s
            [Ff_ct(t,x),f_ct(t,x),Re_ct(t,x),flow_pattern_ct(t,x)]=Friction_coiledtubing(rho_f_ct(t,x),V_f_ct(t,x),mu_f_ct(t,x),D_ct_i(t,x),h_ct); % �������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
            P_ct(t,x)=-rho_f_ct(t,x)*V_f_ct(t,x)^2+P_ct(t,x+1)+rho_f_ct(t,x+1)*V_f_ct(t,x+1)^2-((rho_f_ct(t,x)*g*cosd(theta(x))-Ff_ct(t,x)+rho_f_ct(t,x+1)*g*cosd(theta(x+1))-Ff_ct(t,x+1))*dx(x))/2; % ����ѹ����Pa
            
            err_DriPipePressure=abs(P_ct(t,x)-P_ct_ass(t,x))/P_ct_ass(t,x); % �������ѹ������ֵ�����ֵ֮���������
            P_ct_ass(t,x)=P_ct(t,x); % �µĹ���ѹ������ֵ��Pa
        end
    end
    
    % ��Nx_Lp-t+2��Nx���ռ�ڵ㴦��ز�������
    for x=Nx_Lp-t+2:1:Nx
        alpha_g_ct(t,x)=alpha_g_a(t,x); % �������庬��
        alpha_l_ct(t,x)=alpha_l_a(t,x); % ���ڻ�Һ����
        gamma_g_ct(t,x)=gamma_g_a(t,x); % ������ĭ����
        gamma_l_ct(t,x)=1-gamma_g_ct(t,x); % ����Һ��������
        rho_f_ct(t,x)=rho_m(t,x); % ���������ܶȣ�kg/m^3
        mu_f_ct(t,t)=mu_m(t,x); % ��������ճ�ȣ�Pa*s
        V_f_ct(t,t)=V_m(t,x); % �����������٣�sm/s
        Ff_ct(t,x)=Ff_a(t,x); % �������嵥λ����Ħ��ѹ����Pa/m��
        f_ct(t,x)=f_a(t,x); % �������巶��Ħ������
        Re_ct(t,x)=Re_a(t,x); % ����������ŵ����������̬
        flow_pattern_ct(t,x)=flow_pattern_a(t,x); % ����������̬
        P_ct(t,x)=P_a(t,x); % ����ѹ����Pa
    end
end

%% ��1��Ntʱ��ڵ��̹ܶγ���ѹ������ѹ���㣨�������ڣ�
for t=1:1:Nt
    V_f_0(t)=Qm_f_0/(A_ct(1,1)*rho_f_ct(t,1)); % ��������������٣�m/s
    
    [Ff_ct_0(t),f_ct_0(t),Re_ct_0(t),flow_pattern_ct_0(t)]=Friction_coiledtubing(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct); % ����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_coil(t),f_ct_coil(t),Re_ct_coil(t),flow_pattern_ct_coil(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_r_i); % �̹ܶ����嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    [Ff_ct_goose(t),f_ct_goose(t),Re_ct_goose(t),flow_pattern_ct_goose(t)]=Friction_bentpipe(rho_f_ct(t,1),V_f_0(t),mu_f_ct(t,1),D_ct_i_3,h_ct,D_goose); % �����������嵥λ����Ħ��ѹ����Pa/m��������Ħ�����ӡ���ŵ����������̬
    
    delta_P_wg(t)=Ff_ct_0(t)*L_wg; % ���ڵ�ע��ͷ������Ħ��ѹ����Pa
    delta_P_goose(t)=Ff_ct_goose(t)*L_goose; % ��������Ħ��ѹ����Pa
    delta_P_gr(t)=Ff_ct_0(t)*L_gr; % ����������Ͳ��Ħ��ѹ����Pa
    delta_P_coil(t)=Ff_ct_coil(t)*L_reel(t); % �̹ܶ�Ħ��ѹ����Pa
    
    P_coil(t)=P_ct(t,1)-rho_l_0*g*L_wg+delta_P_wg(t)+rho_l_0*g*L_gr*cosd(theta_gr)+delta_P_goose(t)+delta_P_gr(t); % �̹ܶγ���ѹ����Pa
    P_pump(t)=P_coil(t)+delta_P_coil(t); % ��ѹ��Pa
end

%% ���㾮���ۻ���ɰ��
M_w_tem(1)=0; % ����˲ʱ��ɰ����ֵ��kg
M_w_tot(1)=0; % �����ۻ���ɰ����ֵ��kg
for t=2:1:Nt
    M_w_tem(t)=alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮��˲ʱ��ɰ����kg
    M_w_tot(t)=M_w_tot(t-1)+alpha_s(t,1)*0.25*pi*(D_t_i(t,1)^2-D_ct_o(t,1)^2)*rho_s*V_s(t,1)*dt(t-1); % ��t��ʱ��ڵ㾮���ۻ���ɰ����kg
end

%% ����ĩ̬��Ͳ���ɰŨ��
alpha_s_max=0;
for x=1:1:Nx_Dsp_pen(2)
    if alpha_s(Nt,x) >= alpha_s_max
        alpha_s_max=alpha_s(Nt,x);
    end
end

%% ���㻷����ĭ����ƽ��ֵV_l_mre
V_f_mre=0; % ���շ���ƽ��ֵ��m/s
for t=1:1:Nt
    for x=1:1:Nx
        V_f_mre=V_f_mre+V_f_a(t,x)/(Nt*Nx); % ���շ���ƽ��ֵ��m/s
    end
end

%% ����������ĩ��ƽ��ֵVsr_mre
num=1;
for t=1:1:Nt
    for x=1:1:Nx
        if Vsr(t,x)>0
            VSR(num)=Vsr(t,x); % ��Vsr��������ȡֵ����ĳ���ĩ�٣�m/s
            num=num+1;
        else
        end
    end
end

Vsr_mre=0; % �������ĩ��ƽ��ֵ��m/s
for x=1:1:num-1
    Vsr_mre=Vsr_mre+VSR(x)/(num-1); % �������ĩ��ƽ��ֵ��m/s
end

%% �ж��Ƿ���Ч��ɰ�����շ���ƽ��ֵ����2��ɰ������ĩ��ƽ��ֵ�����Ƿ���ɳ�ɰ
fprintf("������������̣�\n");
if V_f_mre > 2*Vsr_mre
    fprintf("Valid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
    
    if alpha_s_max == 0
        fprintf("Sand Cleanout Finished!\n"); % �����д���˵����ǰ�����¡���ɳ�ɰ��
    else
        fprintf("Sand Cleanout UnFinished!\n"); % �����д���˵����ǰ�����¡�δ��ɳ�ɰ��
    end
    
    else
    fprintf("InValid Sand Cleanout!\n"); % �����д���˵����ǰ������Ϊ����Ч��ɰ��
end

%% ��ͲECD����
for t=1:1:Nt
    for x=2:1:Nx
        ECD_a(t,x)=P_a(x)/(g*Depth(x)); % ����ECD��kg/m^3
        ECD_ct(t,x)=P_ct(x)/(g*Depth(x)); % ����ECD��kg/m^3
    end
    ECD_a(t,1)=ECD_a(t,2)-((ECD_a(t,3)-ECD_a(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
    ECD_ct(t,1)=ECD_ct(t,2)-((ECD_ct(t,3)-ECD_ct(t,2))/dx(2))*dx(1); % ����ECD��kg/m^3
end

%% ���ݴ洢�������������
ANS_Nt_8=Nt; % ʱ��ڵ���
ANS_Time_8=Time(1:Nt); % ���������ʱ�䣬s
ANS_alpha_g_a_8=alpha_g_a(1:Nt,:); % ���������������
ANS_alpha_l_a_8=alpha_l_a(1:Nt,:); % ����Һ���������
ANS_alpha_s_8=alpha_s(1:Nt,:); % ���ຬ��
ANS_P_ct_8=P_ct(1:Nt,:); % ����ѹ����Pa��
ANS_delta_P_SWT_8=delta_P_SWT(1:Nt); % ��ɰ����ѹ����Pa��
ANS_P_a_8=P_a(1:Nt,:); % ����ѹ����Pa��
ANS_T_a_8=T_a(1:Nt,:); % �����¶ȣ��棩
ANS_T_ct_8=T_ct(1:Nt,:); % �����¶ȣ��棩
ANS_P_coil_8=P_coil(1:Nt); % �̹ܶγ���ѹ����Pa��
ANS_P_pump_8=P_pump(1:Nt); % ��ѹ��Pa��
ANS_M_w_tem_8=M_w_tem(1:Nt); % ����˲ʱ��ɰ����kg��
ANS_M_w_tot_8=M_w_tot(1:Nt); % �����ۻ���ɰ����kg��
ANS_Va_s_8=Va_s(1:Nt,:); % ��м����ٶȣ�m/s��
ANS_Va_f_a_8=Va_f_a(1:Nt,:); % ������ĭ������٣�m/s��
ANS_V_s_8=V_s(1:Nt,:); % ��м�����ٶȣ�m/s��
ANS_V_f_a_8=V_f_a(1:Nt,:); % ������ĭ���٣�m/s��
ANS_alpha_g_ct_8=alpha_g_ct(1:Nt,:); % ���������������
ANS_alpha_l_ct_8=alpha_l_ct(1:Nt,:); % ����Һ���������
ANS_gamma_g_ct_8=gamma_g_ct(1:Nt,:); % ������ĭ����
ANS_gamma_g_a_8=gamma_g_a(1:Nt,:); % ������ĭ����
ANS_ECD_a_8=ECD_a(1:Nt,:); % ����ECD��kg/m^3��
ANS_L_coil_8=L_coil(1:Nt); % ���������m��



%% ����
Nt=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7+ANS_Nt_8; % ʱ��ڵ���
for t=1:1:ANS_Nt_1
    Time(t)=ANS_Time_1(t); % ��ʱ�䣬s
    M_w_tot(t)=ANS_M_w_tot_1(t); % �����ۻ���ɰ����kg��
end
for t=ANS_Nt_1+1:1:ANS_Nt_1+ANS_Nt_2
    Time(t)=Time(ANS_Nt_1)+ANS_Time_2(t-ANS_Nt_1); % ��ʱ�䣬s
    M_w_tot(t)=M_w_tot(ANS_Nt_1)+ANS_M_w_tot_2(t-ANS_Nt_1); % �����ۻ���ɰ����kg��
end
for t=ANS_Nt_1+ANS_Nt_2+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2)+ANS_Time_3(t-ANS_Nt_1-ANS_Nt_2); % ��ʱ�䣬s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2)+ANS_M_w_tot_3(t-ANS_Nt_1-ANS_Nt_2); % �����ۻ���ɰ����kg��
end

for t=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3)+ANS_Time_4(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3)); % ��ʱ�䣬s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3)+ANS_M_w_tot_4(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3)); % �����ۻ���ɰ����kg��
end
for t=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4)+ANS_Time_5(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4)); % ��ʱ�䣬s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4)+ANS_M_w_tot_5(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4)); % �����ۻ���ɰ����kg��
end
for t=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5)+ANS_Time_6(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5)); % ��ʱ�䣬s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5)+ANS_M_w_tot_6(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5)); % �����ۻ���ɰ����kg��
end
for t=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6)+ANS_Time_7(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6)); % ��ʱ�䣬s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6)+ANS_M_w_tot_7(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6)); % �����ۻ���ɰ����kg��
end
for t=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7+1:1:ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7+ANS_Nt_8
    Time(t)=Time(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7)+ANS_Time_8(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7)); % ��ʱ�䣬s
    M_w_tot(t)=M_w_tot(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7)+ANS_M_w_tot_8(t-(ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7)); % �����ۻ���ɰ����kg��
end


alpha_g_a=[ANS_alpha_g_a_1;ANS_alpha_g_a_2;ANS_alpha_g_a_3;ANS_alpha_g_a_4;ANS_alpha_g_a_5;ANS_alpha_g_a_6;ANS_alpha_g_a_7;ANS_alpha_g_a_8]; % ���������������
alpha_l_a=[ANS_alpha_l_a_1;ANS_alpha_l_a_2;ANS_alpha_l_a_3;ANS_alpha_l_a_4;ANS_alpha_l_a_5;ANS_alpha_l_a_6;ANS_alpha_l_a_7;ANS_alpha_l_a_8]; % ����Һ���������
alpha_s=[ANS_alpha_s_1;ANS_alpha_s_2;ANS_alpha_s_3;ANS_alpha_s_4;ANS_alpha_s_5;ANS_alpha_s_6;ANS_alpha_s_7;ANS_alpha_s_8]; % ���ຬ��
P_ct=[ANS_P_ct_1;ANS_P_ct_2;ANS_P_ct_3;ANS_P_ct_4;ANS_P_ct_5;ANS_P_ct_6;ANS_P_ct_7;ANS_P_ct_8]; % ����ѹ����Pa��
delta_P_SWT=[ANS_delta_P_SWT_1,ANS_delta_P_SWT_2,ANS_delta_P_SWT_3,ANS_delta_P_SWT_4,ANS_delta_P_SWT_5,ANS_delta_P_SWT_6,ANS_delta_P_SWT_7,ANS_delta_P_SWT_8]; % ��ɰ����ѹ����Pa��
P_a=[ANS_P_a_1;ANS_P_a_2;ANS_P_a_3;ANS_P_a_4;ANS_P_a_5;ANS_P_a_6;ANS_P_a_7;ANS_P_a_8]; % ����ѹ����Pa��
T_a=[ANS_T_a_1;ANS_T_a_2;ANS_T_a_3;ANS_T_a_4;ANS_T_a_5;ANS_T_a_6;ANS_T_a_7;ANS_T_a_8]; % �����¶ȣ��棩
T_ct=[ANS_T_ct_1;ANS_T_ct_2;ANS_T_ct_3;ANS_T_ct_4;ANS_T_ct_5;ANS_T_ct_6;ANS_T_ct_7;ANS_T_ct_8]; % �����¶ȣ��棩
P_coil=[ANS_P_coil_1,ANS_P_coil_2,ANS_P_coil_3,ANS_P_coil_4,ANS_P_coil_5,ANS_P_coil_6,ANS_P_coil_7,ANS_P_coil_8]; % �̹ܶγ���ѹ����Pa��
P_pump=[ANS_P_pump_1,ANS_P_pump_2,ANS_P_pump_3,ANS_P_pump_4,ANS_P_pump_5,ANS_P_pump_6,ANS_P_pump_7,ANS_P_pump_8]; % ��ѹ��Pa��
M_w_tem=[ANS_M_w_tem_1,ANS_M_w_tem_2,ANS_M_w_tem_3,ANS_M_w_tem_4,ANS_M_w_tem_5,ANS_M_w_tem_6,ANS_M_w_tem_7,ANS_M_w_tem_8]; % ����˲ʱ��ɰ����kg��
Va_s=[ANS_Va_s_1;ANS_Va_s_2;ANS_Va_s_3;ANS_Va_s_4;ANS_Va_s_5;ANS_Va_s_6;ANS_Va_s_7;ANS_Va_s_8]; % ��м����ٶȣ�m/s��
Va_f_a=[ANS_Va_f_a_1;ANS_Va_f_a_2;ANS_Va_f_a_3;ANS_Va_f_a_4;ANS_Va_f_a_5;ANS_Va_f_a_6;ANS_Va_f_a_7;ANS_Va_f_a_8]; % ������ĭ������٣�m/s��
V_s=[ANS_V_s_1;ANS_V_s_2;ANS_V_s_3;ANS_V_s_4;ANS_V_s_5;ANS_V_s_6;ANS_V_s_7;ANS_V_s_8]; % ��м�����ٶȣ�m/s��
V_f_a=[ANS_V_f_a_1;ANS_V_f_a_2;ANS_V_f_a_3;ANS_V_f_a_4;ANS_V_f_a_5;ANS_V_f_a_6;ANS_V_f_a_7;ANS_V_f_a_8]; % ������ĭ���٣�m/s��
alpha_g_ct=[ANS_alpha_g_ct_1;ANS_alpha_g_ct_2;ANS_alpha_g_ct_3;ANS_alpha_g_ct_4;ANS_alpha_g_ct_5;ANS_alpha_g_ct_6;ANS_alpha_g_ct_7;ANS_alpha_g_ct_8]; % ���������������
alpha_l_ct=[ANS_alpha_l_ct_1;ANS_alpha_l_ct_2;ANS_alpha_l_ct_3;ANS_alpha_l_ct_4;ANS_alpha_l_ct_5;ANS_alpha_l_ct_6;ANS_alpha_l_ct_7;ANS_alpha_l_ct_8]; % ����Һ���������
gamma_g_ct=[ANS_gamma_g_ct_1;ANS_gamma_g_ct_2;ANS_gamma_g_ct_3;ANS_gamma_g_ct_4;ANS_gamma_g_ct_5;ANS_gamma_g_ct_6;ANS_gamma_g_ct_7;ANS_gamma_g_ct_8]; % ������ĭ����
gamma_g_a=[ANS_gamma_g_a_1;ANS_gamma_g_a_2;ANS_gamma_g_a_3;ANS_gamma_g_a_4;ANS_gamma_g_a_5;ANS_gamma_g_a_6;ANS_gamma_g_a_7;ANS_gamma_g_a_8]; % ������ĭ����
ECD_a=[ANS_ECD_a_1;ANS_ECD_a_2;ANS_ECD_a_3;ANS_ECD_a_4;ANS_ECD_a_5;ANS_ECD_a_6;ANS_ECD_a_7;ANS_ECD_a_8]; % ����ECD��kg/m^3��
L_coil=[ANS_L_coil_1,ANS_L_coil_2,ANS_L_coil_3,ANS_L_coil_4,ANS_L_coil_5,ANS_L_coil_6,ANS_L_coil_7,ANS_L_coil_8]; % ���������m��

%% ��ͼ
%%
fig1 = figure; % ��Ͳѹ��VS����VSʱ��
for t = 1:1:Nt
    plot(P_a(t,:)/10^6,Depth,P_ct(t,:)/10^6,Depth,[0,60],[L_coil(t),L_coil(t)],'--','LineWidth',2);
    legend('����ѹ��','����ѹ��','����������','FontName','����','Location','Best');
    xlabel('ѹ����MPa��','FontName','����','FontSize',12);
    ylabel('���m��','FontName','����','FontSize',12);
    set(gca,'fontsize',12);
    set(gca,'YDir','reverse');
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    if t==1
        title({['��ҵ״̬������������ɰ����'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['��ҵ״̬����ϴ���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['��ҵ״̬�����϶���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['��ҵ״̬������ѭ��1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['��ҵ״̬������������1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['��ҵ״̬����ϴ���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['��ҵ״̬�����϶���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['��ҵ״̬������ѭ��2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    else
        title({['��ҵ״̬�����������'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    end
    
    frame1 = getframe(fig1);
    im{t} = frame2im(frame1);
end

filename1 = '��Ͳѹ��VS����VSʱ��.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename1,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename1,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig2 = figure; % ��ĭɰ���ٶ�VS����VSʱ��
for t = 1:1:Nt
    subplot(1,2,1);
    plot(V_s(t,:),Depth); % ��м�����ٶȣ�m/s��
    xlabel('��м�����ٶȣ�m/s��','FontName','����','FontSize',10);
    ylabel('���m��','FontName','����','FontSize',10);
    set(gca,'fontsize',12);
    set(gca,'YDir','reverse');
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    if t==1
        title({['��ҵ״̬������������ɰ����'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['��ҵ״̬����ϴ���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['��ҵ״̬�����϶���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['��ҵ״̬������ѭ��1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['��ҵ״̬������������1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['��ҵ״̬����ϴ���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['��ҵ״̬�����϶���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['��ҵ״̬������ѭ��2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    else
        title({['��ҵ״̬�����������'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    end
    
    subplot(1,2,2);
    plot(V_f_a(t,:),Depth); % ������ĭ���٣�m/s��
    xlabel('������ĭ���٣�m/s��','FontName','����','FontSize',10);
    ylabel('���m��','FontName','����','FontSize',10);
    set(gca,'fontsize',12);
    set(gca,'YDir','reverse');
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    frame2 = getframe(fig2);
    im{t} = frame2im(frame2);
end

filename2 = '��ĭɰ���ٶ�VS����VSʱ��.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename2,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename2,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig3 = figure; % ɰŨ��VS����VSʱ��
for t = 1:1:Nt
    alpha_s1(t,:)=alpha_s(t,:);
    ind=find(alpha_s(t,:)==0);
    alpha_s1(t,ind)=NaN;
    plot(Depth,alpha_s1(t,:),[L_coil(t),L_coil(t)],[0,0.6],'--','LineWidth',2);
    legend('ɰŨ��','����������','FontName','����','Location','northwest');
    
    ylabel('ɰŨ��','FontName','����','FontSize',12);
    xlabel('���m��','FontName','����','FontSize',12);
    set(gca,'fontsize',12);
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    axis([0,6200,0,0.6]);
    
    if t==1
        title({['��ҵ״̬������������ɰ����'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['��ҵ״̬����ϴ���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['��ҵ״̬�����϶���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['��ҵ״̬������ѭ��1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['��ҵ״̬������������1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['��ҵ״̬����ϴ���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['��ҵ״̬�����϶���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['��ҵ״̬������ѭ��2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    else
        title({['��ҵ״̬�����������'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    end
    
    frame3 = getframe(fig3);
    im{t} = frame2im(frame3);
end

filename3 = 'ɰŨ��VS����VSʱ��.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename3,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename3,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig4 = figure; % ��ѹVS����VSʱ��
for t = 1:1:Nt
    plot(Time(t)/60,P_pump(t)/10^6,'-o');
    hold on;
    xlabel('ʱ�䣨min��','FontName','����','FontSize',12);
    ylabel('��ѹ��MPa��','FontName','����','FontSize',12);
    set(gca,'fontsize',12);
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    if t==1
        title({['��ҵ״̬������������ɰ����'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['��ҵ״̬����ϴ���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['��ҵ״̬�����϶���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['��ҵ״̬������ѭ��1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['��ҵ״̬������������1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['��ҵ״̬����ϴ���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['��ҵ״̬�����϶���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['��ҵ״̬������ѭ��2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    else
        title({['��ҵ״̬�����������'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    end
    
    frame4 = getframe(fig4);
    im{t} = frame2im(frame4);
end

filename4 = '��ѹVS����VSʱ��.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename4,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename4,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig5 = figure; % �̹ܶγ���ѹ��VS����VSʱ��
for t = 1:1:Nt
    plot(Time(t)/60,P_ct(t,1)/10^6,'-o');
    hold on;
    xlabel('ʱ�䣨min��','FontName','����','FontSize',12);
    ylabel('�̹ܶγ���ѹ����MPa��','FontName','����','FontSize',12);
    set(gca,'fontsize',12);
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    if t==1
        title({['��ҵ״̬������������ɰ����'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['��ҵ״̬����ϴ���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['��ҵ״̬�����϶���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['��ҵ״̬������ѭ��1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['��ҵ״̬������������1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['��ҵ״̬����ϴ���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['��ҵ״̬�����϶���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['��ҵ״̬������ѭ��2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    else
        title({['��ҵ״̬�����������'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    end
    
    frame5 = getframe(fig5);
    im{t} = frame2im(frame5);
end

filename5 = '�̹ܶγ���ѹ��VS����VSʱ��.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename5,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename5,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig6 = figure; % ���վ���ѹ��VS����VSʱ��
for t = 1:1:Nt
    plot(Time(t)/60,P_a(t,Nx)/10^6,'-o');
    hold on;
    xlabel('ʱ�䣨min��','FontName','����','FontSize',12);
    ylabel('���վ���ѹ����MPa��','FontName','����','FontSize',12);
    set(gca,'fontsize',12);
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    if t==1
        title({['��ҵ״̬������������ɰ����'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['��ҵ״̬����ϴ���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['��ҵ״̬�����϶���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['��ҵ״̬������ѭ��1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['��ҵ״̬������������1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['��ҵ״̬����ϴ���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['��ҵ״̬�����϶���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['��ҵ״̬������ѭ��2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    else
        title({['��ҵ״̬�����������'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    end
    
    frame6 = getframe(fig6);
    im{t} = frame2im(frame6);
end

filename6 = '���վ���ѹ��VS����VSʱ��.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename6,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename6,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig7 = figure; % ����ECDVS����VSʱ��
for t = 1:1:Nt
    plot(P_a(t,13:Nx)./(g*Depth(13:Nx)),Depth(13:Nx),[max(P_a(t,13:Nx)./(g*Depth(13:Nx)))-30,max(P_a(t,13:Nx)./(g*Depth(13:Nx)))+30],[L_coil(t),L_coil(t)],'--','LineWidth',2);
    legend('����ECD','����������','FontName','����','Location','Best');
    xlabel('����ECD��kg/m3��','FontName','����','FontSize',12);
    ylabel('���m��','FontName','����','FontSize',12);
    set(gca,'fontsize',12);
    set(gca,'YDir','reverse');
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    if t==1
        title({['��ҵ״̬������������ɰ����'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['��ҵ״̬����ϴ���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['��ҵ״̬�����϶���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['��ҵ״̬������ѭ��1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['��ҵ״̬������������1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['��ҵ״̬����ϴ���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['��ҵ״̬�����϶���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['��ҵ״̬������ѭ��2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    else
        title({['��ҵ״̬�����������'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    end
    
    frame7 = getframe(fig7);
    im{t} = frame2im(frame7);
end

filename7 = '����ECDVS����VSʱ��.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 1
        imwrite(A,map,filename7,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename7,'gif','WriteMode','append','DelayTime',0.11);
    end
end

%%
fig8 = figure; % ��Ͳ��ĭ����VS����VSʱ��
for t = 1:1:Nt
    plot(gamma_g_a(t,:),Depth,gamma_g_ct(t,:),Depth,[0,1],[L_coil(t),L_coil(t)],'--','LineWidth',2);
    legend('������ĭ����','������ĭ����','����������','FontName','����','Location','Best');
    xlabel('��ĭ����','FontName','����','FontSize',12);
    ylabel('���m��','FontName','����','FontSize',12);
    set(gca,'fontsize',12);
    set(gca,'YDir','reverse');
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    if t==1
        title({['��ҵ״̬������������ɰ����'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>1 && t<=ANS_Nt_1
        title({['��ҵ״̬����ϴ���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+1 && t<=ANS_Nt_1+ANS_Nt_2
        title({['��ҵ״̬�����϶���1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3
        title({['��ҵ״̬������ѭ��1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4
        title({['��ҵ״̬������������1'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5
        title({['��ҵ״̬����ϴ���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6
        title({['��ҵ״̬�����϶���2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    elseif t>=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+1 && t<=ANS_Nt_1+ANS_Nt_2+ANS_Nt_3+ANS_Nt_4+ANS_Nt_5+ANS_Nt_6+ANS_Nt_7
        title({['��ҵ״̬������ѭ��2'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    else
        title({['��ҵ״̬�����������'],['ʱ�䣺',num2str(round(Time(t)/60)),' min    ','���������',num2str(round(L_coil(t))),' m'],['��ѹ��',num2str(roundn(P_pump(t)/10^6,-2)),' MPa']},'Color','r','FontName','����','FontSize',12);
        pause(0.1);
    end
    
    frame8 = getframe(fig8);
    im{t} = frame2im(frame8);
end

filename8 = '��Ͳ��ĭ����VS����VSʱ��.gif';
for t = 1:1:Nt
    [A,map] = rgb2ind(im{t},256);
    if t == 0
        imwrite(A,map,filename8,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename8,'gif','WriteMode','append','DelayTime',0.11);
    end
end