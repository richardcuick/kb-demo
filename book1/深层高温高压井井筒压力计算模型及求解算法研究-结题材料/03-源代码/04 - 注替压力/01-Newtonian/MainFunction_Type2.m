%% ��������������㣨����������������ѹ����
close all;clear all;clc;
%% �������ݼ�Ԥ�������ʵ�λ�ƣ�
T_0=20; % �������ʣ��ܶȡ���������������¶ȣ���
P_0=0.1*10^5; % �������ʣ��ܶȡ��������������ѹ����Pa
data_rheology_01=csvread('RheologyData_01.csv',1,0); % �����һ������T_0��P_0�������Բ�������
data_rheology_02=csvread('RheologyData_02.csv',1,0); % ����ڶ�������T_0��P_0�������Բ�������
[mu_0_1]=Rheology_Newtonian(data_rheology_01); % �����һ������T_0��P_0��ճ�ȣ�Pa*s��
[mu_0_2]=Rheology_Newtonian(data_rheology_02); % ����ڶ�������T_0��P_0��ճ�ȣ�Pa*s��
rho_0_1=1000; % ��һ������T_0��P_0���ܶȣ�kg/m^3
rho_0_2=1500; % �ڶ�������T_0��P_0���ܶȣ�kg/m^3

g=9.81; % �������ٶȣ�m/s^2
epsilon_e=1*10^(-4); % ���������
epsilon_t=1*10^5; % ��������������

welldepth=4000; % ���m
nx=100; % �ռ�������
Nx=nx+1; % �ռ�ڵ���
for x=1:1:60
    dx(x)=20;%welldepth/nx; % �ռ䲽����m
end
for x=61:1:Nx-1
    dx(x)=(welldepth-20*60)/(Nx-60-1); % �ռ䲽����m
end

Depth(1)=0;
for x=2:1:Nx
    Depth(x)=Depth(x-1)+dx(x-1); % ���m
end

for x=1:1:Nx
    theta(x)=0; % ��б�ǣ���
end

for x=1:1:80
    D_w(x)=140*10^(-3); % ����ֱ����m
end
for x=81:1:Nx
    D_w(x)=140*10^(-3); % ����ֱ����m
end

for x=1:1:80
    D_d_o(x)=114.3*10^(-3); % �����⾶��m
end
for x=81:1:Nx
    D_d_o(x)=114.3*10^(-3); % �����⾶��m
end

for x=1:1:80
    D_d_i(x)=100.53*10^(-3); % �����ھ���m
end
for x=81:1:Nx
    D_d_i(x)=100.53*10^(-3); % �����ھ���m
end

for x=1:1:Nx
    A_d(x)=1/4*pi*D_d_i(x)^2; % ���ڽ������m^2
    A_a(x)=1/4*pi*(D_w(x)^2-D_d_o(x)^2); % ���ս������m^2
end

nt=nx; % ʱ��������
Nt=Nx; % ʱ��ڵ���

Qv_0=1/60; % ��������m^3/s

for x=1:1:Nx
    V_d_0(x)=Qv_0/A_d(x); % ���ڸ��ռ�ڵ��������٣�m/s
    V_a_0(x)=Qv_0/A_a(x); % ���ո��ռ�ڵ��������٣�m/s
end

for x=1:1:Nx-1
    dt_d(x)=dx(x)/V_d_0(x); % ����������������������ʱ�䣬s
    dt_a(x)=dx(x)/V_a_0(x); % ����������������������ʱ�䣬s
end

for t=1:1:(Nt-1)
    dt(t)=dt_d(t); % ʱ�䲽����s
end
for t=Nt:1:(2*Nt-2)
    dt(t)=dt_a(2*Nt-t-1); % ʱ�䲽����s
end

t_inject_d(1)=0;
t_inject_a(1)=0;
for x=2:1:Nx
    t_inject_d(x)=t_inject_d(x-1)+dt_d(x-1); % ���������������ʱ����s
    t_inject_a(x)=t_inject_a(x-1)+dt_a(x-1); % ���滷����������ʱ����s
end
t_total=t_inject_d(Nx)+t_inject_a(Nx); % ��ʱ����s

%%
for t=1:1:(2*Nt-1)
    for x=1:1:Nx
        rho_a_0(t,x)=rho_0_1; % ��ʼʱ�̻��յ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ��ܶȣ�kg/m^3
        rho_d_0(t,x)=rho_0_1; % ��ʼʱ�̹��ڵ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ��ܶȣ�kg/m^3
        mu_a_0(t,x)=mu_0_1; % ��ʼʱ�̻��յ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ�ճ�ȣ�Pa*s
        mu_d_0(t,x)=mu_0_1; % ��ʼʱ�̹��ڵ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ�ճ�ȣ�Pa*s
    end
end

for t=1:1:Nt
    for x=1:1:Nx
        if x<=t
            rho_d_0(t,x)=rho_0_2; % ���ڵ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ��ܶȣ�kg/m^3
            mu_d_0(t,x)=mu_0_2; % ���ڵ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ�ճ�ȣ�Pa*s
        else
            rho_d_0(t,x)=rho_0_1; % ���ڵ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ��ܶȣ�kg/m^3
            mu_d_0(t,x)=mu_0_1; % ���ڵ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ�ճ�ȣ�Pa*s
        end
    end
end

for t=Nt:1:(2*Nt-1)
    for x=Nx:-1:1
        if x>=(2*Nt-t) %% Nx-x<=t-Nt
            rho_a_0(t,x)=rho_0_2; % ���յ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ��ܶȣ�kg/m^3
            mu_a_0(t,x)=mu_0_2; % ���յ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ�ճ�ȣ�Pa*s
        else
            rho_a_0(t,x)=rho_0_1; % ���յ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ��ܶȣ�kg/m^3
            mu_a_0(t,x)=mu_0_1; % ���յ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ�ճ�ȣ�Pa*s
        end
        rho_d_0(t,x)=rho_0_2; % ���ڵ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ��ܶȣ�kg/m^3
        mu_d_0(t,x)=mu_0_2; % ���ڵ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ�ճ�ȣ�Pa*s
    end
end

for t=1:1:(2*Nt-1)
    for x=1:1:Nx
        Qm_d_0(t,x)=Qv_0*rho_d_0(t,x); % ���ڵ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ�����������kg/s
        Qm_a_0(t,x)=Qv_0*rho_a_0(t,x); % ���յ�t��ʱ��ڵ��x���ռ�ڵ㴦����T_0��P_0�µ�����������kg/s
    end
end

C=0.95; % ��������ϵ��
d_nozzle=12.7*10^(-3); % ����ֱ����m
n_nozzle=3; % �������

InPressure=45*10^6; % ���ѹ����Pa

for t=1:1:Nt
    IP(t)=InPressure; % ���ѹ����Pa
end
for t=(Nt+1):1:(2*Nt-1)
    IP(t)=InPressure; % ���ѹ����Pa
end

%% �¶ȼ���
[T_a_T,T_d_T,Depth_T,Time_T]=Temperature(welldepth,(rho_0_1+rho_0_2)/2,(mu_0_1+mu_0_2)/2,Qv_0); % �����¶ȳ������������¶ȣ��棩�������¶ȣ��棩���ռ�ڵ㣨m����ʱ��ڵ㣨s����

%% ��ֵ����ռ�ڵ�Depth��Ӧ�Ļ����¶�T_A�͹����¶�T_D
[r,c]=size(T_a_T);
for t=1:1:c
    for x=1:1:Nx
        for i=1:1:r-1
            if Depth(x)>=Depth_T(i) && Depth(x)<Depth_T(i+1)
                T_A(x,t)=T_a_T(i,t)+(T_a_T(i+1,t)-T_a_T(i,t))*(Depth(x)-Depth_T(i))/(Depth_T(i+1)-Depth_T(i)); % �����꾮Һ�¶ȣ���
                T_D(x,t)=T_d_T(i,t)+(T_d_T(i+1,t)-T_d_T(i,t))*(Depth(x)-Depth_T(i))/(Depth_T(i+1)-Depth_T(i)); % �����꾮Һ�¶ȣ���
            elseif Depth(x)==Depth_T(i+1)
                T_A(x,t)=T_a_T(i+1,t); % �����꾮Һ�¶ȣ���
                T_D(x,t)=T_d_T(i+1,t); % �����꾮Һ�¶ȣ���
            end
        end
    end
end

for t=1:1:(2*Nt-1)
    T_a(t,:)=T_A(:,c); % ���������¶ȣ�����ÿ��ʱ�̶�һ��������
    T_d(t,:)=T_D(:,c); % ���������¶ȣ�����ÿ��ʱ�̶�һ��������
end

%% ��1��ʱ��ڵ���ز������㣨���ڣ�
% ��1���ռ�ڵ㴦��ز�������
P_d(1,1)=InPressure; % �������ѹ����Pa
rho_d(1,1)=Density_TP(rho_d_0(1,1),T_0,P_0,T_d(1,1),P_d(1,1)); % ������������ܶȣ�kg/m^3
Qv_d(1,1)=Qm_d_0(1,1)/rho_d(1,1); % ����������������m^3/s
V_d(1,1)=Qv_d(1,1)/A_d(1); % ��������������٣�m/s
[mu_d(1,1)]=Rheology_TP(mu_d_0(1,1),T_0,P_0,T_d(1,1),P_d(1,1)); % �����������ճ�ȣ�Pa*s
[Ff_d(1,1),flow_pattern_d(1,1)]=Friction_drillpipe(rho_d(1,1),V_d(1,1),mu_d(1,1),D_d_i(1)); % ������ڵ�λ����Ħ��ѹ����Pa/m����������̬������1��������2������3��
Fa_d(1,1)=0; % ������ڵ�λ���ȼ��ٶ�ѹ����Pa/m

% ��2��Nx���ռ�ڵ㴦��ز�������
for x=2:1:Nx
    P_d_ass(1,x)=P_d(1,x-1)+rho_d(1,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
    
    ERROR_PipPressure=1; % ����ѹ��������
    COUNT_PipPressure=0; % ��������
    while abs(ERROR_PipPressure)>epsilon_e && COUNT_PipPressure<epsilon_t
        COUNT_PipPressure=COUNT_PipPressure+1;
        
        rho_d(1,x)=Density_TP(rho_d_0(1,x),T_0,P_0,T_d(1,x),P_d_ass(1,x)); % ���������ܶȣ�kg/m^3
        Qv_d(1,x)=Qm_d_0(1,x)/rho_d(1,x); % �������������m^3/s
        V_d(1,x)=Qv_d(1,x)/A_d(x); % �����������٣�m/s
        [mu_d(1,x)]=Rheology_TP(mu_d_0(1,x),T_0,P_0,T_d(1,x),P_d_ass(1,x)); % ��������ճ�ȣ�Pa*s
        [Ff_d(1,x),flow_pattern_d(1,x)]=Friction_drillpipe(rho_d(1,x),V_d(1,x),mu_d(1,x),D_d_i(x)); % ���ڵ�λ����Ħ��ѹ����Pa/m����������̬������1��������2������3��
        Fa_d(1,x)=abs(0.5*(rho_d(1,x)*V_d(1,x)^2-rho_d(1,x-1)*V_d(1,x-1)^2)); % ���ڵ�λ���ȼ��ٶ�ѹ����Pa/m
        P_d(1,x)=-rho_d(1,x)*V_d(1,x)^2+rho_d(1,x-1)*V_d(1,x-1)^2+P_d(1,x-1)+(((rho_d(1,x)*g*cosd(theta(x))-Ff_d(1,x)-Fa_d(1,x))+(rho_d(1,x-1)*g*cosd(theta(x-1))-Ff_d(1,x-1)-Fa_d(1,x-1)))/2)*dx(x-1); % ����ѹ������ֵ��Pa
        
        ERROR_PipPressure=abs(P_d(1,x)-P_d_ass(1,x))/P_d_ass(1,x); % ����ѹ��������
        P_d_ass(1,x)=P_d(1,x); % �µĹ���ѹ������ֵ��Pa
    end
end

%% ��2����2*Nt-1����ʱ��ڵ���ز������㣨���ڣ�
for t=2:1:(2*Nt-1)
    P_d(t,Nx)=P_d(t-1,Nx); % ���ھ���ѹ������ֵ��Pa
    
    ERROR_PipInPressure=1; % �������ѹ��������
    COUNT_PipInPressure=0; % ��������
    while abs(ERROR_PipInPressure)>epsilon_e && COUNT_PipInPressure<epsilon_t
        COUNT_PipInPressure=COUNT_PipInPressure+1;
        
        % ��Nx���ռ�ڵ㴦��ز�������
        rho_d(t,Nx)=Density_TP(rho_d_0(t,Nx),T_0,P_0,T_d(t,Nx),P_d(t,Nx)); % ���ھ��������ܶȣ�kg/m^3
        Qv_d(t,Nx)=Qm_d_0(t,Nx)/rho_d(t,Nx); % ���ھ������������m^3/s
        V_d(t,Nx)=Qv_d(t,Nx)/A_d(Nx); % ���ھ����������٣�m/s
        [mu_d(t,Nx)]=Rheology_TP(mu_d_0(t,Nx),T_0,P_0,T_d(t,Nx),P_d(t,Nx)); % ���ھ�������ճ�ȣ�Pa*s��
        [Ff_d(t,Nx),flow_pattern_d(t,Nx)]=Friction_drillpipe(rho_d(t,Nx),V_d(t,Nx),mu_d(t,Nx),D_d_i(Nx)); % ���ھ��׵�λ����Ħ��ѹ����Pa/m����������̬������1��������2������3��
        Fa_d(t,Nx)=0; % ���ھ��׵�λ���ȼ��ٶ�ѹ����Pa/m
        
        % ��Nx-1��1���ռ�ڵ㴦��ز�������
        for x=Nx-1:-1:1
            P_d_ass(t,x)=P_d(t,x+1)-rho_d(t,x+1)*g*cosd(theta(x+1))*dx(x); % ����ѹ������ֵ��Pa
            
            ERROR_PipPressure=1; % ����ѹ��������
            COUNT_PipPressure=0; % ��������
            while abs(ERROR_PipPressure)>epsilon_e && COUNT_PipPressure<epsilon_t
                COUNT_PipPressure=COUNT_PipPressure+1;
                
                rho_d(t,x)=Density_TP(rho_d_0(t,x),T_0,P_0,T_d(t,x),P_d_ass(t,x)); % ���������ܶȣ�kg/m^3
                Qv_d(t,x)=Qm_d_0(t,x)/rho_d(t,x); % �������������m^3/s
                V_d(t,x)=Qv_d(t,x)/A_d(x); % �����������٣�m/s
                [mu_d(t,x)]=Rheology_TP(mu_d_0(t,x),T_0,P_0,T_d(t,x),P_d_ass(t,x)); % ��������ճ�ȣ�Pa*s��
                [Ff_d(t,x),flow_pattern_d(t,x)]=Friction_drillpipe(rho_d(t,x),V_d(t,x),mu_d(t,x),D_d_i(x)); % ���ڵ�λ����Ħ��ѹ����Pa/m����������̬������1��������2������3��
                Fa_d(t,x)=abs(0.5*(rho_d(t,x+1)*V_d(t,x+1)^2-rho_d(t,x)*V_d(t,x)^2)); % ���ڵ�λ���ȼ��ٶ�ѹ����Pa/m
                M1=(((rho_d(t,x)*V_d(t,x)+rho_d(t,x+1)*V_d(t,x+1))-(rho_d(t-1,x)*V_d(t-1,x)+rho_d(t-1,x+1)*V_d(t-1,x+1)))*dx(x))/(2*dt(t-1));
                M2=((rho_d(t,x+1)*V_d(t,x+1)^2+rho_d(t-1,x+1)*V_d(t-1,x+1)^2)-(rho_d(t,x)*V_d(t,x)^2+rho_d(t-1,x)*V_d(t-1,x)^2))/2;
                M3=-(((rho_d(t,x)*g*cosd(theta(x))-Ff_d(t,x)-Fa_d(t,x))+(rho_d(t,x+1)*g*cosd(theta(x+1))-Ff_d(t,x+1)-Fa_d(t,x+1))+(rho_d(t-1,x)*g*cosd(theta(x))-Ff_d(t-1,x)-Fa_d(t-1,x))+(rho_d(t-1,x+1)*g*cosd(theta(x+1))-Ff_d(t-1,x+1)-Fa_d(t-1,x+1)))*dx(x))/4;
                P_d(t,x)=P_d(t,x+1)+M1+M2+M3; % ����ѹ������ֵ��Pa
                
                ERROR_PipPressure=abs(P_d(t,x)-P_d_ass(t,x))/P_d_ass(t,x); % ����ѹ��������
                P_d_ass(t,x)=P_d(t,x); % �µĹ���ѹ������ֵ��Pa
            end
        end
        
        ERROR_PipInPressure=abs(P_d(t,1)-IP(t))/IP(t); % �������ѹ��������
        P_d(t,Nx)=P_d(t,Nx)-(P_d(t,1)-IP(t)); % �µľ���ѹ������ֵ��Pa
    end
end

%% ��1����2*Nt-1����ʱ��ڵ���ͷѹ������
for t=1:1:(2*Nt-1)
%     [delta_P_bit(t),V_nozzle(t)]=PressureDrop_Bit(C,d_nozzle,n_nozzle,Qv_d(t,Nx),rho_d(t,Nx)); % ��ͷѹ����Pa�����������٣�m/s��
    delta_P_bit(t)=0; % ��ͷѹ����Pa
end

%% ��1��ʱ��ڵ���ز������㣨���գ�
% ��Nx���ռ�ڵ㴦��ز�������
P_a(1,Nx)=P_d(1,Nx)-delta_P_bit(1); % ���վ���ѹ����Pa
rho_a(1,Nx)=Density_TP(rho_a_0(1,Nx),T_0,P_0,T_a(1,Nx),P_a(1,Nx)); % ���վ��������ܶȣ�kg/m^3
Qv_a(1,Nx)=Qm_a_0(1,Nx)/rho_a(1,Nx); % ���վ������������m^3/s
V_a(1,Nx)=Qv_a(1,Nx)/A_a(Nx); % ���վ����������٣�m/s
[mu_a(1,Nx)]=Rheology_TP(mu_a_0(1,Nx),T_0,P_0,T_a(1,Nx),P_a(1,Nx)); % ���վ�������ճ�ȣ�Pa*s��
[Ff_a(1,Nx),flow_pattern_a(1,Nx)]=Friction_annulus(rho_a(1,Nx),V_a(1,Nx),mu_a(1,Nx),D_w(Nx),D_d_o(Nx)); % ���վ��׵�λ����Ħ��ѹ����Pa/m����������̬������1��������2������3��
Fa_a(1,Nx)=0; % ���վ��׵�λ���ȼ��ٶ�ѹ����Pa/m

% ��Nx-1��1���ռ�ڵ㴦��ز�������
for x=Nx-1:-1:1
    P_a_ass(1,x)=P_a(1,x+1)-rho_a(1,x+1)*g*cosd(theta(x+1))*dx(x); % ����ѹ������ֵ��Pa
    
    ERROR_AnnPressure=1; % ����ѹ��������
    COUNT_AnnPressure=0; % ��������
    while abs(ERROR_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
        COUNT_AnnPressure=COUNT_AnnPressure+1;
        
        rho_a(1,x)=Density_TP(rho_a_0(1,x),T_0,P_0,T_a(1,x),P_a_ass(1,x)); % ���������ܶȣ�kg/m^3
        Qv_a(1,x)=Qm_a_0(1,x)/rho_a(1,x); % ���ճ������������m^3/s
        V_a(1,x)=Qv_a(1,x)/A_a(x); % �����������٣�m/s
        [mu_a(1,x)]=Rheology_TP(mu_a_0(1,x),T_0,P_0,T_a(1,x),P_a_ass(1,x)); % ��������ճ�ȣ�Pa*s��
        [Ff_a(1,x),flow_pattern_a(1,x)]=Friction_annulus(rho_a(1,x),V_a(1,x),mu_a(1,x),D_w(x),D_d_o(x)); % ���յ�λ����Ħ��ѹ����Pa/m����������̬������1��������2������3��
        Fa_a(1,x)=abs(0.5*(rho_a(1,x)*V_a(1,x)^2-rho_a(1,x+1)*V_a(1,x+1)^2)); % ���յ�λ���ȼ��ٶ�ѹ����Pa/m
        P_a(1,x)=-rho_a(1,x)*V_a(1,x)^2+rho_a(1,x+1)*V_a(1,x+1)^2+P_a(1,x+1)+(((-rho_a(1,x)*g*cosd(theta(x))-Ff_a(1,x)-Fa_a(1,x))+(-rho_a(1,x+1)*g*cosd(theta(x+1))-Ff_a(1,x+1)-Fa_a(1,x+1)))*dx(x))/2; % ����ѹ������ֵ��Pa
        
        ERROR_AnnPressure=abs(P_a(1,x)-P_a_ass(1,x))/P_a_ass(1,x); % ����ѹ��������
        P_a_ass(1,x)=P_a(1,x); % �µĻ���ѹ������ֵ��Pa
    end
end

%% ��2����2*Nt-1����ʱ��ڵ���ز������㣨���գ�
for t=2:1:(2*Nt-1)
    P_a_REAL(t,Nx)=P_d(t,Nx)-delta_P_bit(t); % ���վ���ѹ����Pa
    P_a(t,1)=P_a(t-1,1); % ���ճ���ѹ������ֵ��Pa
    
    ERROR_AnnInPressure=1; % ���վ���ѹ��������
    COUNT_AnnInPressure=0; % ��������
    while abs(ERROR_AnnInPressure)>epsilon_e && COUNT_AnnInPressure<epsilon_t
        COUNT_AnnInPressure=COUNT_AnnInPressure+1;
        
        % ��1���ռ�ڵ㴦��ز�������
        rho_a(t,1)=Density_TP(rho_a_0(t,1),T_0,P_0,T_a(t,1),P_a(t,1)); % ���ճ��������ܶȣ�kg/m^3
        Qv_a(t,1)=Qm_a_0(t,1)/rho_a(t,1); % ���ճ������������m^3/s
        V_a(t,1)=Qv_a(t,1)/A_a(1); % ���ճ����������٣�m/s
        [mu_a(t,1)]=Rheology_TP(mu_a_0(t,1),T_0,P_0,T_a(t,1),P_a(t,1)); % ���ճ�������ճ�ȣ�Pa*s��
        [Ff_a(t,1),flow_pattern_a(t,1)]=Friction_annulus(rho_a(t,1),V_a(t,1),mu_a(t,1),D_w(1),D_d_o(1)); % ���ճ��ڵ�λ����Ħ��ѹ����Pa/m����������̬������1��������2������3��
        Fa_a(t,1)=0; % ���ճ��ڵ�λ���ȼ��ٶ�ѹ����Pa/m
        
        % ��2��Nx���ռ�ڵ㴦��ز�������
        for x=2:1:Nx
            P_a_ass(t,x)=P_a(t,x-1)+rho_a(t,x-1)*g*cosd(theta(x-1))*dx(x-1); % ����ѹ������ֵ��Pa
            
            ERROR_AnnPressure=1; % ����ѹ��������
            COUNT_AnnPressure=0; % ��������
            while abs(ERROR_AnnPressure)>epsilon_e && COUNT_AnnPressure<epsilon_t
                COUNT_AnnPressure=COUNT_AnnPressure+1;
                
                rho_a(t,x)=Density_TP(rho_a_0(t,x),T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ���������ܶȣ�kg/m^3
                Qv_a(t,x)=Qm_a_0(t,x)/rho_a(t,x); % �������������m^3/s
                V_a(t,x)=Qv_a(t,x)/A_a(x); % �����������٣�m/s
                [mu_a(t,x)]=Rheology_TP(mu_a_0(t,x),T_0,P_0,T_a(t,x),P_a_ass(t,x)); % ��������ճ�ȣ�Pa*s��
                [Ff_a(t,x),flow_pattern_a(t,x)]=Friction_annulus(rho_a(t,x),V_a(t,x),mu_a(t,x),D_w(x),D_d_o(x)); % �������嵥λ����Ħ��ѹ����Pa/m����������̬������1��������2������3��
                Fa_a(t,x)=abs(0.5*(rho_a(t,x-1)*V_a(t,x-1)^2-rho_a(t,x)*V_a(t,x)^2)); % ���յ�λ���ȼ��ٶ�ѹ����Pa/m
                M1=(((rho_a(t,x-1)*V_a(t,x-1)+rho_a(t,x)*V_a(t,x))-(rho_a(t-1,x-1)*V_a(t-1,x-1)+rho_a(t-1,x)*V_a(t-1,x)))*dx(x-1))/(2*dt(t-1));
                M2=((rho_a(t-1,x-1)*V_a(t-1,x-1)^2+rho_a(t,x-1)*V_a(t,x-1)^2)-(rho_a(t-1,x)*V_a(t-1,x)^2+rho_a(t,x)*V_a(t,x)^2))/2;
                M3=-(((-rho_a(t,x)*g*cosd(theta(x))-Ff_a(t,x)-Fa_a(t,x))+(-rho_a(t,x-1)*g*cosd(theta(x-1))-Ff_a(t,x-1)-Fa_a(t,x-1))+(-rho_a(t-1,x)*g*cosd(theta(x))-Ff_a(t-1,x)-Fa_a(t-1,x))+(-rho_a(t-1,x-1)*g*cosd(theta(x-1))-Ff_a(t-1,x-1)-Fa_a(t-1,x-1)))*dx(x-1))/4;
                P_a(t,x)=P_a(t,x-1)+M1+M2+M3; % ����ѹ������ֵ��Pa
                
                ERROR_AnnPressure=abs(P_a(t,x)-P_a_ass(t,x))/P_a_ass(t,x); % ����ѹ��������
                P_a_ass(t,x)=P_a(t,x); % �µĻ���ѹ������ֵ��Pa
            end
        end
        
        ERROR_AnnInPressure=abs(P_a(t,Nx)-P_a_REAL(t,Nx))/P_a_REAL(t,Nx); % ���վ���ѹ��������
        P_a(t,1)=P_a(t,1)-(P_a(t,Nx)-P_a_REAL(t,Nx)); % �µĻ��ճ���ѹ������ֵ��Pa
    end
end

%% Ħ��ѹ������
for t=1:1:(2*Nt-1)
    FrictionalPressureDrop_annulus(t,1)=0; % FrictionalPressureDrop_annulusΪ����Ħ��ѹ����Pa
    FrictionalPressureDrop_drillpipe(t,1)=0; % FrictionalPressureDrop_drillpipeΪ���Ħ��ѹ����Pa
    for x=2:1:Nx
        FrictionalPressureDrop_annulus(t,x)=FrictionalPressureDrop_annulus(t,x-1)+Ff_a(t,x)*dx(x-1); % �ؾ���Ի���Ħ��ѹ���ۼ���ͣ�Pa
        FrictionalPressureDrop_drillpipe(t,x)=FrictionalPressureDrop_drillpipe(t,x-1)+Ff_d(t,x)*dx(x-1); % �ؾ�������Ħ��ѹ���ۼ���ͣ�Pa
    end
end

%% �������
fprintf('��ʼʱ�����ѹ����MPa����%8.5f\n',P_d(1,1)/(10^6)); % ��ʼʱ�����ѹ����MPa
fprintf('��ʼʱ�̻��վ���ѹ����MPa����%8.5f\n',P_a(1,Nx)/(10^6)); % ��ʼʱ�̻��վ���ѹ����MPa
fprintf('��ʼʱ������Ħ��ѹ����MPa����%8.5f\n',FrictionalPressureDrop_drillpipe(1,Nx)/(10^6)); % ��ʼʱ������Ħ��ѹ����MPa
fprintf('��ʼʱ�̻���Ħ��ѹ����MPa����%8.5f\n',FrictionalPressureDrop_annulus(1,Nx)/(10^6)); % ��ʼʱ�̻���Ħ��ѹ����MPa
fprintf('��ʼʱ����ͷѹ����MPa����%8.5f\n',delta_P_bit(1)/(10^6)); % ��ʼʱ����ͷѹ����MPa

%% ����ͼ��
Time(1)=0;
for t=2:1:Nt
    Time(t)=Time(t-1)+dt_d(t-1);
end
for t=Nt+1:1:(2*Nt-1)
    Time(t)=Time(t-1)+dt_a(2*Nt-t);
end

figure(1); % ���չ���ѹ��VS����VSʱ��
for t=1:1:(2*Nt-1)
    plot(P_a(t,:)/10^6,Depth,'--',P_d(t,:)/10^6,Depth); % ���ա�����ѹ����MPa��
    xlabel('ѹ����MPa��','FontName','����','FontSize',16);
    ylabel('���m��','FontName','����','FontSize',16);
    title([num2str(Time(t)/60),' min']);
    legend('����ѹ��','����ѹ��');
    set(gca,'fontsize',16);
    set(gca,'YDir','reverse');
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    pause(0.001);
    
    % ����gif��ͼ
    drawnow;
    F1=getframe(gcf);
    I1=frame2im(F1);
    [I1,map1]=rgb2ind(I1,256);
    if t == 1
        imwrite(I1,map1,'���չ���ѹ��VS����VSʱ��.gif','gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I1,map1,'���չ���ѹ��VS����VSʱ��.gif','gif','WriteMode','append','DelayTime',0.2);
    end
end

figure(2); % ����ѹ��VSʱ��
for t=1:1:(2*Nt-1)
    scatter(Time(t)/60,P_d(t,1)/10^6); % ����ѹ����MPa��
    hold on;
    xlabel('ʱ�䣨min��','FontName','����','FontSize',16);
    ylabel('����ѹ����MPa��','FontName','����','FontSize',16);
    set(gca,'fontsize',16);
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    pause(0.001);
    
    % ����gif��ͼ
    drawnow;
    F2=getframe(gcf);
    I2=frame2im(F2);
    [I2,map2]=rgb2ind(I2,256);
    if t == 1
        imwrite(I2,map2,'����ѹ��VSʱ��.gif','gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I2,map2,'����ѹ��VSʱ��.gif','gif','WriteMode','append','DelayTime',0.2);
    end
end

figure(3); % ����ѹ��VSʱ��
for t=1:1:(2*Nt-1)
    scatter(Time(t)/60,P_a(t,1)/10^6); % ����ѹ����MPa��
    hold on;
    xlabel('ʱ�䣨min��','FontName','����','FontSize',16);
    ylabel('����ѹ����MPa��','FontName','����','FontSize',16);
    set(gca,'fontsize',16);
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    pause(0.001);
    
    % ����gif��ͼ
    drawnow;
    F3=getframe(gcf);
    I3=frame2im(F3);
    [I3,map3]=rgb2ind(I3,256);
    if t == 1
        imwrite(I3,map3,'����ѹ��VSʱ��.gif','gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I3,map3,'����ѹ��VSʱ��.gif','gif','WriteMode','append','DelayTime',0.2);
    end
end

figure(4); % ��Ͳ�¶�VS����VSʱ��
for t = 1:1:length(Time_T)
    plot(T_d_T(:,t),Depth_T,T_a_T(:,t),Depth_T,'--'); % ��Ͳ�¶ȣ���
    xlabel('�¶ȣ��棩','FontName','����','FontSize',10);
    ylabel('���m��','FontName','����','FontSize',10);
    title([num2str(Time_T(t)/60),' min']);
    set(gca,'YDir','reverse');
    legend('������¶�','�������¶�','Location','Best');
    set(gca,'FontSize',14);
    box on; % ��ʾ������ı߿�
    grid on; % ��ʾ���������������
    grid minor; % ��ʾ������Ĵ�������
    
    pause(0.001);
    
    % ����gif��ͼ
    drawnow;
    F4=getframe(gcf);
    I4=frame2im(F4);
    [I4,map4]=rgb2ind(I4,256);
    if t == 1
        imwrite(I4,map4,'��Ͳ�¶�VS����VSʱ��.gif','gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(I4,map4,'��Ͳ�¶�VS����VSʱ��.gif','gif','WriteMode','append','DelayTime',0.2);
    end
end