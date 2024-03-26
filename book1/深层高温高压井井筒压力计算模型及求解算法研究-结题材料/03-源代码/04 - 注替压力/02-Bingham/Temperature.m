function [T_Annulus,T_Drillpipe,Depth,Time]=Temperature(welldepth,rho_0,mu_0,Qv_0)
% ���������welldepthΪ���m����rho_0Ϊ�꾮Һ�ܶȣ�kg/m^3����mu_0Ϊ�꾮Һճ�ȣ�Pa��s����Qv_0Ϊ�꾮Һ������m^3/s��
% ���������T_AnnulusΪ�����꾮Һ�¶ȣ��棩��T_DrillpipeΪ�����꾮Һ�¶ȣ��棩��DepthΪ�ռ�ڵ㣨m����TimeΪʱ��ڵ㣨s��

delta_t=3600;           % ʱ�䲽����s
timesteps=120;          % ѭ������������ģ��ʱ��Ϊdelta_t��timesteps

dz=40;                  % ��dzΪ��Ҫ������Ԫ���ȣ�m

casing_number=3;        % �׹ܲ���������㡢�����������׹ܲ�����������β�ܣ�

Depth_weiguan=2505;     % β��������ȣ�m

% weiguan=0;              % ������β�ܣ��򻷿��¶�ֻ��һ�������
weiguan=1;              % ����β�ܣ����轫�����¶ȳ���Ϊ���������

Rheology_type=1;        % 1Ϊţ�����壬2Ϊ��������
% Rheology_type=2;        % 1Ϊţ�����壬2Ϊ��������

%% ����ṹ�������񻮷�
nobox_r=200;    % ��������Ԫ��

for i=1:casing_number   % �׹���ţ���Ŵ����������ӣ���casing_number=1Ϊ�����׹�
    Depth_casing(i)=welldepth-4000/4*i;     % �׹����m
    if i==casing_number                     % ������׹ܣ�����׹ܣ�
        Depth_fluid(i)=0;                   % ����׹�ˮ�෵�m��һ�㷵�����棩
    else
        Depth_fluid(i)=Depth_casing(i)-150; % �Ǳ���׹�ˮ�෵�m��һ�㷵���ܷ������50~200�ף�
    end
end

%% �Ƿ����β�ܣ�����β�ܺͲ�����β�ܵ�ϵ�����������ܴ����ʹ�ò�ͬ��ģ��
if weiguan==0 % ������β��
    section_number=casing_number+1;             % ������=�׹ܲ���+1���ؾ���ֶΣ�ÿһ����dzΪ��Ҫ����Ԫ���Ȼ���������ŴӾ����򾮿�����
    if casing_number==0                         % ����׹ܲ���=0�������۾����򾮶���Ϊ1
        section_length1(1)=welldepth;           % ��һ���εĶγ�Ϊ���m
    else
        for i=1:section_number
            if i==1
                section_length1(i)=welldepth-Depth_casing(i);
            elseif i==section_number
                section_length1(i)=Depth_casing(i-1);
            else
                section_length1(i)=Depth_casing(i-1)-Depth_casing(i);
            end
        end
    end
else % ����β��
    section_number=casing_number+2;             % ������=�׹ܲ���+2���ؾ���ֶΣ�ÿһ����dzΪ��Ҫ����Ԫ���Ȼ���������ŴӾ����򾮿�����
    if casing_number==0
        fprintf('β�ܴ��ڵ�ʱ���������׹ܣ����׹���>=1')
    elseif casing_number==1
        section_length1(1)=welldepth-Depth_casing(1);
        section_length1(2)=Depth_casing(1)-Depth_weiguan;
        section_length1(3)=Depth_weiguan;
    elseif casing_number==2
        section_length1(1)=welldepth-Depth_casing(1);
        section_length1(2)=Depth_casing(1)-Depth_weiguan;
        section_length1(3)=Depth_weiguan-Depth_casing(2);
        section_length1(4)=Depth_casing(2);
    else
        for i=1:section_number
            if i==1
                section_length1(1)=welldepth-Depth_casing(1);
            elseif i==2
                section_length1(2)=Depth_casing(1)-Depth_weiguan;
            elseif i==3
                section_length1(3)=Depth_weiguan-Depth_casing(2);
            elseif i==section_number
                section_length1(i)=Depth_casing(i-2);
            else
                section_length1(i)=Depth_casing(i-2)-Depth_casing(i-1);
            end
        end      
    end
end

section_length1(find(section_length1==0))=[];   % ȥ�������е�0Ԫ�أ����ڶγ�Ϊ0�Ĳ�����Ϊһ�λ�������
section_length=section_length1(end:-1:1);       % ���γ��������򣬼���ŴӾ����򾮵�����
section_number=length(section_length);          % ���ڿ���ȥ���������е�0Ԫ�أ������¼���γ���

for i=1:section_number
    n(i)=floor(section_length(i)/dz);           % ����ÿһ���б�dz������
    m(i)=mod(section_length(i),dz);             % ����ÿһ���б�dz�����������
end

delta_z=[];                                     % �ܵ�����Ԫ��������
for i=1:section_number
    int1=[];                                    % ���ε�����Ԫ��������
    int1(1:n(i))=dz;                            % �ö�n(i)��dz������
    if m(i)==0
        int1=int1;                              % �������Ϊ0�����������������������0�������뵽����Ԫ��������
    else
        int1=[int1,m(i)];                       % ���������Ϊ0�������û�����������������뵽����Ԫ��������
    end
    delta_z=[delta_z,int1];                     % ����һ�ε�����Ԫ������������һ�ε�����Ԫ����ƴ����һ��
end

nobox_z=length(delta_z)+1;                      % ����߽���=������+1

Depth=zeros(nobox_z,1);                         % ����������m
Depth(1)=0;                                     % ������ȣ�m
for i=2:nobox_z
    Depth(i)=Depth(i-1)+delta_z(i-1);
end

Weiguan_i=find(Depth==Depth_weiguan);           % �ж�β�ܶ�Ӧ�Ľڵ���

%% ���������ֵ
[rou_m,Cpm,lambdam]=deal(rho_0,1600,1.75);       % �꾮Һ���ܶȣ�kg/m^3���������ݣ�J/(kg���棩��������ϵ����W/(m���棩��
[rou_s,Cps,lambdas]=deal(8000,400,43.75);       % �ֲģ���˻��׹ܣ����ܶȣ�kg/m^3���������ݣ�J/(kg���棩��������ϵ����W/(m���棩��
[rou_f,Cpf,lambdaf]=deal(2640,800,2.25);        % �ز���ܶȣ�kg/m^3���������ݣ�J/(kg���棩��������ϵ����W/(m���棩��
[rou_c,Cpc,lambdac]=deal(2140,2000,0.7);        % ˮ����ܶȣ�kg/m^3���������ݣ�J/(kg���棩��������ϵ����W/(m���棩��

Tm0=20;         % �꾮Һע���¶ȣ���
Tf0=15;         % �ر��¶ȣ���
Tfg=0.03;     % �����ݶȣ���/m

r=zeros(nobox_z,nobox_r+1);     % �뾶�����ؾ���;�����������

if weiguan==0   % ������β��
    r(:,1)=0.10053/2;                % ����ڱڰ뾶��m
    r(:,2)=0.1143/2;                % �����ڰ뾶��m
    r(:,3)=0.14/2;                % ��Ͳ�ڱڰ뾶��m
    r(:,4)=0.1851/2;                % �����4����ڰ뾶��m
    r(:,5)=0.2051/2;                % �����5����ڰ뾶��m
    r(:,6)=0.2251/2;                % �����6����ڰ뾶��m
    r(:,7)=0.2451/2;                % �����7����ڰ뾶��m
    r(:,8)=0.2651/2;                % �����8����ڰ뾶��m
    r(:,9)=0.2851/2;                % �����9����ڰ뾶��m

    for i=10:nobox_r+1
        delta_r(:,i)=0.02;              % �����i�����i-1��֮��İ뾶�m
        r(:,i)=r(:,i-1)+delta_r(:,i);   % �����i����ڰ뾶��m
    end
                 
    A_annulus=pi*(r(:,3).^2-r(:,2).^2);     % �����ڽ������m^2
else    % ����β��
    r(:,1)=0.0943/2;        % ����ڱڰ뾶��m
    r(:,2)=0.1143/2;        % �����ڰ뾶��m
    r(:,3)=0.1451/2;        % β���ڱڰ뾶��m
    r(:,4)=0.1551/2;        % β��ˮ�໷�ڱڰ뾶��m
    r(:,5)=0.1651/2;        % ���հ뾶��m
    r(:,6)=0.1851/2;        % �Ͳ��׹���ڰ뾶��m
    r(:,7)=0.2051/2;        % �����׹��ڱڰ뾶��m
    r(:,8)=0.2251/2;        % �����׹���ڰ뾶��m
    r(:,9)=0.2451/2;        % ����׹��ڱڰ뾶��m
    r(:,10)=0.2651/2;       % ����׹���ڰ뾶��m
    r(:,11)=0.2851/2;       % ���ˮ�໷��ڰ뾶��m
    
    for i=12:nobox_r+1
        delta_r(:,i)=0.02;  % �����i�����i-1��֮��İ뾶�m
        r(:,i)=r(:,i-1)+delta_r(:,i);   % �����i����ڰ뾶��m
    end
    
    A_annulus=zeros(nobox_z,1);     % ���㻷���ڽ������m^2
    for j=1:nobox_z
        if Depth(j)<Depth_weiguan
            A_annulus(j)=pi*(r(j,5)^2-r(j,2)^2);
        else
            A_annulus(j)=pi*(r(j,3)^2-r(j,2)^2);
        end
        
    end
end

A_drillpipe=pi*r(:,1).^2;       % ��������ڽ������m^2

q=Qv_0;                      % �꾮Һ������m^3/s

V_drillpipe=q./A_drillpipe;             % ������꾮Һ���٣�m/s
V_annulus=q./A_annulus;                 % �������꾮Һ���٣�m/s

if Rheology_type==1
    miu_m=mu_0;     % �꾮Һ�ȣ�Pa��s
    miu_drillpipe(1:nobox_z,1)=miu_m;       % ������꾮Һ�ȣ�Pa��s
    miu_annulus(1:nobox_z,1)=miu_m;         % �������꾮Һ�ȣ�Pa��s
else
    K_power=0.8;    % ����������ϵ��,Pa��s^n
    n_power=0.6;    % ������������ָ��
    K_drillpipe(1:nobox_z,1)=K_power;
    n_drillpipe(1:nobox_z,1)=n_power;
    miu_drillpipe(1:nobox_z,1)=1/2*1022^n_power*K_power/0.511/1000;
    miu_annulus=miu_drillpipe;
end

Rho=zeros(nobox_z,nobox_r);             % �ܶȾ����ؾ���;�����
Cp=zeros(nobox_z,nobox_r);              % �����ݾ����ؾ���;�����
lambda=zeros(nobox_z,nobox_r+1);        % ����ϵ�������ؾ���;�����

Rho(:,1)=rou_m;                         % ��һ����ܶ�
Cp(:,1)=Cpm;                            % ��һ��ı�����
lambda(:,1)=lambdam;                    % ��һ��ĵ���ϵ��

Rho(:,2)=rou_s;                         % �ڶ�����ܶ�
Cp(:,2)=Cps;                            % �ڶ���ı�����
lambda(:,2)=lambdas;                    % �ڶ���ĵ���ϵ��

Rho(:,3)=rou_m;                         % ��������ܶ�
Cp(:,3)=Cpm;                            % ������ı�����
lambda(:,3)=lambdam;                    % ������ĵ���ϵ��

lambda(:,nobox_r+1)=lambdaf;        % �ڼ������һ��ز���¶�ʱ����Ҫ�õ�nobox_r+1�ĵ���ϵ��

if weiguan==0           % ������β��

    if casing_number==0                     % ���۾������4�㼰֮��Ĳ�����Ϊ�ز�
        Rho(:,4:nobox_r)=rou_f;             % 4�㼰֮��ĵز��ܶ�
        Cp(:,4:nobox_r)=Cpf;                % 4�㼰֮��ĵز������
        lambda(:,4:nobox_r)=lambdaf;        % 4�㼰֮��ĵز㵼��ϵ��
    else
        for i=1:casing_number
            for j=1:nobox_z
                if Depth(j)<Depth_casing(i) % С���׹������ֵΪ�׹ܲ���
                    Rho(j,3+2*i-1)=rou_s;       % �ܶ�
                    Cp(j,3+2*i-1)=Cps;          % ������
                    lambda(j,3+2*i-1)=lambdas;  % ����ϵ��
                elseif Depth(j)>=Depth_casing(i) && Depth(j)<=welldepth  % �����׹������ֵΪ�ز����
                    Rho(j,3+2*i-1)=rou_f;       % �ܶ�
                    Cp(j,3+2*i-1)=Cpf;          % ������
                    lambda(j,3+2*i-1)=lambdaf;  % ����ϵ��
                end

                if Depth(j)<Depth_fluid(i)  % С��ˮ�໷�����ֵΪ�꾮Һ����
                    Rho(j,3+2*i)=rou_m;         % �ܶ�
                    Cp(j,3+2*i)=Cpm;            % ������
                    lambda(j,3+2*i)=lambdam;    % ����ϵ��
                elseif Depth(j)>=Depth_fluid(i) & Depth(j)<Depth_casing(i)  % ����ˮ�෵����С���׹������ֵΪˮ�����
                    Rho(j,3+2*i)=rou_c;         % �ܶ�
                    Cp(j,3+2*i)=Cpc;            % ������
                    lambda(j,3+2*i)=lambdac;    % ����ϵ��
                elseif Depth(j)>=Depth_casing(i) & Depth(j)<=welldepth      % �����׹������ֵΪ�ز����
                    Rho(j,3+2*i)=rou_f;         % �ܶ�
                    Cp(j,3+2*i)=Cpf;            % ������
                    lambda(j,3+2*i)=lambdaf;    % ����ϵ��
                end
            end
        end
    end
    for i=3+2*casing_number+1:nobox_r   % ��ȥ�׹ܺ�ˮ�࣬ʣ��ľ����ȫ����ֵΪ�ز����
        for j=1:nobox_z
            Rho(j,i)=rou_f;             % �ܶ�
            Cp(j,i)=Cpf;                % ������
            lambda(j,i)=lambdaf;        % ����ϵ��
        end
    end

    if Rheology_type==1
        Re_drillpipe=Rho(:,1).*V_drillpipe.*(2*r(:,1))./miu_m;                  % �������ŵ��
        Re_annulus=Rho(:,3).*V_annulus.*(2*(r(:,3)-r(:,2)))/miu_m;              % ��������ŵ��
        h_drillpipe=0.027.*lambda(:,1).*Re_drillpipe.^0.8.*(miu_drillpipe.*Cp(:,1)./lambda(:,1)).^0.4./(2*r(:,1))/5;        % �꾮Һ������ڱڶ�������ϵ����W/(m^2���棩   
        h_annulus=0.027*lambda(:,3).*Re_annulus.^0.8.*(miu_drillpipe.*Cp(:,3)./lambda(:,3)).^0.4./(2*(r(:,3)-r(:,2)))/5;	% �꾮Һ�������ڶ�������ϵ����W/(m^2���棩
        h1=h_drillpipe;                                                                                                     % �꾮Һ������ڱڶ�������ϵ����W/(m^2���棩    
        h2=h_annulus;                                                                                                       % �꾮Һ�������ڶ�������ϵ����W/(m^2���棩
        h3=h_annulus;                                                                                                       % �꾮Һ�뾮Ͳ�ڱڶ�������ϵ����W/(m^2���棩
    elseif Rheology_type==2
        Re_drillpipe=Rho(:,1).*V_drillpipe.^(2-n_power).*(2*r(:,1)).^n_power./(8^(n_power-1).*K_power.*((1+3*n_power)./(4*n_power)).^n_power);                  % �������ŵ��
        delta_correct=(3*n_power+1)/(4*n_power);
        h_drillpipe=0.091*Re_drillpipe.^0.87.*(miu_drillpipe.*Cp(:,1)./lambda(:,1)).^(1/3)*delta_correct./(2*r(:,1))/5;
        Re_annulus=zeros(nobox_z,1);
        h_annulus=zeros(nobox_z,1);
        for j=1:nobox_z
            if Depth(j)<Depth_weiguan
                Re_annulus(j)=Rho(j,3)*V_annulus(j)^(2-n_power)*(2*(r(j,5)-r(j,2)))^n_power/(8^(n_power-1)*K_power*((1+3*n_power)/(4*n_power))^n_power);              % ��������ŵ��
                h_annulus(j)=0.091*Re_annulus(j)^0.87*(miu_annulus(j)*Cp(j,3)./lambda(j,3))^(1/3)*delta_correct/(2*(r(j,5)-r(j,2)))/5;  % �꾮Һ�������ڶ�������ϵ����W/(m^2���棩
            else
                Re_annulus(j)=Rho(j,3)*V_annulus(j)^(2-n_power)*(2*(r(j,3)-r(j,2)))^n_power/(8^(n_power-1)*K_power*((1+3*n_power)/(4*n_power))^n_power);              % ��������ŵ��
                h_annulus(j)=0.091*Re_annulus(j)^0.87*(miu_annulus(j)*Cp(j,3)./lambda(j,3))^(1/3)*delta_correct/(2*(r(j,3)-r(j,2)))/5;	% �꾮Һ�������ڶ�������ϵ����W/(m^2���棩
            end
        end
    end
    
    h1=h_drillpipe;                                                                                                     % �꾮Һ������ڱڶ�������ϵ����W/(m^2���棩    
    h2=h_annulus;                                                                                                       % �꾮Һ�������ڶ�������ϵ����W/(m^2���棩
    h3=h_annulus;                                                                                                       % �꾮Һ�뾮Ͳ�ڱڶ�������ϵ����W/(m^2���棩
    
    %% ϵ���������
    A=zeros(nobox_z,nobox_r);   % ϵ������A
    B=zeros(nobox_z,nobox_r);   % ϵ������B
    C=zeros(nobox_z,nobox_r);   % ϵ������C
    D=zeros(nobox_z,nobox_r);   % ϵ������D
    E=zeros(nobox_z,nobox_r);   % ϵ������E
    F=zeros(nobox_z,nobox_r);   % ϵ������F

    for j=1:nobox_z
        if j==1     % ��������ϵ���ļ���
            A(j,1)=-Rho(j,1)*q*Cp(j,1)/delta_z(j)-2*pi*r(j,1)*h1(j,1)-Rho(j,1)*Cp(j,1)*pi*r(j,1)^2/delta_t;             % ����������꾮Һ�¶ȳ������ϵ��A
            B(j,1)=Rho(j,1)*q*Cp(j,1)/delta_z(j);                                                                       % ����������꾮Һ�¶ȳ������ϵ��B                                        
            C(j,1)=0;                                                                                                   % ����������꾮Һ�¶ȳ������ϵ��C
            D(j,1)=0;                                                                                                   % ����������꾮Һ�¶ȳ������ϵ��D
            E(j,1)=2*pi*r(j,1)*h1(j,1);                                                                                 % ����������꾮Һ�¶ȳ������ϵ��E
            F(j,1)=-Rho(j,1)*Cp(j,1)*pi*r(j,1)^2/delta_t;                                                               % ����������꾮Һ�¶ȳ������ϵ��F

            A(j,2)=-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j)+delta_z(j))*delta_z(j))-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j)+delta_z(j))*delta_z(j))-2*pi*r(j,1)*h1(j,1)-2*pi*r(j,2)*h2(j,1)-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;    % ������˱����¶ȳ������ϵ��A
            B(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j)+delta_z(j))*delta_z(j));                           % ������˱����¶ȳ������ϵ��B
            C(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j)+delta_z(j))*delta_z(j));                           % ������˱����¶ȳ������ϵ��C
            D(j,2)=2*pi*r(j,1)*h1(j,1);                                                                                 % ������˱����¶ȳ������ϵ��D
            E(j,2)=2*pi*r(j,2)*h2(j,1);                                                                                 % ������˱����¶ȳ������ϵ��E
            F(j,2)=-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;                                                    % ������˱����¶ȳ������ϵ��F

            A(j,3)=-Rho(j,3)*q*Cp(j,3)/delta_z(j)-2*pi*r(j,2)*h2(j,1)-2*pi*r(j,3)*h3(j,1)-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;  % ���㻷�����꾮Һ�¶ȳ������ϵ��A
            B(j,3)=0;                                                                                                   % ���㻷�����꾮Һ�¶ȳ������ϵ��B
            C(j,3)=Rho(j,3)*q*Cp(j,3)/delta_z(j);                                                                       % ���㻷�����꾮Һ�¶ȳ������ϵ��C
            D(j,3)=2*pi*r(j,2)*h2(j,1);                                                                                 % ���㻷�����꾮Һ�¶ȳ������ϵ��D
            E(j,3)=2*pi*r(j,3)*h3(j,1);                                                                                 % ���㻷�����꾮Һ�¶ȳ������ϵ��E
            F(j,3)=-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;                                                    % ���㻷�����꾮Һ�¶ȳ������ϵ��F

            A(j,4)=-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j)+delta_z(j))*delta_z(j))-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j)+delta_z(j))*delta_z(j))-2*pi*r(j,3)*h3(j,1)-pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)))-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;  % ������Ĳ��¶ȳ������ϵ��A
            B(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j)+delta_z(j))*delta_z(j));                           % ������Ĳ��¶ȳ������ϵ��B
            C(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j)+delta_z(j))*delta_z(j));                           % ������Ĳ��¶ȳ������ϵ��C
            D(j,4)=2*pi*r(j,3)*h3(j,1);                                                                                 % ������Ĳ��¶ȳ������ϵ��D
            E(j,4)=pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)));   % ������Ĳ��¶ȳ������ϵ��E
            F(j,4)=-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;                                                    % ������Ĳ��¶ȳ������ϵ��F
        else
            A(j,1)=-Rho(j,1)*q*Cp(j,1)/delta_z(j-1)-2*pi*r(j,1)*h1(j,1)-Rho(j,1)*Cp(j,1)*pi*r(j,1)^2/delta_t;           % ����������꾮Һ�¶ȳ������ϵ��A
            B(j,1)=Rho(j,1)*q*Cp(j,1)/delta_z(j-1);                                                                     % ����������꾮Һ�¶ȳ������ϵ��B
            C(j,1)=0;                                                                                                   % ����������꾮Һ�¶ȳ������ϵ��C
            D(j,1)=0;                                                                                                   % ����������꾮Һ�¶ȳ������ϵ��D
            E(j,1)=2*pi*r(j,1)*h1(j,1);                                                                                 % ����������꾮Һ�¶ȳ������ϵ��E
            F(j,1)=-Rho(j,1)*Cp(j,1)*pi*r(j,1)^2/delta_t;                                                               % ����������꾮Һ�¶ȳ������ϵ��F

            if j==nobox_z   % ��������ϵ���ļ���
                A(j,2)=-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-2*pi*r(j,1)*h1(j,1)-2*pi*r(j,2)*h2(j,1)-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;    % ������˱����¶ȳ������ϵ��A
                B(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));                 % ������˱����¶ȳ������ϵ��B
                C(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));                 % ������˱����¶ȳ������ϵ��C
                D(j,2)=2*pi*r(j,1)*h1(j,1);                                                                             % ������˱����¶ȳ������ϵ��D
                E(j,2)=2*pi*r(j,2)*h2(j,1);                                                                             % ������˱����¶ȳ������ϵ��E
                F(j,2)=-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;                                                % ������˱����¶ȳ������ϵ��F

                A(j,3)=-Rho(j,3)*q*Cp(j,3)/delta_z(j-1)-2*pi*r(j,2)*h2(j,1)-2*pi*r(j,3)*h3(j,1)-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;    % ���㻷�����꾮Һ�¶ȳ������ϵ��A
                B(j,3)=0;                                                                                               % ���㻷�����꾮Һ�¶ȳ������ϵ��B
                C(j,3)=Rho(j,3)*q*Cp(j,3)/delta_z(j-1);                                                                 % ���㻷�����꾮Һ�¶ȳ������ϵ��C
                D(j,3)=2*pi*r(j,2)*h2(j,1);                                                                             % ���㻷�����꾮Һ�¶ȳ������ϵ��D
                E(j,3)=2*pi*r(j,3)*h3(j,1);                                                                             % ���㻷�����꾮Һ�¶ȳ������ϵ��E
                F(j,3)=-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;                                                % ���㻷�����꾮Һ�¶ȳ������ϵ��F

                A(j,4)=-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-2*pi*r(j,3)*h3(j,1)-pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)))-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;  % ������Ĳ��¶ȳ������ϵ��A
                B(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));                 % ������Ĳ��¶ȳ������ϵ��B
                C(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));                 % ������Ĳ��¶ȳ������ϵ��C
                D(j,4)=2*pi*r(j,3)*h3(j,1);                                                                             % ������Ĳ��¶ȳ������ϵ��D
                E(j,4)=pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)));   % ������Ĳ��¶ȳ������ϵ��E
                F(j,4)=-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;                                                % ������Ĳ��¶ȳ������ϵ��F
            else    % �м�����ϵ���ļ���
                A(j,2)=-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j))-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j-1))-2*pi*r(j,1)*h1(j,1)-2*pi*r(j,2)*h2(j,1)-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;
                B(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j-1));
                C(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j));
                D(j,2)=2*pi*r(j,1)*h1(j,1);
                E(j,2)=2*pi*r(j,2)*h2(j,1);
                F(j,2)=-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;

                A(j,3)=-Rho(j,3)*q*Cp(j,3)/delta_z(j)-2*pi*r(j,2)*h2(j,1)-2*pi*r(j,3)*h3(j,1)-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;
                B(j,3)=0;
                C(j,3)=Rho(j,3)*q*Cp(j,3)/delta_z(j);
                D(j,3)=2*pi*r(j,2)*h2(j,1);
                E(j,3)=2*pi*r(j,3)*h3(j,1);
                F(j,3)=-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;

                A(j,4)=-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j))-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j-1))-2*pi*r(j,3)*h3(j,1)-pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)))-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;
                B(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j-1));
                C(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j));
                D(j,4)=2*pi*r(j,3)*h3(j,1);
                E(j,4)=pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)));
                F(j,4)=-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;
            end
        end
    end

    for i=5:nobox_r     % �����5�㼰֮������ϵ������
        for j=1:nobox_z
            if j==1     % ��������ϵ���ļ���
                A(j,i)=-log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)))-log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)))-2*lambda(j,i)/((delta_z(j)+delta_z(j))*delta_z(j))-2*lambda(j,i)/((delta_z(j)+delta_z(j))*delta_z(j))-Rho(j,i)*Cp(j,i)/delta_t;
                B(j,i)=2*lambda(j,i)/((delta_z(j)+delta_z(j))*delta_z(j));
                C(j,i)=2*lambda(j,i)/((delta_z(j)+delta_z(j))*delta_z(j));
                D(j,i)=log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)));
                E(j,i)=log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)));
                F(j,i)=-Rho(j,i)*Cp(j,i)/delta_t;
            elseif j==nobox_z       % ��������ϵ���ļ���
                A(j,i)=-log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)))-log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)))-2*lambda(j,i)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-2*lambda(j,i)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-Rho(j,i)*Cp(j,i)/delta_t;
                B(j,i)=2*lambda(j,i)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));
                C(j,i)=2*lambda(j,i)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));
                D(j,i)=log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)));
                E(j,i)=log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)));
                F(j,i)=-Rho(j,i)*Cp(j,i)/delta_t;
            else        % �м�����ϵ���ļ���
                A(j,i)=-log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)))-log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)))-2*lambda(j,i)/((delta_z(j-1)+delta_z(j))*delta_z(j))-2*lambda(j,i)/((delta_z(j-1)+delta_z(j))*delta_z(j-1))-Rho(j,i)*Cp(j,i)/delta_t;
                B(j,i)=2*lambda(j,i)/((delta_z(j-1)+delta_z(j))*delta_z(j-1));
                C(j,i)=2*lambda(j,i)/((delta_z(j-1)+delta_z(j))*delta_z(j));
                D(j,i)=log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)));
                E(j,i)=log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)));
                F(j,i)=-Rho(j,i)*Cp(j,i)/delta_t;
            end
        end
    end

    A99=zeros(nobox_r*nobox_z,1);   % ϵ������A�洢������A99��
    B99=zeros(nobox_r*nobox_z,1);   % ϵ������B�洢������B99��
    C99=zeros(nobox_r*nobox_z,1);   % ϵ������C�洢������C99��
    D99=zeros(nobox_r*nobox_z,1);   % ϵ������D�洢������D99��
    E99=zeros(nobox_r*nobox_z,1);   % ϵ������E�洢������E99��
    F99=zeros(nobox_r*nobox_z,1);   % ϵ������F�洢������F99��

    % ϵ������A�洢������A99��
    A99(1:nobox_z)=A(:,1);
    A99(nobox_z+1)=A(1,2)+B(1,2);
    A99(nobox_z+2:2*nobox_z-1)=A(2:nobox_z-1,2);
    A99(2*nobox_z)=A(nobox_z,2)+C(nobox_z,2);
    A99(2*nobox_z+1:3*nobox_z-1)=A(1:nobox_z-1,3);
    A99(3*nobox_z)=A(nobox_z,3)+C(nobox_z,3);
    for i=4:nobox_r
        A99((i-1)*nobox_z+1)=A(1,i)+B(1,i);
        A99((i-1)*nobox_z+2:i*nobox_z-1)=A(2:nobox_z-1,i);
        A99(i*nobox_z)=A(nobox_z,i)+C(nobox_z,i);
    end

    % ϵ������B��C��D��E��F�洢������B99��C99��D99��E99��F99��
    for i=1:nobox_r
        B99((i-1)*nobox_z+1:i*nobox_z)=B(:,i);
        C99((i-1)*nobox_z+1:i*nobox_z)=C(:,i);
        D99((i-1)*nobox_z+1:i*nobox_z)=D(:,i);
        E99((i-1)*nobox_z+1:i*nobox_z)=E(:,i);
        F99((i-1)*nobox_z+1:i*nobox_z)=F(:,i);
    end

    %% ���¶ȳ���������ֵ
    X=zeros(nobox_r*nobox_z,1);     % �¶ȳ�����
    % ������꾮Һ����˱ڡ��������꾮Һ�¶�
    X(1)=Tm0;           % �����¶�Ϊ�꾮Һע���¶�
    X(1+nobox_z)=Tm0;   % �����¶�Ϊ�꾮Һע���¶�
    X(1+2*nobox_z)=Tm0; % �����¶�Ϊ�꾮Һע���¶�
    for j=2:nobox_z
        X(j)=X(j-1)+delta_z(j-1)*Tfg;                       % ������꾮Һ�¶��Ե����ݶ�����
        X(j+nobox_z)=X(j+nobox_z-1)+delta_z(j-1)*Tfg;       % ��˱��¶��Ե����ݶ�����
        X(j+2*nobox_z)=X(j+2*nobox_z-1)+delta_z(j-1)*Tfg;   % �������꾮Һ�¶��Եز��ݶ�����
    end

    for i=4:nobox_r     % ��4�㼰֮������¶ȳ�ֵ
        X((i-1)*nobox_z+1)=Tf0; % �����¶�Ϊ�ر��¶�
        for j=2:nobox_z
            X((i-1)*nobox_z+j)=X((i-1)*nobox_z+j-1)+delta_z(j-1)*Tfg;   % �¶��Ե����ݶ�����
        end 
    end

    Tf=zeros(nobox_z,1);    % ԭʼ�ز��¶�
    Tf(1)=Tf0;              % �����¶�Ϊ�ر��¶�
    for j=2:nobox_z
        Tf(j)=Tf(j-1)+delta_z(j-1)*Tfg; % �¶��Ե����ݶ�����
    end

    X0=zeros(nobox_r*nobox_z,1);    % �洢ÿ�μ�����¶ȳ����Աȵ����¶ȳ��Ĳ�ֵ
    BC=zeros(nobox_r*nobox_z,1);    % �¶ȿ��Ʒ��̵��Ҷ�������
    X1=zeros(nobox_r*nobox_z,timesteps);	% �洢ÿһ��ʱ�䲽������¶ȳ�

    for nt=1:timesteps

        %% �¶ȿ��Ʒ��̵��Ҷ�����������
        BC(1)=-B99(1)*Tm0+F99(1)*X(1);
        for i=2:(nobox_r-1)*nobox_z
            BC(i)=F99(i)*X(i);
        end
        for i=(nobox_r-1)*nobox_z+1:nobox_r*nobox_z
            BC(i)=F99(i)*X(i)-E99(i)*Tf(i-(nobox_r-1)*nobox_z);
        end

        max1=1;         % ��������ֵ��Ϊ�˽���ѭ��
        err=0.001;      % ���������
        k=0;            % ��������������

        while max1>err
            k=k+1;
            max1=0;     % ����ѭ����max1��ֵΪ0����������õ���max1Ϊ���������еõ���

            X0=X;       % ��X��ֵ��X0�洢

            %% �����µ�X
            X(1)=(BC(1)-E99(1)*X(1+nobox_z))/A99(1);
            for j=2:nobox_z
                X(j)=(BC(j)-B99(j)*X(j-1)-E99(j)*X(j+nobox_z))/A99(j);
            end

            X(2*nobox_z)=X(nobox_z);
            for j=2*nobox_z-1:-1:nobox_z+2
                X(j)=(BC(j)-B99(j)*X(j-1)-C99(j)*X(j+1)-D99(j)*X(j-nobox_z)-E99(j)*X(j+nobox_z))/A99(j);
            end
            X(nobox_z+1)=(BC(nobox_z+1)-C99(nobox_z+1)*X(nobox_z+2)-D99(nobox_z+1)*X(1)-E99(nobox_z+1)*X(2*nobox_z+1))/A99(nobox_z+1);

            X(3*nobox_z)=X(nobox_z);
            for j=3*nobox_z-1:-1:2*nobox_z+1
                X(j)=(BC(j)-C99(j)*X(j+1)-D99(j)*X(j-nobox_z)-E99(j)*X(j+nobox_z))/A99(j);
            end

            for i=4:nobox_r-1
                X((i-1)*nobox_z+1)=(BC((i-1)*nobox_z+1)-C99((i-1)*nobox_z+1)*X((i-1)*nobox_z+2)-D99((i-1)*nobox_z+1)*X((i-2)*nobox_z+1)-E99((i-1)*nobox_z+1)*X(i*nobox_z+1))/A99((i-1)*nobox_z+1);
                for j=(i-1)*nobox_z+2:i*nobox_z-1
                    X(j)=(BC(j)-B99(j)*X(j-1)-C99(j)*X(j+1)-D99(j)*X(j-nobox_z)-E99(j)*X(j+nobox_z))/A99(j);
                end
                X(i*nobox_z)=(BC(i*nobox_z)-B99(i*nobox_z)*X(i*nobox_z-1)-D99(i*nobox_z)*X((i-1)*nobox_z)-E99(i*nobox_z)*X((i+1)*nobox_z))/A99(i*nobox_z);
            end

            X((nobox_r-1)*nobox_z+1)=(BC((nobox_r-1)*nobox_z+1)-C99((nobox_r-1)*nobox_z+1)*X((nobox_r-1)*nobox_z+2)-D99((nobox_r-1)*nobox_z+1)*X((nobox_r-2)*nobox_z+1))/A99((nobox_r-1)*nobox_z+1);
            for j=(nobox_r-1)*nobox_z+2:nobox_r*nobox_z-1
                X(j)=(BC(j)-B99(j)*X(j-1)-C99(j)*X(j+1)-D99(j)*X(j-nobox_z))/A99(j);
            end
            X(nobox_r*nobox_z)=(BC(nobox_r*nobox_z)-B99(nobox_r*nobox_z)*X(nobox_r*nobox_z-1)-D99(nobox_r*nobox_z)*X((nobox_r-1)*nobox_z))/A99(nobox_r*nobox_z);

            %% �Ա��¼����X�����ǰ��X֮��������õ�������ֵ
            for i=1:nobox_r*nobox_z
                m=abs(X0(i)-X(i));
                if m>max1
                    max1=m;
                end
            end
        end
        X1(1:nobox_r*nobox_z,nt)=X(1:nobox_r*nobox_z);          % ����ǰʱ�䲽���������õ����¶ȳ���ֵ��X1
        T_Drillpipe(1:nobox_z,nt)=X1(1:nobox_z,nt);             % ��X1����ȡ��������꾮Һ�¶�
        T_Annulus(1:nobox_z,nt)=X1(2*nobox_z+1:3*nobox_z,nt);   % ��X1����ȡ���������꾮Һ�¶�
    end
else
    
    for j=1:nobox_z
        if Depth(j)<Depth_weiguan
            Rho(j,4)=rou_m;
            Cp(j,4)=Cpm;
            lambda(j,4)=lambdam;
            Rho(j,5)=rou_m;
            Cp(j,5)=Cpm;
            lambda(j,5)=lambdam;
        else
            Rho(j,4)=rou_s;
            Cp(j,4)=Cps;
            lambda(j,4)=lambdas;
            Rho(j,5)=rou_c;
            Cp(j,5)=Cpc;
            lambda(j,5)=lambdac;
        end
    end
    
    if casing_number==0
        Rho(:,6:nobox_r)=rou_f;
        Cp(:,6:nobox_r)=Cpf;
        lambda(:,6:nobox_r)=lambdaf;
    else
        for i=1:casing_number
            for j=1:nobox_z
                if Depth(j)<Depth_casing(i)
                    Rho(j,5+2*i-1)=rou_s;
                    Cp(j,5+2*i-1)=Cps;
                    lambda(j,5+2*i-1)=lambdas;
                elseif Depth(j)>=Depth_casing(i) && Depth(j)<=welldepth
                    Rho(j,5+2*i-1)=rou_f;
                    Cp(j,5+2*i-1)=Cpf;
                    lambda(j,5+2*i-1)=lambdaf;
                end
                
                if Depth(j)<Depth_fluid(i)
                    Rho(j,5+2*i)=rou_m;
                    Cp(j,5+2*i)=Cpm;
                    lambda(j,5+2*i)=lambdam;
                elseif Depth(j)>=Depth_fluid(i) && Depth(j)<Depth_casing(i)
                    Rho(j,5+2*i)=rou_c;
                    Cp(j,5+2*i)=Cpc;
                    lambda(j,5+2*i)=lambdac;
                elseif Depth(j)>=Depth_casing(i) && Depth(j)<=welldepth
                    Rho(j,5+2*i)=rou_f;
                    Cp(j,5+2*i)=Cpf;
                    lambda(j,5+2*i)=lambdaf;
                end
            end
        end
        for i=5+2*casing_number+1:nobox_r
            for j=1:nobox_z
                Rho(j,i)=rou_f;
                Cp(j,i)=Cpf;
                lambda(j,i)=lambdaf;
            end
        end
    end
    
    if Rheology_type==1
        Re_drillpipe=Rho(:,1).*V_drillpipe.*(2*r(:,1))./miu_m;                  % �������ŵ��
        h_drillpipe=0.027.*lambda(:,1).*Re_drillpipe.^0.8.*(miu_drillpipe.*Cp(:,1)./lambda(:,1)).^0.4./(2*r(:,1))/5;        % �꾮Һ������ڱڶ�������ϵ����W/(m^2���棩   
        Re_annulus=zeros(nobox_z,1);
        h_annulus=zeros(nobox_z,1);
        for j=1:nobox_z
            if Depth(j)<Depth_weiguan
                Re_annulus(j)=Rho(j,3)*V_annulus(j)*(2*(r(j,5)-r(j,2)))/miu_m;              % ��������ŵ��
                h_annulus(j)=0.027*lambda(j,3)*Re_annulus(j)^0.8*(miu_drillpipe(j)*Cp(j,3)/lambda(j,3))^0.4/(2*(r(j,5)-r(j,2)))/5;	% �꾮Һ�������ڶ�������ϵ����W/(m^2���棩
            else
                Re_annulus(j)=Rho(j,3)*V_annulus(j)*(2*(r(j,3)-r(j,2)))/miu_m;              % ��������ŵ��
                h_annulus(j)=0.027*lambda(j,3)*Re_annulus(j)^0.8*(miu_drillpipe(j)*Cp(j,3)/lambda(j,3))^0.4/(2*(r(j,3)-r(j,2)))/5;	% �꾮Һ�������ڶ�������ϵ����W/(m^2���棩
            end
        end
    elseif Rheology_type==2
        Re_drillpipe=Rho(:,1).*V_drillpipe.^(2-n_power).*(2*r(:,1)).^n_power./(8^(n_power-1).*K_power.*((1+3*n_power)./(4*n_power)).^n_power);                  % �������ŵ��
        delta_correct=(3*n_power+1)/(4*n_power);
        h_drillpipe=0.091*Re_drillpipe.^0.87.*(miu_drillpipe.*Cp(:,1)./lambda(:,1)).^(1/3)*delta_correct./(2*r(:,1))/5;
        Re_annulus=zeros(nobox_z,1);
        h_annulus=zeros(nobox_z,1);
        for j=1:nobox_z
            if Depth(j)<Depth_weiguan
                Re_annulus(j)=Rho(j,3)*V_annulus(j)^(2-n_power)*(2*(r(j,5)-r(j,2)))^n_power/(8^(n_power-1)*K_power*((1+3*n_power)/(4*n_power))^n_power);              % ��������ŵ��
                h_annulus(j)=0.091*Re_annulus(j)^0.87*(miu_annulus(j)*Cp(j,3)./lambda(j,3))^(1/3)*delta_correct/(2*(r(j,5)-r(j,2)))/5;
            else
                Re_annulus(j)=Rho(j,3)*V_annulus(j)^(2-n_power)*(2*(r(j,3)-r(j,2)))^n_power/(8^(n_power-1)*K_power*((1+3*n_power)/(4*n_power))^n_power);              % ��������ŵ��
                h_annulus(j)=0.091*Re_annulus(j)^0.87*(miu_annulus(j)*Cp(j,3)./lambda(j,3))^(1/3)*delta_correct/(2*(r(j,3)-r(j,2)))/5;	% �꾮Һ�������ڶ�������ϵ����W/(m^2���棩
            end
        end
    end
    h1=h_drillpipe;                                                                                                     % �꾮Һ������ڱڶ�������ϵ����W/(m^2���棩    
    h2=h_annulus;                                                                                                       % �꾮Һ�������ڶ�������ϵ����W/(m^2���棩
    h3=h_annulus;                                                                                                       % �꾮Һ�뾮Ͳ�ڱڶ�������ϵ����W/(m^2���棩

    %% ϵ���������
    A=zeros(nobox_z,nobox_r);   % ϵ������A
    B=zeros(nobox_z,nobox_r);   % ϵ������B
    C=zeros(nobox_z,nobox_r);   % ϵ������C
    D=zeros(nobox_z,nobox_r);   % ϵ������D
    E=zeros(nobox_z,nobox_r);   % ϵ������E
    F=zeros(nobox_z,nobox_r);   % ϵ������F
    
    for j=1:nobox_z
        if j==1     % ��������ϵ���ļ���
            A(j,1)=-Rho(j,1)*q*Cp(j,1)/delta_z(j)-2*pi*r(j,1)*h1(j,1)-Rho(j,1)*Cp(j,1)*pi*r(j,1)^2/delta_t;
            B(j,1)=Rho(j,1)*q*Cp(j,1)/delta_z(j);
            C(j,1)=0;
            D(j,1)=0;
            E(j,1)=2*pi*r(j,1)*h1(j,1);
            F(j,1)=-Rho(j,1)*Cp(j,1)*pi*r(j,1)^2/delta_t;

            A(j,2)=-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j)+delta_z(j))*delta_z(j))-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j)+delta_z(j))*delta_z(j))-2*pi*r(j,1)*h1(j,1)-2*pi*r(j,2)*h2(j,1)-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;
            B(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j)+delta_z(j))*delta_z(j));
            C(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j)+delta_z(j))*delta_z(j));
            D(j,2)=2*pi*r(j,1)*h1(j,1);
            E(j,2)=2*pi*r(j,2)*h2(j,1);
            F(j,2)=-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;

            A(j,3)=-Rho(j,3)*q*Cp(j,3)/delta_z(j)-2*pi*r(j,2)*h2(j,1)-2*pi*r(j,3)*h3(j,1)-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;
            B(j,3)=0;
            C(j,3)=Rho(j,3)*q*Cp(j,3)/delta_z(j);
            D(j,3)=2*pi*r(j,2)*h2(j,1);
            E(j,3)=2*pi*r(j,3)*h3(j,1);
            F(j,3)=-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;

            A(j,4)=-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j)+delta_z(j))*delta_z(j))-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j)+delta_z(j))*delta_z(j))-2*pi*r(j,3)*h3(j,1)-pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)))-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;
            B(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j)+delta_z(j))*delta_z(j));
            C(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j)+delta_z(j))*delta_z(j));
            D(j,4)=2*pi*r(j,3)*h3(j,1);
            E(j,4)=pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)));
            F(j,4)=-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;
        else
            A(j,1)=-Rho(j,1)*q*Cp(j,1)/delta_z(j-1)-2*pi*r(j,1)*h1(j,1)-Rho(j,1)*Cp(j,1)*pi*r(j,1)^2/delta_t;
            B(j,1)=Rho(j,1)*q*Cp(j,1)/delta_z(j-1);
            C(j,1)=0;
            D(j,1)=0;
            E(j,1)=2*pi*r(j,1)*h1(j,1);
            F(j,1)=-Rho(j,1)*Cp(j,1)*pi*r(j,1)^2/delta_t;

            if j==nobox_z   % ��������ϵ���ļ���
                A(j,2)=-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-2*pi*r(j,1)*h1(j,1)-2*pi*r(j,2)*h2(j,1)-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;
                B(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));
                C(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));
                D(j,2)=2*pi*r(j,1)*h1(j,1);
                E(j,2)=2*pi*r(j,2)*h2(j,1);
                F(j,2)=-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;

                A(j,3)=-Rho(j,3)*q*Cp(j,3)/delta_z(j-1)-2*pi*r(j,2)*h2(j,1)-2*pi*r(j,3)*h3(j,1)-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;
                B(j,3)=0;
                C(j,3)=Rho(j,3)*q*Cp(j,3)/delta_z(j-1);
                D(j,3)=2*pi*r(j,2)*h2(j,1);
                E(j,3)=2*pi*r(j,3)*h3(j,1);
                F(j,3)=-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;

                A(j,4)=-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-2*pi*r(j,3)*h3(j,1)-pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)))-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;
                B(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));
                C(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));
                D(j,4)=2*pi*r(j,3)*h3(j,1);
                E(j,4)=pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)));
                F(j,4)=-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;
            else    % �м�����ϵ���ļ���
                A(j,2)=-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j))-2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j-1))-2*pi*r(j,1)*h1(j,1)-2*pi*r(j,2)*h2(j,1)-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;
                B(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j-1));
                C(j,2)=2*lambda(j,2)*pi*(r(j,2)^2-r(j,1)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j));
                D(j,2)=2*pi*r(j,1)*h1(j,1);
                E(j,2)=2*pi*r(j,2)*h2(j,1);
                F(j,2)=-Rho(j,2)*Cp(j,2)*pi*(r(j,2)^2-r(j,1)^2)/delta_t;

                A(j,3)=-Rho(j,3)*q*Cp(j,3)/delta_z(j)-2*pi*r(j,2)*h2(j,1)-2*pi*r(j,3)*h3(j,1)-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;
                B(j,3)=0;
                C(j,3)=Rho(j,3)*q*Cp(j,3)/delta_z(j);
                D(j,3)=2*pi*r(j,2)*h2(j,1);
                E(j,3)=2*pi*r(j,3)*h3(j,1);
                F(j,3)=-Rho(j,3)*Cp(j,3)*pi*(r(j,3)^2-r(j,2)^2)/delta_t;

                A(j,4)=-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j))-2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j-1))-2*pi*r(j,3)*h3(j,1)-pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)))-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;
                B(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j-1));
                C(j,4)=2*lambda(j,4)*pi*(r(j,4)^2-r(j,3)^2)/((delta_z(j-1)+delta_z(j))*delta_z(j));
                D(j,4)=2*pi*r(j,3)*h3(j,1);
                E(j,4)=pi*(r(j,4)^2-r(j,3)^2)*log(r(j,5)/r(j,3))/(r(j,4)*(r(j,4)-r(j,3))*log(r(j,5)/r(j,4))*(log(r(j,4)/r(j,3))/lambda(j,4)+log(r(j,5)/r(j,4))/lambda(j,5)));
                F(j,4)=-Rho(j,4)*Cp(j,4)*pi*(r(j,4)^2-r(j,3)^2)/delta_t;
            end
        end
    end

    for i=5:nobox_r     % �����5�㼰֮������ϵ������
        for j=1:nobox_z
            if j==1     % ��������ϵ���ļ���
                A(j,i)=-log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)))-log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)))-2*lambda(j,i)/((delta_z(j)+delta_z(j))*delta_z(j))-2*lambda(j,i)/((delta_z(j)+delta_z(j))*delta_z(j))-Rho(j,i)*Cp(j,i)/delta_t;
                B(j,i)=2*lambda(j,i)/((delta_z(j)+delta_z(j))*delta_z(j));
                C(j,i)=2*lambda(j,i)/((delta_z(j)+delta_z(j))*delta_z(j));
                D(j,i)=log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)));
                E(j,i)=log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)));
                F(j,i)=-Rho(j,i)*Cp(j,i)/delta_t;
            elseif j==nobox_z       % ��������ϵ���ļ���
                A(j,i)=-log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)))-log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)))-2*lambda(j,i)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-2*lambda(j,i)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1))-Rho(j,i)*Cp(j,i)/delta_t;
                B(j,i)=2*lambda(j,i)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));
                C(j,i)=2*lambda(j,i)/((delta_z(j-1)+delta_z(j-1))*delta_z(j-1));
                D(j,i)=log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)));
                E(j,i)=log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)));
                F(j,i)=-Rho(j,i)*Cp(j,i)/delta_t;
            else        % �м�����ϵ���ļ���
                A(j,i)=-log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)))-log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)))-2*lambda(j,i)/((delta_z(j-1)+delta_z(j))*delta_z(j))-2*lambda(j,i)/((delta_z(j-1)+delta_z(j))*delta_z(j-1))-Rho(j,i)*Cp(j,i)/delta_t;
                B(j,i)=2*lambda(j,i)/((delta_z(j-1)+delta_z(j))*delta_z(j-1));
                C(j,i)=2*lambda(j,i)/((delta_z(j-1)+delta_z(j))*delta_z(j));
                D(j,i)=log(r(j,i)/r(j,i-2))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i)/r(j,i-1))*(log(r(j,i-1)/r(j,i-2))/lambda(j,i-1)+log(r(j,i)/r(j,i-1))/lambda(j,i)));
                E(j,i)=log(r(j,i+1)/r(j,i-1))/(r(j,i)*(r(j,i)-r(j,i-1))*log(r(j,i+1)/r(j,i))*(log(r(j,i)/r(j,i-1))/lambda(j,i)+log(r(j,i+1)/r(j,i))/lambda(j,i+1)));
                F(j,i)=-Rho(j,i)*Cp(j,i)/delta_t;
            end
        end
    end

    A99=zeros(nobox_r*nobox_z,1);   % ϵ������A�洢������A99��
    B99=zeros(nobox_r*nobox_z,1);   % ϵ������B�洢������B99��
    C99=zeros(nobox_r*nobox_z,1);   % ϵ������C�洢������C99��
    D99=zeros(nobox_r*nobox_z,1);   % ϵ������D�洢������D99��
    E99=zeros(nobox_r*nobox_z,1);   % ϵ������E�洢������E99��
    F99=zeros(nobox_r*nobox_z,1);   % ϵ������F�洢������F99��

    % ϵ������A�洢������A99��
    A99(1:nobox_z)=A(:,1);
    A99(nobox_z+1)=A(1,2)+B(1,2);
    A99(nobox_z+2:2*nobox_z-1)=A(2:nobox_z-1,2);
    A99(2*nobox_z)=A(nobox_z,2)+C(nobox_z,2);
    A99(2*nobox_z+1:3*nobox_z-1)=A(1:nobox_z-1,3);
    A99(3*nobox_z)=A(nobox_z,3)+C(nobox_z,3);
    for i=4:nobox_r
        A99((i-1)*nobox_z+1)=A(1,i)+B(1,i);
        A99((i-1)*nobox_z+2:i*nobox_z-1)=A(2:nobox_z-1,i);
        A99(i*nobox_z)=A(nobox_z,i)+C(nobox_z,i);
    end

    % ϵ������B��C��D��E��F�洢������B99��C99��D99��E99��F99��
    for i=1:nobox_r
        B99((i-1)*nobox_z+1:i*nobox_z)=B(:,i);
        C99((i-1)*nobox_z+1:i*nobox_z)=C(:,i);
        D99((i-1)*nobox_z+1:i*nobox_z)=D(:,i);
        E99((i-1)*nobox_z+1:i*nobox_z)=E(:,i);
        F99((i-1)*nobox_z+1:i*nobox_z)=F(:,i);
    end

    %% ���¶ȳ���������ֵ
    X=zeros(nobox_r*nobox_z,1);     % �¶ȳ�����
    % ������꾮Һ����˱ڡ��������꾮Һ�¶�
    X(1)=Tm0;           % �����¶�Ϊ�꾮Һע���¶�
    X(1+nobox_z)=Tm0;   % �����¶�Ϊ�꾮Һע���¶�
    X(1+2*nobox_z)=Tm0; % �����¶�Ϊ�꾮Һע���¶�
    for j=2:nobox_z
        X(j)=X(j-1)+delta_z(j-1)*Tfg;                       % ������꾮Һ�¶��Ե����ݶ�����
        X(j+nobox_z)=X(j+nobox_z-1)+delta_z(j-1)*Tfg;       % ��˱��¶��Ե����ݶ�����
        X(j+2*nobox_z)=X(j+2*nobox_z-1)+delta_z(j-1)*Tfg;   % �������꾮Һ�¶��Եز��ݶ�����
    end

    for i=4:nobox_r     % ��4�㼰֮������¶ȳ�ֵ
        X((i-1)*nobox_z+1)=Tf0; % �����¶�Ϊ�ر��¶�
        for j=2:nobox_z
            X((i-1)*nobox_z+j)=X((i-1)*nobox_z+j-1)+delta_z(j-1)*Tfg;   % �¶��Ե����ݶ�����
        end 
    end

    Tf=zeros(nobox_z,1);    % ԭʼ�ز��¶�
    Tf(1)=Tf0;              % �����¶�Ϊ�ر��¶�
    for j=2:nobox_z
        Tf(j)=Tf(j-1)+delta_z(j-1)*Tfg; % �¶��Ե����ݶ�����
    end

    X0=zeros(nobox_r*nobox_z,1);    % �洢ÿ�μ�����¶ȳ����Աȵ����¶ȳ��Ĳ�ֵ
    BC=zeros(nobox_r*nobox_z,1);    % �¶ȿ��Ʒ��̵��Ҷ�������
    X1=zeros(nobox_r*nobox_z,timesteps);	% �洢ÿһ��ʱ�䲽������¶ȳ�
    
    for nt=1:timesteps

        %% �¶ȿ��Ʒ��̵��Ҷ�����������
        BC(1)=-B99(1)*Tm0+F99(1)*X(1);
        for i=2:(nobox_r-1)*nobox_z
            BC(i)=F99(i)*X(i);
        end
        for i=(nobox_r-1)*nobox_z+1:nobox_r*nobox_z
            BC(i)=F99(i)*X(i)-E99(i)*Tf(i-(nobox_r-1)*nobox_z);
        end

        max1=1;         % ��������ֵ��Ϊ�˽���ѭ��
        err=0.001;      % ���������
        k=0;            % ��������������

        while max1>err
            k=k+1;
            max1=0;     % ����ѭ����max1��ֵΪ0����������õ���max1Ϊ���������еõ���

            X0=X;       % ��X��ֵ��X0�洢

            %% �����µ�X
            X(1)=(BC(1)-E99(1)*X(1+nobox_z))/A99(1);
            for j=2:nobox_z
                X(j)=(BC(j)-B99(j)*X(j-1)-E99(j)*X(j+nobox_z))/A99(j);
            end

            X(2*nobox_z)=X(nobox_z);
            for j=2*nobox_z-1:-1:nobox_z+2
                X(j)=(BC(j)-B99(j)*X(j-1)-C99(j)*X(j+1)-D99(j)*X(j-nobox_z)-E99(j)*X(j+nobox_z))/A99(j);
            end
            X(nobox_z+1)=(BC(nobox_z+1)-C99(nobox_z+1)*X(nobox_z+2)-D99(nobox_z+1)*X(1)-E99(nobox_z+1)*X(2*nobox_z+1))/A99(nobox_z+1);

            X(3*nobox_z)=X(nobox_z);
            for j=3*nobox_z-1:-1:2*nobox_z+1
                X(j)=(BC(j)-C99(j)*X(j+1)-D99(j)*X(j-nobox_z)-E99(j)*X(j+nobox_z))/A99(j);
            end

            X(3*nobox_z+1:3*nobox_z+Weiguan_i)=X(2*nobox_z+1:2*nobox_z+Weiguan_i);
            for j=Weiguan_i+1:nobox_z-1
                X(3*nobox_z+j)=(BC(3*nobox_z+j)-B99(3*nobox_z+j)*X(3*nobox_z+j-1)-C99(3*nobox_z+j)*X(3*nobox_z+j+1)-D99(3*nobox_z+j)*X(2*nobox_z+j)-E99(3*nobox_z+j)*X(4*nobox_z+j))/A99(3*nobox_z+j);
            end
            X(4*nobox_z)=(BC(4*nobox_z)-B99(4*nobox_z)*X(4*nobox_z-1)-D99(4*nobox_z)*X(3*nobox_z)-E99(4*nobox_z)*X(5*nobox_z))/A99(4*nobox_z);
            
            X(4*nobox_z+1:4*nobox_z+Weiguan_i)=X(2*nobox_z+1:2*nobox_z+Weiguan_i);
            X(4*nobox_z+1:4*nobox_z+Weiguan_i)=X(2*nobox_z+1:2*nobox_z+Weiguan_i);
            for j=Weiguan_i+1:nobox_z-1
                X(4*nobox_z+j)=(BC(4*nobox_z+j)-B99(4*nobox_z+j)*X(4*nobox_z+j-1)-C99(4*nobox_z+j)*X(4*nobox_z+j+1)-D99(4*nobox_z+j)*X(3*nobox_z+j)-E99(4*nobox_z+j)*X(5*nobox_z+j))/A99(4*nobox_z+j);
            end
            X(5*nobox_z)=(BC(5*nobox_z)-B99(5*nobox_z)*X(5*nobox_z-1)-D99(5*nobox_z)*X(4*nobox_z)-E99(5*nobox_z)*X(6*nobox_z))/A99(5*nobox_z);
         
            for i=6:nobox_r-1
                X((i-1)*nobox_z+1)=(BC((i-1)*nobox_z+1)-C99((i-1)*nobox_z+1)*X((i-1)*nobox_z+2)-D99((i-1)*nobox_z+1)*X((i-2)*nobox_z+1)-E99((i-1)*nobox_z+1)*X(i*nobox_z+1))/A99((i-1)*nobox_z+1);
                for j=(i-1)*nobox_z+2:i*nobox_z-1
                    X(j)=(BC(j)-B99(j)*X(j-1)-C99(j)*X(j+1)-D99(j)*X(j-nobox_z)-E99(j)*X(j+nobox_z))/A99(j);
                end
                X(i*nobox_z)=(BC(i*nobox_z)-B99(i*nobox_z)*X(i*nobox_z-1)-D99(i*nobox_z)*X((i-1)*nobox_z)-E99(i*nobox_z)*X((i+1)*nobox_z))/A99(i*nobox_z);
            end

            X((nobox_r-1)*nobox_z+1)=(BC((nobox_r-1)*nobox_z+1)-C99((nobox_r-1)*nobox_z+1)*X((nobox_r-1)*nobox_z+2)-D99((nobox_r-1)*nobox_z+1)*X((nobox_r-2)*nobox_z+1))/A99((nobox_r-1)*nobox_z+1);
            for j=(nobox_r-1)*nobox_z+2:nobox_r*nobox_z-1
                X(j)=(BC(j)-B99(j)*X(j-1)-C99(j)*X(j+1)-D99(j)*X(j-nobox_z))/A99(j);
            end
            X(nobox_r*nobox_z)=(BC(nobox_r*nobox_z)-B99(nobox_r*nobox_z)*X(nobox_r*nobox_z-1)-D99(nobox_r*nobox_z)*X((nobox_r-1)*nobox_z))/A99(nobox_r*nobox_z);

            %% �Ա��¼����X�����ǰ��X֮��������õ�������ֵ
            for i=1:nobox_r*nobox_z
                m=abs(X0(i)-X(i));
                if m>max1
                    max1=m;
                end
            end
        end
        X1(1:nobox_r*nobox_z,nt)=X(1:nobox_r*nobox_z);          % ����ǰʱ�䲽���������õ����¶ȳ���ֵ��X1
        T_Drillpipe(1:nobox_z,nt)=X1(1:nobox_z,nt);             % ��X1����ȡ��������꾮Һ�¶�
        T_Annulus(1:nobox_z,nt)=X1(2*nobox_z+1:3*nobox_z,nt);   % ��X1����ȡ���������꾮Һ�¶�
    end
end

for t=1:1:timesteps
    Time(t)=t*delta_t; % ��ʱ��ڵ��Ӧʱ�䣬s
end

end