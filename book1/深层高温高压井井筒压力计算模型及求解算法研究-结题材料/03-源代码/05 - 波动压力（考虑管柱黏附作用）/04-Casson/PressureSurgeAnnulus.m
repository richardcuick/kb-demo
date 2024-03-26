function [PressureSurge_a,flow_pattern_a]=PressureSurgeAnnulus(rho_a,V_a,eta_infinity_a,tau_y_c_a,D_w,D_d_o)
    % ���ݻ��սṹ���꾮Һ���ʺ��������������㲨��ѹ��
    % ���������rho_aΪ�����꾮Һ�ܶȣ�kg/m^3����V_aΪ�����꾮Һ���٣�m/s����eta_infinity_aΪ�����꾮Һ��ɭ�ȣ�Pa��s����tau_y_c_aΪ�����꾮Һ��ɭ����Ӧ����Pa����D_wΪ��Ͳ�ھ���m����D_d_oΪ�����⾶��m����QmΪ�꾮Һ����������kg/s��
    % ���������PressureSurge_aΪ�����е�λ�����ϵĲ���ѹ����Pa/m����flow_pattern_aΪ��������
    
    Qv=V_a*(0.25*pi*(D_w^2-D_d_o^2)); % �꾮Һ���������m^3/s
    D_h=D_w-D_d_o; % ˮ��ֱ����m
    
    %% ʹ�õ�����������洦�ļ���Ӧ��
    tau_w_old=100; % �ٶ����洦����Ӧ���ĳ�ʼֵ
    loop=1;
    err_tau=1;
    while abs(err_tau)>1e-4 && loop<=1000 % ������ֹ����
        int=(6*Qv*eta_infinity_a)/(pi*(D_w/2+D_d_o/2)*(D_h/2)^2)-(3*tau_y_c_a)/2+(5*tau_y_c_a^3)/(2*tau_w_old^2)+(12*tau_y_c_a^0.5*(tau_w_old^2.5-tau_y_c_a^2.5))/(5*tau_w_old^2);
        tau_w_new=int;
        err_tau=abs(tau_w_old-tau_w_new)/tau_w_old; % �������Ӧ�����ֵ
        tau_w_old=tau_w_new;
        loop=loop+1;
    end
    tau_w=tau_w_new; % �������Ӧ����Pa
    
    gammadot_w=(tau_w-2*((tau_w*tau_y_c_a)^0.5)+tau_y_c_a)/eta_infinity_a; % ����������ʣ�s^-1
    
    nprime=(4*V_a/D_h)/(gammadot_w-8*V_a/D_h); % ��������ָ��
    
    Re_gcl_a=3470-1370*nprime; % �����ٽ������ŵ��
    Re_gct_a=4270-1370*nprime; % �����ٽ������ŵ��
    
    D_eff=(2*nprime*D_h)/(2*nprime+1); % ��Ч�ܾ�
    mu_aw=tau_w/gammadot_w; % �������Ӧ����Pa��
    
    Re_g_a=(rho_a*D_eff*V_a)/mu_aw; % ������ŵ��
    
    %% ���ݻ����в�ͬ�����ͣ����㻷�����꾮Һ�ڵ�λ�����ϵĲ���ѹ��   
    if Re_g_a<=Re_gcl_a % ����
        f_a=16/Re_g_a; % ������������£��������꾮ҺĦ������
        PressureSurge_a=2*f_a*rho_a*V_a^2/D_h; % ������������£��������꾮Һ�ڵ�λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_a=1; % flow_pattern_aΪ�������ͣ�1��ʾ����
    elseif Re_g_a>Re_gct_a % ����
        f_a_old=100; % ������ֵ
        loop=1;
        err_tau=1;
        while abs(err_tau)>1e-4 && loop<=1000 % ������ֹ����
            int=(4/nprime^0.75)*log10(Re_g_a*f_a_old^(1-nprime/2))-0.395/nprime^1.2;
            f_a_c_new=(int^-1)^2;
            err_tau=abs(f_a_old-f_a_c_new)/f_a_old;
            f_a_old=f_a_c_new;
            loop=loop+1;
        end
        f_a=f_a_c_new; % �������������£��������꾮ҺĦ������
        PressureSurge_a=2*f_a*rho_a*V_a^2/D_h; % �������������£��������꾮Һ�ڵ�λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_a=3; % flow_pattern_aΪ�������ͣ�3��ʾ����
    else
        Re_g_a=Re_gcl_a;
        f_a1=16/Re_g_a;
        
        Re_g_a=Re_gct_a;
        f_a_old=100; % ������ֵ
        loop=1;
        err_tau=1;
        while abs(err_tau)>1e-4&loop<=1000 % ������ֹ����
            int=(4/nprime^0.75)*log10(Re_g_a*f_a_old^(1-nprime/2))-0.395/nprime^1.2;
            f_a_c_new=(int^-1)^2;
            err_tau=abs(f_a_old-f_a_c_new)/f_a_old;
            f_a_old=f_a_c_new;
            loop=loop+1;
        end
        f_a2=f_a_c_new;

        f_a=f_a1+((Re_g_a-Re_gcl_a)*(f_a2-f_a1))/(Re_gct_a-Re_gcl_a); % ��������������£��������꾮ҺĦ������
        PressureSurge_a=2*f_a*rho_a*V_a^2/D_h; % ��������������£��������꾮Һ�ڵ�λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_a=2; % flow_pattern_aΪ�������ͣ�2��ʾ������
    end
end