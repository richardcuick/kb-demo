function [PressureSurge_d,flow_pattern_d]=PressureSurgeDrillPipe(rho_d,V_d,eta_infinity_d,tau_y_c_d,D_d_i)
    % ������˽ṹ���꾮Һ���ʺ������������������Ħ��ѹ��������
    % ���������rho_dΪ������꾮Һ�ܶȣ�kg/m^3����V_dΪ������꾮Һ���٣�m/s����A_dΪ����ڽ������m^2����eta_infinity_dΪ������꾮Һ�Ŀ�ɭ�ȣ�Pa��s����tau_y_c_dΪ������꾮Һ�Ŀ�ɭ����Ӧ����Pa����D_d_iΪ����ھ���m��
    % ���������Ff_dΪ����е�λ�����ϵĲ���ѹ����Pa/m����flow_pattern_dΪ������꾮Һ����
    
    %% ʹ�õ�����������洦������Ӧ��
    Qv=V_d*(0.25*pi*D_d_i^2); % �꾮Һ���������m^3/s
    tau_w_old=100; % ������ֵ
    loop=1;
    err_tau=1;
    while abs(err_tau)>1e-4&&loop<=1000 % ������ֹ����
        int=(4*Qv*eta_infinity_d)/(pi*(D_d_i/2)^3)-(4*tau_y_c_d)/3+(16*(tau_y_c_d*tau_w_old)^0.5)/7+(tau_y_c_d^4)/(21*tau_w_old^3);
        tau_w_new=int;
        err_tau=abs(tau_w_old-tau_w_new)/tau_w_old;
        tau_w_old=tau_w_new;
        loop=loop+1;
    end
    tau_w=tau_w_new; % �������Ӧ����Pa
    
    gammadot_w=(tau_w-2*((tau_w*tau_y_c_d)^0.5)+tau_y_c_d)/eta_infinity_d; % ����������ʣ�s^-1
    
    nprime=(8*V_d/D_d_i)/(4*gammadot_w-24*V_d/D_d_i); % ��������ָ��
    
    Re_gcl_d=3470-1370*nprime; % �����ٽ������ŵ��
    Re_gct_d=4270-1370*nprime; % �����ٽ������ŵ��
    
    D_eff=(4*nprime*D_d_i)/(3*nprime+1); % ��Ч�ܾ�
    mu_aw=tau_w/gammadot_w; % ����������ʣ�s^-1
    
    Re_g_d=(rho_d*D_eff*V_d)/mu_aw; % ������ŵ��
    
    %% ��������в�ͬ�����ͣ�����������꾮Һ�ڵ�λ�����ϵĲ���ѹ��
    if Re_g_d==0 % ���徲ֹ
        f_d=0; % �������徲ֹ�����£�������꾮ҺĦ������
        PressureSurge_d=0; % �������徲ֹ�����£�������꾮Һ�ڵ�λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_d=0; % flow_pattern_dΪ������ͣ�0��ʾ���徲ֹ
    elseif Re_g_d>0 && Re_g_d<=Re_gcl_d % ����
        f_d=16/Re_g_d; % ������������£�������꾮ҺĦ������
        PressureSurge_d=2*f_d*rho_d*V_d^2/D_d_i; % ������������£�������꾮Һ�ڵ�λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_d=1; % flow_pattern_dΪ������ͣ�1��ʾ����
    elseif Re_g_d>Re_gct_d % ����
        f_d_old=100; % ������ֵ
        err_f=1;
        loop=1;
        while abs(err_f)>1e-4&&loop<=1000 % ������ֹ����
            int=(4/nprime^0.75)*log10(Re_g_d*f_d_old^(1-nprime/2))-0.395/nprime^1.2;
            f_d_new=(int^-1)^2;
            err_f=abs(f_d_old-f_d_new)/f_d_old;
            f_d_old=f_d_new;
            loop=loop+1;
        end
        f_d=f_d_new; % �������Ӧ����Pa
        PressureSurge_d=2*f_d*rho_d*V_d^2/D_d_i; % ������������£�������꾮Һ�ڵ�λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_d=3;
    else
        Re_g_d=Re_gcl_d;
        f_d1=16/Re_g_d;
        
        Re_g_d=Re_gct_d;
        f_d_old=100; % ������ֵ
        loop=1;
        err_f=1;
        while abs(err_f)>1e-4&&loop<=1000 % ������ֹ����
            int=(4/nprime^0.75)*log10(Re_g_d*f_d_old^(1-nprime/2))-0.395/nprime^1.2;
            f_d_new=(int^-1)^2;
            err_f=abs(f_d_old-f_d_new)/f_d_old;
            f_d_old=f_d_new;
            loop=loop+1;
        end
        f_d2=f_d_new; % �������Ӧ����Pa
        
        f_d=f_d1+((Re_g_d-Re_gcl_d)*(f_d2-f_d1))/(Re_gct_d-Re_gcl_d);
        PressureSurge_d=2*f_d*rho_d*V_d^2/D_d_i; % ������������£�������꾮Һ�ڵ�λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_d=2;
    end
end