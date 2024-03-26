function [PressureSurge_d,flow_pattern_d]=PressureSurgeDrillPipe(rho_d,V_d,mu_p_d,tau_y,D_d_i)
    % ���㲨��ѹ��
    % ���������rho_dΪ������꾮Һ�ܶȣ�kg/m^3����V_dΪ������꾮Һ���٣�m/s����mu_p_dΪ������꾮Һ�������ȣ�Pa��s����tau_yΪ������꾮Һ������Ӧ����Pa����D_d_iΪ����ھ���m��
    % ���������PressureSurge_dΪ����е�λ�����ϵĲ���ѹ����Pa/m����flow_pattern_dΪ������꾮Һ����
    
    %% ʹ�õ�����������洦������Ӧ��
    Qv=V_d*(0.25*pi*D_d_i^2); % �꾮Һ���������m^3/s
    tau_w_old=100; % �ٶ����洦����Ӧ���ĳ�ʼֵ
    err_tau=1; % �������洦����Ӧ���ĳ�ʼ��ʹ��������ѭ��
    while abs(err_tau)>1e-4 % ���洦����Ӧ���������жϱ�׼
        tau_w_new=4*mu_p_d*Qv/(pi*(D_d_i/2)^3)+4/3*tau_y-1/3*tau_y^4/tau_w_old^3; % �����µı�������Ӧ��ֵ
        err_tau=abs(tau_w_new-tau_w_old)/tau_w_old; % ��������Ӧ�����ֵ
        tau_w_old=tau_w_new; % ���¼���ı��洦����Ӧ��ֵ��ֵ������ı��洦����Ӧ��ֵ
    end
    
    tau_w=tau_w_new; % �õ���������Ӧ������ֵ
    
    int1=tau_y/tau_w; % int1Ϊ����Ӧ���ͱ�������Ӧ���ı�ֵ�������������
    He=16800*int1/(1-int1)^3; % ����Hedstrom��
    
   %% ��������в�ͬ�����ͣ�����������꾮Һ�ڵ�λ�����ϵĲ���ѹ��
    Re_d_c=(1-4/3*int1+1/3*int1^4)/(8*int1)*He; % Re_d_cΪ�ٽ������ŵ��
    Re_d=rho_d*V_d*D_d_i/mu_p_d; % Re_dΪ�����ŵ��
    
    if Re_d==0 % ���徲ֹ
        PressureSurge_d=0; % �������徲ֹ�����£�������꾮Һ��λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_d=0; % ������ͣ�0��ʾ���徲ֹ
    elseif Re_d>0 && Re_d<=Re_d_c % ����
        PressureSurge_d=6895*mu_p_d*V_d/(216*D_d_i^2)+5.33355*tau_y/D_d_i; % ������������£�������꾮Һ�ڵ�λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_d=1; % ������ͣ�1��ʾ����
    else % ����
        PressureSurge_d=0.158278*rho_d^0.75*V_d^1.75*mu_p_d^0.25/D_d_i^1.25; % �������������£�������꾮Һ�ڵ�λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_d=3; % ������ͣ�3��ʾ����
    end
end