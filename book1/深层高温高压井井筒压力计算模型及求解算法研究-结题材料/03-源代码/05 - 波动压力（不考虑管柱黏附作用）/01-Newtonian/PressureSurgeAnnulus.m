function [PressureSurge_a,flow_pattern_a]=PressureSurgeAnnulus(rho_a,V_a,mu_a,D_w,D_d_o)
    % ���㲨��ѹ��
    % ���������rho_aΪ�������꾮Һ�ܶȣ�kg/m^3����V_aΪ�������꾮Һ���٣�m/s����mu_aΪ�������꾮Һ���ȣ�Pa��s����D_wΪ��Ͳ�ھ���m����D_d_oΪ�����⾶��m��
    % ���������PressureSurge_aΪ�����е�λ�����ϵĲ���ѹ����Pa/m����flow_pattern_aΪ��������
    
    Re_a=0.816*rho_a*V_a*(D_w-D_d_o)/mu_a; % ���㻷����ŵ��
    if Re_a==0 % ���徲ֹ
        f_a=0; % �������徲ֹ�����·���Ħ������
        flow_pattern_a=0; % �������ͣ�0��ʾ���徲ֹ
    elseif Re_a>0 && Re_a<=2100 % ����
        f_a=16/Re_a; % ������������·���Ħ������
        flow_pattern_a=1; % �������ͣ�1��ʾ����
    else % ����
        f_a=0.0791/Re_a^0.25; % �������������·���Ħ������
        flow_pattern_a=3; % �������ͣ�3��ʾ����
    end
    PressureSurge_a=2*f_a*V_a^2*rho_a/(0.816*(D_w-D_d_o)); % ���㻷�����꾮Һ��λ�����ϵĲ���ѹ����Pa/m
end