function [PressureSurge_d,flow_pattern_d]=PressureSurgeDrillpipe(rho_d,V_d,mu_d,D_d_i)
    % ���㲨��ѹ��
    % ���������rho_dΪ������꾮Һ�ܶȣ�kg/m^3����V_dΪ������꾮Һ���٣�m/s����mu_dΪ������꾮Һ�������ȣ�Pa��s����D_d_iΪ����ھ���m��
    % ���������PressureSurge_dΪ����е�λ�����ϵĲ���ѹ����Pa/m����flow_pattern_dΪ������꾮Һ����
    
    Re_d=rho_d*V_d*D_d_i/mu_d; % ���������ŵ��
    if Re_d==0 % ���徲ֹ
        f_d=0; % �������徲ֹ�����·���Ħ������
        flow_pattern_d=0; % �������ͣ�0��ʾ���徲ֹ
    elseif Re_d>0 && Re_d<=2100 % ����
        f_d=16/Re_d; % ������������·���Ħ������
        flow_pattern_d=1; % �������ͣ�1��ʾ����
    else % ����
        f_d=0.0791/Re_d^0.25; % �������������·���Ħ������
        flow_pattern_d=3; % �������ͣ�3��ʾ����
    end
    PressureSurge_d=2*f_d*V_d^2*rho_d/D_d_i; % ����������꾮Һ��λ�����ϵĲ���ѹ����Pa/m
end