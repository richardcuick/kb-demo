function [Ff_d,flow_pattern_d]=Friction_drillpipe(rho_d,V_d,mu_d,D_di)
    % ������˽ṹ���꾮Һ���ʺ������������������Ħ��ѹ��������
    % ���������rho_dΪ������꾮Һ�ܶȣ�kg/m^3����V_dΪ������꾮Һ���٣�m/s����mu_dΪ������꾮Һ�������ȣ�Pa��s����D_d_iΪ����ھ���m��
    % ���������Ff_dΪ����е�λ�����ϵ�Ħ��ѹ����Pa/m����flow_pattern_dΪ������꾮Һ����
    
    Re_d=rho_d*V_d*D_di/mu_d; % ���������ŵ��
    if Re_d<=2100 % ����
        f_d=16/Re_d; % ������������·���Ħ������
        flow_pattern_d=1; % flow_pattern_dΪ������ͣ�1��ʾ����
    else % ����
        f_d=0.0791/Re_d^0.25; % �������������·���Ħ������
        flow_pattern_d=3; % flow_pattern_dΪ������ͣ�3��ʾ����
    end
    Ff_d=2*f_d*V_d^2*rho_d/D_di; % ���㻷�����꾮Һ�ڵ�λ�����ϵ�Ħ��ѹ��
end