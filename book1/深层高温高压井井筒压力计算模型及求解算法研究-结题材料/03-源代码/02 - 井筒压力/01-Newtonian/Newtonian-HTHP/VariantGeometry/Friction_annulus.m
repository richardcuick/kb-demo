function [Ff_a,flow_pattern_a]=Friction_annulus(rho_a,V_a,mu_a,D_w,D_d_o)
    % ���ݻ��սṹ���꾮Һ���ʺ��������������㻷��Ħ��ѹ��������
    % ���������rho_aΪ�������꾮Һ�ܶȣ�kg/m^3����V_aΪ�������꾮Һ���٣�m/s����mu_aΪ�������꾮Һ���ȣ�Pa��s����D_wΪ��Ͳ�ھ���m����D_d_oΪ����⾶��m��
    % ���������Ff_aΪ�����е�λ�����ϵ�Ħ��ѹ����Pa/m����flow_pattern_aΪ�������꾮Һ����
    
    Re_a=0.816*rho_a*V_a*(D_w-D_d_o)/mu_a; % ���㻷����ŵ��
    if Re_a<=2100 % ����
        f_a=16/Re_a; % ������������·���Ħ������
        flow_pattern_a=1; % flow_pattern_aΪ�������ͣ�1��ʾ����
    else % ����
        f_a=0.0791/Re_a^0.25; % �������������·���Ħ������
        flow_pattern_a=3; % flow_pattern_aΪ�������ͣ�3��ʾ����
    end
    Ff_a=2*f_a*V_a^2*rho_a/(0.816*(D_w-D_d_o)); % ���㻷�����꾮Һ�ڵ�λ�����ϵ�Ħ��ѹ��
end