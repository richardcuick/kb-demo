function [delta_P_bit,V_nozzle]=PressureDrop_Bit(C,d_nozzle,n_nozzle,Qv,rho_0)
    % ������ͷ���������������ֱ������������ϵ����������ͷ��������������ͷѹ��
    % ���������CΪ��������ϵ����d_nozzleΪ����ֱ����m����m_nozzleΪ���������QvΪ������ͷ�����������m^3/s����rho_0Ϊ�꾮Һ�ܶȣ�kg/m^3��
    % ���������delta_P_bitΪ��ͷѹ����Pa��
    
    g=9.81; % gΪ�������ٶȣ�m/s^2��
    V_nozzle=Qv/(n_nozzle*(pi/4*d_nozzle^2)); % V_nozzleΪ�꾮Һ������������٣�m/s��
    delta_P_bit=rho_0*V_nozzle^2/(0.2039*g*C^2); % delta_P_bitΪ��ͷѹ����Pa��
end