function delta_P_bit=Bit_PressureDrop(C,d_nozzle,m_nozzle,QL,rho_0)
    % ������ͷ���������������ֱ������������ϵ����������ͷ��������������ͷѹ��
    % ���������CΪ��������ϵ����d_nozzleΪ����ֱ����m����m_nozzleΪ���������QLΪ������ͷ�����������m^3/s����rho_0Ϊ�꾮�ܶȣ�kg/m^3��
    % ���������delta_P_bitΪ��ͷѹ����Pa��
    
    g=9.81; % gΪ�������ٶȣ�m/s^2��
    V_b=QL/(m_nozzle*(pi/4*d_nozzle^2)); % V_bΪ�꾮Һ������������٣�m/s��
    delta_P_bit=rho_0*V_b^2/(0.2039*g*C^2); % delta_P_bitΪ��ͷѹ����Pa��
end