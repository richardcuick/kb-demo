function [PressureSurge_d,flow_pattern_d]=PressureSurgeDrillPipe(rho_d,V_d,K_d,n_d,D_d_i)
    % ���㲨��ѹ��
    % ���������rho_dΪ�������꾮Һ�ܶȣ�kg/m^3����V_dΪ������꾮Һ���٣�m/s����K_dΪ������꾮Һ�ĳ��ϵ����Pa��s^n����n_dΪ������꾮Һ������ָ����D_d_iΪ����ھ���m��
    % ���������PressureSurge_dΪ����е�λ�����ϵĲ���ѹ����Pa/m����flow_pattern_dΪ������꾮Һ����
    
    Re_d=8^(1-n_d)*rho_d*D_d_i^n_d*V_d^(2-n_d)/(K_d*((3*n_d+1)/(4*n_d))^n_d); % �����ŵ��
    
    Re_c_l=3470-1370*n_d; % �����ٽ���ŵ��
    Re_c_t=4270-1370*n_d; % �����ٽ���ŵ��
    
    if Re_d==0 % ���徲ֹ
        f_d=0; % ����Ħ������
        PressureSurge_d=0; % ��λ���Ȳ���ѹ����Pa/m
        flow_pattern_d=0; % �������ͣ�0��ʾ���徲ֹ
    elseif Re_d>0 && Re_d<=Re_c_l % ����
        f_d=16/Re_d; % ����Ħ������
        PressureSurge_d=2*f_d*rho_d*V_d^2/D_d_i; % ��λ���Ȳ���ѹ����Pa/m
        flow_pattern_d=1; % ������ͣ�1��ʾ����
    elseif Re_d>=Re_c_t % ����
        a=0.02*log10(n_d)+0.0786; % ϵ��
        b=0.25-0.143*log10(n_d); % ϵ��
        f_d=a/Re_d^b; % ����Ħ������
        PressureSurge_d=2*f_d*rho_d*V_d^2/D_d_i; % ��λ���Ȳ���ѹ����Pa/m
        flow_pattern_d=3; % ������ͣ�3��ʾ����
    else % ������
        a=0.02*log10(n_d)+0.0786; % ϵ��
        b=0.25-0.143*log10(n_d); % ϵ��
        f_c_l=16/Re_c_l; % �ٽ��������Ħ������
        f_c_t=a/Re_c_t^b; % �ٽ���������Ħ������
        f_d=f_c_l+(Re_d-Re_c_l)/(Re_c_t-Re_c_l)*(f_c_t-f_c_l); % ����Ħ������
        PressureSurge_d=2*f_d*rho_d*V_d^2/D_d_i; % ��λ���Ȳ���ѹ����Pa/m
        flow_pattern_d=2; % ������ͣ�2��ʾ������
    end
end