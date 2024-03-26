function [PressureSurge_a,flow_pattern_a]=PressureSurgeAnnulus(rho_a,V_a,mu_p_a,tau_y,D_w,D_d_o)
    % ���㲨��ѹ��
    % ���������rho_aΪ�����꾮Һ�ܶȣ�kg/m^3����V_aΪ�����꾮Һ���٣�m/s����mu_p_aΪ�����꾮Һ������ճ�ȣ�Pa��s����tau_yΪ�����꾮Һ����Ӧ����Pa����D_wΪ��Ͳ�ھ���m����D_d_oΪ�����⾶��m����QmΪ�꾮Һ����������kg/s��
    % ���������PressureSurge_aΪ�����е�λ�����ϵĲ���ѹ����Pa/m����flow_pattern_aΪ��������
    
    %% ʹ�õ�����������洦������Ӧ��
    tau_w_old=100; % �ٶ����洦����Ӧ���ĳ�ʼֵ
    err_tau=1; % �������洦����Ӧ���ĳ�ʼ��ʹ��������ѭ��
    while abs(err_tau)>1e-4 % ���洦����Ӧ���������жϱ�׼
        tau_w_new=8*mu_p_a*(V_a*(0.25*pi*(D_w^2-D_d_o^2)))/(pi*((D_w-D_d_o)/2)^2*(D_w+D_d_o)/2)+3/2*tau_y-1/2*tau_y^3/tau_w_old^2; % �����µı�������Ӧ��ֵ
        err_tau=abs(tau_w_new-tau_w_old)/tau_w_old; % ��������Ӧ�����ֵ
        tau_w_old=tau_w_new; % ���¼���ı��洦����Ӧ��ֵ��ֵ������ı��洦����Ӧ��ֵ
    end
    
    tau_w=tau_w_new; % �õ���������Ӧ������ֵ
    
    int1=tau_y/tau_w; % int1Ϊ����Ӧ���ͱ�������Ӧ���ı�ֵ�������������
    He=16800*int1/(1-int1)^3; % ����Hedstrom��
    
    %% ���ݻ����в�ͬ�����ͣ����㻷�����꾮Һ�ڵ�λ�����ϵĲ���ѹ��
    Re_a_c=(1-4/3*int1+1/3*int1^4)/(8*int1)*He; % Re_a_cΪ���������ٽ���ŵ��
    Re_a=0.81619*rho_a*V_a*(D_w-D_d_o)/mu_p_a; % Re_aΪ����������ŵ��
    
    if Re_a==0 % ���徲ֹ
        PressureSurge_a=0; % �������徲ֹ�����£������꾮Һ��λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_a=0; % �������ͣ�0��ʾ���徲ֹ
    elseif Re_a>0 && Re_a<=Re_a_c % ����
        PressureSurge_a=6895*mu_p_a*V_a/(216*(D_w-D_d_o)^2)+5.33355*tau_y/(D_w-D_d_o); % ������������£������꾮Һ��λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_a=1; % flow_pattern_aΪ�������ͣ�1��ʾ����
    else % ����
        PressureSurge_a=0.158278*rho_a^0.75*V_a^1.75*mu_p_a^0.25/(D_w-D_d_o)^1.25; % �������������£������꾮Һ��λ�����ϵĲ���ѹ����Pa/m
        flow_pattern_a=3; % flow_pattern_aΪ�������ͣ�3��ʾ����
    end
end