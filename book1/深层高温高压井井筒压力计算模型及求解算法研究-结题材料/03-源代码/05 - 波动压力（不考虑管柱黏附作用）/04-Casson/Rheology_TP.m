function [tau_y_c,eta_infinity]=Rheology_TP(tau_y_c_0,eta_infinity_0,T_0,P_0,T,P)
    % ���ݲ����¶ȡ�ѹ���µĿ�ɭ����Ӧ������ɭճ�ȣ�����Ŀ���¶ȡ�ѹ���µĿ�ɭ����Ӧ������ɭճ��
    % ���������tau_y_c_0Ϊ�����¶ȡ�ѹ���µĿ�ɭ����Ӧ����Pa����eta_infinity_0Ϊ�����¶ȡ�ѹ���µĿ�ɭճ�ȣ�Pa*s����T_0Ϊ�����¶ȣ��棩��P_0Ϊ����ѹ����Pa����TΪĿ���¶ȣ��棩��PΪĿ��ѹ����Pa��
    % ���������tau_y_cΪĿ���¶ȡ�ѹ���µĿ�ɭ����Ӧ����Pa����eta_infinityΪĿ���¶ȡ�ѹ���µĿ�ɭճ�ȣ�Pa*s��
    
    alpha_tau=-0.01094;
    beta_tau=7.51953*10^(-10);
%     tau_y_c=tau_y_c_0*exp(alpha_tau*(T-T_0)+beta_tau*(P-P_0));
    tau_y_c=tau_y_c_0;
    
    alpha_eta=-0.01974;
    beta_eta=0.01233*10^(-6);
%     eta_infinity=eta_infinity_0*exp(alpha_eta*(T-T_0)+beta_eta*(P-P_0));
    eta_infinity=eta_infinity_0;
end