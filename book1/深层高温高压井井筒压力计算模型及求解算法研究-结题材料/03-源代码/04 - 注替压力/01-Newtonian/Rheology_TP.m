function mu=Rheology_TP(mu_0,T_0,P_0,T,P)
    % ���ݲ����¶ȡ�ѹ���Ͳ��������µ�ճ�ȣ�����Ŀ���¶ȡ�ѹ���µ�ճ��
    % ���������mu_0Ϊ�����¶ȡ�ѹ���µ�ճ�ȣ�Pa��s����T_0Ϊ�����¶ȣ��棩��P_0Ϊ����ѹ����Pa����TΪĿ���¶ȣ��棩��PΪĿ��ѹ����Pa��
    % ���������muΪĿ���¶ȡ�ѹ���µ�ճ�ȣ�Pa��s��
    
    alpha=-0.011; % ˮ���꾮Һϵ��ʹ��alpha��beta��gamma
    beta=6.4*10^(-6); % ˮ���꾮Һϵ��ʹ��alpha��beta��gamma
    gamma=0.736; % ˮ���꾮Һϵ��ʹ��alpha��beta��gamma
    mu=gamma*mu_0*exp(alpha*(T-T_0)+beta*(P-P_0)/6895); % ����Ŀ���¶ȡ�ѹ���µ�ճ�ȣ�Pa��s
    
%     alpha_mu=-0.016;
%     beta_mu=8.759*10^(-11);
%     gamma_mu=4.427*10^(-6);
% 
%     mu=mu_0*exp(alpha_mu*(T-T_0)+beta_mu*(P-P_0)+gamma_mu*(T-T_0)^2);
    
%     mu=mu_0;
end