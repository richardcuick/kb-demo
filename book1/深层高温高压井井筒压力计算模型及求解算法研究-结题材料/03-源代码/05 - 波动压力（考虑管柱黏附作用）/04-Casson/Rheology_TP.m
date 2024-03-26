function [tau_y_c,eta_infinity]=Rheology_TP(tau_y_c_0,eta_infinity_0,T_0,P_0,T,P)
    % 根据测试温度、压力下的卡森屈服应力、卡森粘度，计算目标温度、压力下的卡森屈服应力、卡森粘度
    % 输入参数：tau_y_c_0为测试温度、压力下的卡森屈服应力（Pa），eta_infinity_0为测试温度、压力下的卡森粘度（Pa*s），T_0为测试温度（℃），P_0为测试压力（Pa），T为目标温度（℃），P为目标压力（Pa）
    % 输出参数：tau_y_c为目标温度、压力下的卡森屈服应力（Pa），eta_infinity为目标温度、压力下的卡森粘度（Pa*s）
    
    alpha_tau=-0.01094;
    beta_tau=7.51953*10^(-10);
%     tau_y_c=tau_y_c_0*exp(alpha_tau*(T-T_0)+beta_tau*(P-P_0));
    tau_y_c=tau_y_c_0;
    
    alpha_eta=-0.01974;
    beta_eta=0.01233*10^(-6);
%     eta_infinity=eta_infinity_0*exp(alpha_eta*(T-T_0)+beta_eta*(P-P_0));
    eta_infinity=eta_infinity_0;
end