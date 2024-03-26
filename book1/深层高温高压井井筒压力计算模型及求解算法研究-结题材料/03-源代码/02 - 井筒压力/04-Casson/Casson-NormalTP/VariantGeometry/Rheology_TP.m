function [tau_y_c,eta_infinity]=Rheology_TP(tau_y_c_0,eta_infinity_0,T_0,P_0,T,P)
    % 根据测试温度、压力下的卡森屈服应力、卡森粘度，计算目标温度、压力下的卡森屈服应力、卡森粘度
    % 输入参数：tau_y_c_0为测试温度、压力下的卡森屈服应力（Pa），eta_infinity_0为测试温度、压力下的卡森粘度（Pa*s），T_0为测试温度（℃），P_0为测试压力（Pa），T为目标温度（℃），P为目标压力（Pa）
    % 输出参数：tau_y_c为目标温度、压力下的卡森屈服应力（Pa），eta_infinity为目标温度、压力下的卡森粘度（Pa*s）

    tau_y_c=tau_y_c_0;
    
    eta_infinity=eta_infinity_0;
end