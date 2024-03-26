function mu=Rheology_TP(mu_0,T_0,P_0,T,P)
    % 根据测试温度、压力和测试条件下的黏度，计算目标温度、压力下的黏度
    % 输入参数：mu_0为测试温度、压力下的黏度（Pa·s），T_0为测试温度（℃），P_0为测试压力（Pa），T为目标温度（℃），P为目标压力（Pa）
    % 输出参数：mu为目标温度、压力下的黏度（Pa·s）
    
    alpha=-0.011; % 水基钻井液系数使用alpha，beta和gamma
    beta=6.4*10^(-6); % 水基钻井液系数使用alpha，beta和gamma
    gamma=0.736; % 水基钻井液系数使用alpha，beta和gamma
    mu=gamma*mu_0*exp(alpha*(T-T_0)+beta*(P-P_0)/6895); % 计算目标温度、压力下的塑性黏度
end