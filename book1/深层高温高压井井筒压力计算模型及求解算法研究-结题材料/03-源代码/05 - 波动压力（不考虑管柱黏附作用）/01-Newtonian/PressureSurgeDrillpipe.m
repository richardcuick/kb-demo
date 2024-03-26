function [PressureSurge_d,flow_pattern_d]=PressureSurgeDrillpipe(rho_d,V_d,mu_d,D_d_i)
    % 计算波动压力
    % 输入参数：rho_d为钻杆中钻井液密度（kg/m^3），V_d为钻杆中钻井液流速（m/s），mu_d为钻杆中钻井液的塑性黏度（Pa·s），D_d_i为钻杆内径（m）
    % 输出参数：PressureSurge_d为钻杆中单位长度上的波动压力（Pa/m），flow_pattern_d为钻杆中钻井液流型
    
    Re_d=rho_d*V_d*D_d_i/mu_d; % 计算钻杆雷诺数
    if Re_d==0 % 流体静止
        f_d=0; % 计算流体静止条件下范宁摩擦因子
        flow_pattern_d=0; % 环空流型，0表示流体静止
    elseif Re_d>0 && Re_d<=2100 % 层流
        f_d=16/Re_d; % 计算层流条件下范宁摩擦因子
        flow_pattern_d=1; % 管内流型，1表示层流
    else % 湍流
        f_d=0.0791/Re_d^0.25; % 计算湍流条件下范宁摩擦因子
        flow_pattern_d=3; % 管内流型，3表示湍流
    end
    PressureSurge_d=2*f_d*V_d^2*rho_d/D_d_i; % 计算钻杆内钻井液单位长度上的波动压力，Pa/m
end