function [PressureSurge_a,flow_pattern_a]=PressureSurgeAnnulus(rho_a,V_a,mu_a,D_w,D_d_o)
    % 计算波动压力
    % 输入参数：rho_a为环空中钻井液密度（kg/m^3），V_a为环空中钻井液流速（m/s），mu_a为环空中钻井液的黏度（Pa·s），D_w为井筒内径（m），D_d_o为管柱外径（m）
    % 输出参数：PressureSurge_a为环空中单位长度上的波动压力（Pa/m），flow_pattern_a为环空流型
    
    Re_a=0.816*rho_a*V_a*(D_w-D_d_o)/mu_a; % 计算环空雷诺数
    if Re_a==0 % 流体静止
        f_a=0; % 计算流体静止条件下范宁摩擦因子
        flow_pattern_a=0; % 环空流型，0表示流体静止
    elseif Re_a>0 && Re_a<=2100 % 层流
        f_a=16/Re_a; % 计算层流条件下范宁摩擦因子
        flow_pattern_a=1; % 环空流型，1表示层流
    else % 湍流
        f_a=0.0791/Re_a^0.25; % 计算湍流条件下范宁摩擦因子
        flow_pattern_a=3; % 环空流型，3表示湍流
    end
    PressureSurge_a=2*f_a*V_a^2*rho_a/(0.816*(D_w-D_d_o)); % 计算环空内钻井液单位长度上的波动压力，Pa/m
end