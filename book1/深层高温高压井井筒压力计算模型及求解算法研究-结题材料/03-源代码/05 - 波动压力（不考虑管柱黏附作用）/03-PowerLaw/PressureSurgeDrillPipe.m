function [PressureSurge_d,flow_pattern_d]=PressureSurgeDrillPipe(rho_d,V_d,K_d,n_d,D_d_i)
    % 计算波动压力
    % 输入参数：rho_d为环空中钻井液密度（kg/m^3），V_d为钻杆中钻井液流速（m/s），K_d为钻杆中钻井液的稠度系数（Pa·s^n），n_d为钻杆中钻井液的流性指数，D_d_i为钻杆内径（m）
    % 输出参数：PressureSurge_d为钻杆中单位长度上的波动压力（Pa/m），flow_pattern_d为钻杆中钻井液流型
    
    Re_d=8^(1-n_d)*rho_d*D_d_i^n_d*V_d^(2-n_d)/(K_d*((3*n_d+1)/(4*n_d))^n_d); % 钻杆雷诺数
    
    Re_c_l=3470-1370*n_d; % 层流临界雷诺数
    Re_c_t=4270-1370*n_d; % 湍流临界雷诺数
    
    if Re_d==0 % 流体静止
        f_d=0; % 范宁摩擦因子
        PressureSurge_d=0; % 单位长度波动压力，Pa/m
        flow_pattern_d=0; % 环空流型，0表示流体静止
    elseif Re_d>0 && Re_d<=Re_c_l % 层流
        f_d=16/Re_d; % 范宁摩擦因子
        PressureSurge_d=2*f_d*rho_d*V_d^2/D_d_i; % 单位长度波动压力，Pa/m
        flow_pattern_d=1; % 钻杆流型，1表示层流
    elseif Re_d>=Re_c_t % 湍流
        a=0.02*log10(n_d)+0.0786; % 系数
        b=0.25-0.143*log10(n_d); % 系数
        f_d=a/Re_d^b; % 范宁摩擦因子
        PressureSurge_d=2*f_d*rho_d*V_d^2/D_d_i; % 单位长度波动压力，Pa/m
        flow_pattern_d=3; % 钻杆流型，3表示湍流
    else % 过渡流
        a=0.02*log10(n_d)+0.0786; % 系数
        b=0.25-0.143*log10(n_d); % 系数
        f_c_l=16/Re_c_l; % 临界层流范宁摩擦因子
        f_c_t=a/Re_c_t^b; % 临界湍流范宁摩擦因子
        f_d=f_c_l+(Re_d-Re_c_l)/(Re_c_t-Re_c_l)*(f_c_t-f_c_l); % 范宁摩擦因子
        PressureSurge_d=2*f_d*rho_d*V_d^2/D_d_i; % 单位长度波动压力，Pa/m
        flow_pattern_d=2; % 钻杆流型，2表示过渡流
    end
end