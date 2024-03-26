function [PressureSurge_a,flow_pattern_a]=PressureSurgeAnnulus(rho_a,V_a,K_a,n_a,D_w,D_d_o)
    % 计算波动压力
    % 输入参数：rho_a为环空钻井液密度（kg/m^3），V_a为环空钻井液流速（m/s），K_a为环空钻井液稠度系数（Pa·s^n），n_a为环空钻井液流性指数，D_w为井筒内径（m），D_d_o为管柱外径（m）
    % 输出参数：PressureSurge_a为环空中单位长度上的波动压力（Pa/m），flow_pattern_a为环空流型
    
    Re_a=12^(1-n_a)*rho_a*(D_w-D_d_o)^n_a*V_a^(2-n_a)/(K_a*((2*n_a+1)/(3*n_a))^n_a); % 环空雷诺数
    
    Re_c_l=3470-1370*n_a; % 层流临界雷诺数
    Re_c_t=4270-1370*n_a; % 湍流临界雷诺数
    
    if Re_a==0 % 流体静止
        f_a=0; % 范宁摩擦因子
        PressureSurge_a=0; % 单位长度上的波动压力，Pa/m
        flow_pattern_a=0; % 环空流型，0表示流体静止
    elseif Re_a>0 && Re_a<=Re_c_l % 层流
        f_a=24/Re_a; % 范宁摩擦因子
        PressureSurge_a=2*f_a*rho_a*V_a^2/(D_w-D_d_o); % 单位长度上的波动压力，Pa/m
        flow_pattern_a=1; % flow_pattern_a为环空流型，1表示层流
    elseif Re_a>=Re_c_t % 湍流
        a=0.02*log10(n_a)+0.0786; % 系数
        b=0.25-0.143*log10(n_a); % 系数
        f_a=a/Re_a^b; % 范宁摩擦因子
        PressureSurge_a=2*f_a*rho_a*V_a^2/(D_w-D_d_o); % 单位长度上的波动压力，Pa/m
        flow_pattern_a=3; % flow_pattern_a为环空流型，3表示湍流
    else % 过渡流
        a=0.02*log10(n_a)+0.0786; % 系数
        b=0.25-0.143*log10(n_a); % 系数
        f_c_l=24/Re_c_l; % 临界层流范宁摩擦因子
        f_c_t=a/Re_c_t^b; % 临界湍流范宁摩擦因子
        f_a=f_c_l+(Re_a-Re_c_l)/(Re_c_t-Re_c_l)*(f_c_t-f_c_l); % 范宁摩擦因子
        PressureSurge_a=2*f_a*rho_a*V_a^2/(D_w-D_d_o); % 单位长度上的波动压力，Pa/m
        flow_pattern_a=2; % flow_pattern_a为环空流型，2表示过渡流
    end
end