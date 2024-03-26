function [PressureSurge_a,flow_pattern_a]=PressureSurgeAnnulus(rho_a,V_a,eta_infinity_a,tau_y_c_a,D_w,D_d_o)
    % 根据环空结构、钻井液性质和流动参数，计算波动压力
    % 输入参数：rho_a为环空钻井液密度（kg/m^3），V_a为环空钻井液流速（m/s），eta_infinity_a为环空钻井液卡森黏度（Pa·s），tau_y_c_a为环空钻井液卡森屈服应力（Pa），D_w为井筒内径（m），D_d_o为管柱外径（m），Qm为钻井液质量流量（kg/s）
    % 输出参数：PressureSurge_a为环空中单位长度上的波动压力（Pa/m），flow_pattern_a为环空流型
    
    Qv=V_a*(0.25*pi*(D_w^2-D_d_o^2)); % 钻井液体积流量，m^3/s
    D_h=D_w-D_d_o; % 水力直径，m
    
    %% 使用迭代法计算壁面处的剪切应力
    tau_w_old=100; % 假定壁面处剪切应力的初始值
    loop=1;
    err_tau=1;
    while abs(err_tau)>1e-4 && loop<=1000 % 迭代终止条件
        int=(6*Qv*eta_infinity_a)/(pi*(D_w/2+D_d_o/2)*(D_h/2)^2)-(3*tau_y_c_a)/2+(5*tau_y_c_a^3)/(2*tau_w_old^2)+(12*tau_y_c_a^0.5*(tau_w_old^2.5-tau_y_c_a^2.5))/(5*tau_w_old^2);
        tau_w_new=int;
        err_tau=abs(tau_w_old-tau_w_new)/tau_w_old; % 计算剪切应力误差值
        tau_w_old=tau_w_new;
        loop=loop+1;
    end
    tau_w=tau_w_new; % 壁面剪切应力，Pa
    
    gammadot_w=(tau_w-2*((tau_w*tau_y_c_a)^0.5)+tau_y_c_a)/eta_infinity_a; % 壁面剪切速率，s^-1
    
    nprime=(4*V_a/D_h)/(gammadot_w-8*V_a/D_h); % 广义流性指数
    
    Re_gcl_a=3470-1370*nprime; % 层流临界广义雷诺数
    Re_gct_a=4270-1370*nprime; % 湍流临界广义雷诺数
    
    D_eff=(2*nprime*D_h)/(2*nprime+1); % 有效管径
    mu_aw=tau_w/gammadot_w; % 壁面剪切应力（Pa）
    
    Re_g_a=(rho_a*D_eff*V_a)/mu_aw; % 广义雷诺数
    
    %% 根据环空中不同的流型，计算环空中钻井液在单位长度上的波动压力   
    if Re_g_a<=Re_gcl_a % 层流
        f_a=16/Re_g_a; % 计算层流条件下，环空内钻井液摩擦因子
        PressureSurge_a=2*f_a*rho_a*V_a^2/D_h; % 计算层流条件下，环空内钻井液在单位长度上的波动压力，Pa/m
        flow_pattern_a=1; % flow_pattern_a为环空流型，1表示层流
    elseif Re_g_a>Re_gct_a % 湍流
        f_a_old=100; % 迭代初值
        loop=1;
        err_tau=1;
        while abs(err_tau)>1e-4 && loop<=1000 % 迭代终止条件
            int=(4/nprime^0.75)*log10(Re_g_a*f_a_old^(1-nprime/2))-0.395/nprime^1.2;
            f_a_c_new=(int^-1)^2;
            err_tau=abs(f_a_old-f_a_c_new)/f_a_old;
            f_a_old=f_a_c_new;
            loop=loop+1;
        end
        f_a=f_a_c_new; % 计算湍流条件下，环空内钻井液摩擦因子
        PressureSurge_a=2*f_a*rho_a*V_a^2/D_h; % 计算湍流条件下，环空内钻井液在单位长度上的波动压力，Pa/m
        flow_pattern_a=3; % flow_pattern_a为环空流型，3表示湍流
    else
        Re_g_a=Re_gcl_a;
        f_a1=16/Re_g_a;
        
        Re_g_a=Re_gct_a;
        f_a_old=100; % 迭代初值
        loop=1;
        err_tau=1;
        while abs(err_tau)>1e-4&loop<=1000 % 迭代终止条件
            int=(4/nprime^0.75)*log10(Re_g_a*f_a_old^(1-nprime/2))-0.395/nprime^1.2;
            f_a_c_new=(int^-1)^2;
            err_tau=abs(f_a_old-f_a_c_new)/f_a_old;
            f_a_old=f_a_c_new;
            loop=loop+1;
        end
        f_a2=f_a_c_new;

        f_a=f_a1+((Re_g_a-Re_gcl_a)*(f_a2-f_a1))/(Re_gct_a-Re_gcl_a); % 计算过渡流条件下，环空内钻井液摩擦因子
        PressureSurge_a=2*f_a*rho_a*V_a^2/D_h; % 计算过渡流条件下，环空内钻井液在单位长度上的波动压力，Pa/m
        flow_pattern_a=2; % flow_pattern_a为环空流型，2表示过渡流
    end
end