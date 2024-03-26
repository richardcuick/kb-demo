function [Ff_d,flow_pattern_d]=Friction_drillpipe(rho_d,V_d,eta_infinity_d,tau_y_c_d,D_di,ML)
    % 根据钻杆结构、钻井液性质和流动参数，计算钻杆摩擦压降及流型
    % 输入参数：rho_d为钻杆中钻井液密度（kg/m^3），V_d为钻杆中钻井液流速（m/s），A_d为钻杆内截面积（m^2），eta_infinity_d为钻杆中钻井液的卡森黏度（Pa·s），tau_y_c_d为钻杆中钻井液的卡森屈服应力（Pa），D_di为钻杆内径（m），ML为测试温度、压力下的钻井液质量排量（kg/s）
    % 输出参数：Ff_d为钻杆中单位长度上的摩擦压降（Pa/m），flow_pattern_d为钻杆中钻井液流型
    
    QL=ML/rho_d; % 钻井液体积排量=钻井液质量流量/钻井液密度
    
    %% 使用迭代法计算壁面处的屈服应力
    tau_w_old=100; % 迭代初值
    loop=1;
    err_tau=1;
    while abs(err_tau)>1e-4&&loop<=1000 % 迭代终止条件
        int=(4*QL*eta_infinity_d)/(pi*(D_di/2)^3)-(4*tau_y_c_d)/3+(16*(tau_y_c_d*tau_w_old)^0.5)/7+(tau_y_c_d^4)/(21*tau_w_old^3);
        tau_w_new=int;
        err_tau=abs(tau_w_old-tau_w_new)/tau_w_old;
        tau_w_old=tau_w_new;
        loop=loop+1;
    end
    tau_w=tau_w_new; % 壁面剪切应力，Pa
    
    gammadot_w=(tau_w-2*((tau_w*tau_y_c_d)^0.5)+tau_y_c_d)/eta_infinity_d; % 壁面剪切速率，s^-1
    
    nprime=(8*V_d/D_di)/(4*gammadot_w-24*V_d/D_di); % 广义流性指数
    
    Re_gcl_d=3470-1370*nprime; % 层流临界广义雷诺数
    Re_gct_d=4270-1370*nprime; % 湍流临界广义雷诺数
    
    D_eff=(4*nprime*D_di)/(3*nprime+1); % 有效管径
    mu_aw=tau_w/gammadot_w; % 壁面剪切速率，s^-1
    
    Re_g_d=(rho_d*D_eff*V_d)/mu_aw; % 广义雷诺数
    
    %% 根据钻杆中不同的流型，计算钻杆中钻井液在单位长度上的摩擦压降
    % 广义雷诺数<=层流临界广义雷诺数，为层流；
    % 广义雷诺数>湍流临界广义雷诺数，为湍流；
    % 层流临界广义雷诺数<广义雷诺数<=湍流临界广义雷诺数，为过渡流；
    
    if Re_g_d<=Re_gcl_d % 层流
        f_d=16/Re_g_d; % 计算层流条件下，钻杆内钻井液摩擦因子
        Ff_d=2*f_d*rho_d*V_d^2/D_di; % 计算层流条件下，钻杆内钻井液在单位长度上的摩擦压降
        flow_pattern_d=1; % flow_pattern_d为钻杆流型，1表示层流
    elseif Re_g_d>Re_gct_d % 湍流
        f_d_old=100; % 迭代初值
        err_f=1;
        loop=1;
        while abs(err_f)>1e-4&&loop<=1000 % 迭代终止条件
            int=(4/nprime^0.75)*log10(Re_g_d*f_d_old^(1-nprime/2))-0.395/nprime^1.2;
            f_d_new=(int^-1)^2;
            err_f=abs(f_d_old-f_d_new)/f_d_old;
            f_d_old=f_d_new;
            loop=loop+1;
        end
        f_d=f_d_new; % 壁面剪切应力，Pa
        Ff_d=2*f_d*rho_d*V_d^2/D_di; % 计算层流条件下，钻杆内钻井液在单位长度上的摩擦压降
        flow_pattern_d=3;
    else
        Re_g_d=Re_gcl_d;
        f_d1=16/Re_g_d;
        
        Re_g_d=Re_gct_d;
        f_d_old=100; % 迭代初值
        loop=1;
        err_f=1;
        while abs(err_f)>1e-4&&loop<=1000 % 迭代终止条件
            int=(4/nprime^0.75)*log10(Re_g_d*f_d_old^(1-nprime/2))-0.395/nprime^1.2;
            f_d_new=(int^-1)^2;
            err_f=abs(f_d_old-f_d_new)/f_d_old;
            f_d_old=f_d_new;
            loop=loop+1;
        end
        f_d2=f_d_new; % 壁面剪切应力，Pa
        
        f_d=f_d1+((Re_g_d-Re_gcl_d)*(f_d2-f_d1))/(Re_gct_d-Re_gcl_d);
        Ff_d=2*f_d*rho_d*V_d^2/D_di; % 计算层流条件下，钻杆内钻井液在单位长度上的摩擦压降
        flow_pattern_d=2;
    end
end