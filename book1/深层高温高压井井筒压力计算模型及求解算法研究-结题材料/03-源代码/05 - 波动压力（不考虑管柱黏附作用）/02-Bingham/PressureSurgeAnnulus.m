function [PressureSurge_a,flow_pattern_a]=PressureSurgeAnnulus(rho_a,V_a,mu_p_a,tau_y,D_w,D_d_o)
    % 计算波动压力
    % 输入参数：rho_a为环空钻井液密度（kg/m^3），V_a为环空钻井液流速（m/s），mu_p_a为环空钻井液的塑性粘度（Pa·s），tau_y为环空钻井液屈服应力（Pa），D_w为井筒内径（m），D_d_o为管柱外径（m），Qm为钻井液质量流量（kg/s）
    % 输出参数：PressureSurge_a为环空中单位长度上的波动压力（Pa/m），flow_pattern_a为环空流型
    
    %% 使用迭代法计算壁面处的屈服应力
    tau_w_old=100; % 假定壁面处屈服应力的初始值
    err_tau=1; % 给定壁面处屈服应力的初始误差，使其进入迭代循环
    while abs(err_tau)>1e-4 % 壁面处屈服应力收敛的判断标准
        tau_w_new=8*mu_p_a*(V_a*(0.25*pi*(D_w^2-D_d_o^2)))/(pi*((D_w-D_d_o)/2)^2*(D_w+D_d_o)/2)+3/2*tau_y-1/2*tau_y^3/tau_w_old^2; % 计算新的壁面屈服应力值
        err_tau=abs(tau_w_new-tau_w_old)/tau_w_old; % 计算屈服应力误差值
        tau_w_old=tau_w_new; % 将新计算的壁面处屈服应力值赋值给假设的壁面处屈服应力值
    end
    
    tau_w=tau_w_new; % 得到壁面屈服应力的终值
    
    int1=tau_y/tau_w; % int1为屈服应力和壁面屈服应力的比值，方便后续计算
    He=16800*int1/(1-int1)^3; % 计算Hedstrom数
    
    %% 根据环空中不同的流型，计算环空中钻井液在单位长度上的波动压力
    Re_a_c=(1-4/3*int1+1/3*int1^4)/(8*int1)*He; % Re_a_c为环空流体临界雷诺数
    Re_a=0.81619*rho_a*V_a*(D_w-D_d_o)/mu_p_a; % Re_a为环空流体雷诺数
    
    if Re_a==0 % 流体静止
        PressureSurge_a=0; % 计算流体静止条件下，环空钻井液单位长度上的波动压力，Pa/m
        flow_pattern_a=0; % 环空流型，0表示流体静止
    elseif Re_a>0 && Re_a<=Re_a_c % 层流
        PressureSurge_a=6895*mu_p_a*V_a/(216*(D_w-D_d_o)^2)+5.33355*tau_y/(D_w-D_d_o); % 计算层流条件下，环空钻井液单位长度上的波动压力，Pa/m
        flow_pattern_a=1; % flow_pattern_a为环空流型，1表示层流
    else % 湍流
        PressureSurge_a=0.158278*rho_a^0.75*V_a^1.75*mu_p_a^0.25/(D_w-D_d_o)^1.25; % 计算湍流条件下，环空钻井液单位长度上的波动压力，Pa/m
        flow_pattern_a=3; % flow_pattern_a为环空流型，3表示湍流
    end
end