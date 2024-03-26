function [Ff_a,flow_pattern_a]=Friction_annulus(rho_a,V_a,tau_y,mu_p_a,D_w,D_d_o,Qv_a)
    % 根据环空结构、钻井液性质和流动参数，计算环空摩擦压降及流型
    % 输入参数：rho_a为环空中钻井液密度（kg/m^3），V_a为环空中钻井液流速（m/s），mu_p_a为环空中钻井液的塑性黏度（Pa·s），tau_y为环空中钻井液的屈服应力（Pa），D_w为井筒内径（m），D_d_o为钻杆外径（m），Qv_a为钻井液体积流量（m^3/s）
    % 输出参数：Ff_a为环空中单位长度上的摩擦压降（Pa/m），flow_pattern_a为环空中钻井液流型
    
    %% 使用迭代法计算壁面处的屈服应力
    tau_w_old=100; % 假定壁面处屈服应力的初始值
    err_tau=1; % 给定壁面处屈服应力的初始误差，使其进入迭代循环
    while abs(err_tau)>1e-4 % 壁面处屈服应力收敛的判断标准
        tau_w_new=8*mu_p_a*Qv_a/(pi*((D_w-D_d_o)/2)^2*(D_w+D_d_o)/2)+3/2*tau_y-1/2*tau_y^3/tau_w_old^2; % 计算新的壁面屈服应力值
        err_tau=abs(tau_w_new-tau_w_old)/tau_w_old; % 计算屈服应力误差值
        tau_w_old=tau_w_new; % 将新计算的壁面处屈服应力值赋值给假设的壁面处屈服应力值
    end
    
    tau_w=tau_w_new; % 得到壁面屈服应力的终值
    
    int1=tau_y/tau_w; % int1为屈服应力和壁面屈服应力的比值，方便后续计算
    He=16800*int1/(1-int1)^3; % 计算Hedstrom数
    
    %% 根据环空中不同的流型，计算环空中钻井液在单位长度上的摩擦压降
    % 雷诺数<=临界雷诺数，为层流；雷诺数>临界雷诺数，为湍流
    % 对于宾汉流型，不划分过渡流，仅划分为层流及湍流
    Re_a_c=(1-4/3*int1+1/3*int1^4)/(8*int1)*He; % Re_a_c为临界环空雷诺数
    Re_a=0.81619*rho_a*V_a*(D_w-D_d_o)/mu_p_a; % Re_a为环空雷诺数
    
    if Re_a<=Re_a_c % 层流
        Ff_a=6895*mu_p_a*V_a/(216*(D_w-D_d_o)^2)+5.33355*tau_y/(D_w-D_d_o); % 计算层流条件下，环空内钻井液在单位长度上的摩擦压降
        flow_pattern_a=1; % flow_pattern_a为环空流型，1表示层流
    else % 湍流
        Ff_a=0.158278*rho_a^0.75*V_a^1.75*mu_p_a^0.25/(D_w-D_d_o)^1.25; % 计算湍流条件下，环空内钻井液在单位长度上的摩擦压降
        flow_pattern_a=3; % flow_pattern_a为环空流型，3表示湍流
    end
end