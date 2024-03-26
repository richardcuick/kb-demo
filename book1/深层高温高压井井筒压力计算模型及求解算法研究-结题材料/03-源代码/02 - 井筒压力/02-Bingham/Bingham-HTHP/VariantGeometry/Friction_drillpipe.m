function [Ff_d,flow_pattern_d]=Friction_drillpipe(rho_d,V_d,mu_p_d,tau_y,D_d_i,ML)
    % 根据钻杆结构、钻井液性质和流动参数，计算钻杆摩擦压降及流型
    % 输入参数：rho_d为钻杆中钻井液密度（kg/m^3），V_d为钻杆中钻井液流速（m/s），mu_p_d为钻杆中钻井液的塑性黏度（Pa·s），tau_y为钻杆中钻井液的屈服应力（Pa），D_d_i为钻杆内径（m），ML为测试温度、压力下的钻井液质量排量（kg/s）
    % 输出参数：Ff_d为钻杆中单位长度上的摩擦压降（Pa/m），flow_pattern_d为钻杆中钻井液流型
    
    QL=ML/rho_d; % 钻井液体积排量=钻井液质量流量/钻井液密度
    
    %% 使用迭代法计算壁面处的屈服应力
    tau_w_old=100; % 假定壁面处屈服应力的初始值
    err_tau=1; % 给定壁面处屈服应力的初始误差，使其进入迭代循环
    while abs(err_tau)>1e-4 % 壁面处屈服应力收敛的判断标准
        tau_w_new=4*mu_p_d*QL/(pi*(D_d_i/2)^3)+4/3*tau_y-1/3*tau_y^4/tau_w_old^3; % 计算新的壁面屈服应力值
        err_tau=abs(tau_w_new-tau_w_old)/tau_w_old; % 计算屈服应力误差值
        tau_w_old=tau_w_new; % 将新计算的壁面处屈服应力值赋值给假设的壁面处屈服应力值
    end
    
    tau_w=tau_w_new; % 得到壁面屈服应力的终值
    
    int1=tau_y/tau_w; % int1为屈服应力和壁面屈服应力的比值，方便后续计算
    He=16800*int1/(1-int1)^3; % 计算Hedstrom数
    
    %% 根据钻杆中不同的流型，计算钻杆中钻井液在单位长度上的摩擦压降
    % 雷诺数<=临界雷诺数，为层流；雷诺数>临界雷诺数，为湍流
    % 对于宾汉流型，不划分过渡流，仅划分为层流及湍流
    Re_d_c=(1-4/3*int1+1/3*int1^4)/(8*int1)*He; % Re_d_c为临界钻杆雷诺数
    Re_d=rho_d*V_d*D_d_i/mu_p_d; % Re_d为钻杆雷诺数
    
    Ff_d1=6895*mu_p_d*V_d/(216*D_d_i^2)+5.33355*tau_y/D_d_i; % 计算层流条件下，钻杆内钻井液在单位长度上的摩擦压降
    Ff_d2=0.158278*rho_d^0.75*V_d^1.75*mu_p_d^0.25/D_d_i^1.25; % 计算湍流条件下，钻杆内钻井液在单位长度上的摩擦压降

    if Ff_d1>Ff_d2
        flow_pattern_d=1; % flow_pattern_d为钻杆流型，1表示层流
        Ff_d=Ff_d1;
    else
        flow_pattern_d=3; % flow_pattern_d为钻杆流型，3表示湍流
        Ff_d=Ff_d2;
    end
end