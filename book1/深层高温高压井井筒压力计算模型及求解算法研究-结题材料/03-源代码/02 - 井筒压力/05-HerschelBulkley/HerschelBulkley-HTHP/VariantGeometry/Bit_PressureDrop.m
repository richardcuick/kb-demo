function delta_P_bit=Bit_PressureDrop(C,d_nozzle,m_nozzle,QL,rho_0)
    % 根据钻头上喷嘴个数、喷嘴直径、喷嘴流量系数及流经钻头的排量，计算钻头压降
    % 输入参数：C为喷嘴流量系数，d_nozzle为喷嘴直径（m），m_nozzle为喷嘴个数，QL为流经钻头的体积排量（m^3/s），rho_0为钻井密度（kg/m^3）
    % 输出参数：delta_P_bit为钻头压降（Pa）
    
    g=9.81; % g为重力加速度（m/s^2）
    V_b=QL/(m_nozzle*(pi/4*d_nozzle^2)); % V_b为钻井液流经喷嘴的流速（m/s）
    delta_P_bit=rho_0*V_b^2/(0.2039*g*C^2); % delta_P_bit为钻头压降（Pa）
end