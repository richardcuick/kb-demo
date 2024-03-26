function mu=Rheology_Newtonian(data_rheology)
    % 根据测试温度、压力条件下的牛顿流体六速粘度计数据，采用回归法计算得到黏度mu
    % 计算方法：使用最小二乘法进行回归计算
    % 输入参数：data_rheology为测试温度、压力条件下的牛顿流体六速粘度计数据，第一列为转速，第二列为读数
    % 输出参数：mu为黏度（Pa·s）
    
    N=data_rheology(:,1); % N为六速粘度计的转速
    theta=data_rheology(:,2); % theta为六速粘度计的读数
    
    m=length(N); % m为六速粘度计数据的组数，此情况考虑到现场无法输入六速全部读数时，钻井液流变参数仍然可以使用回归方法进行计算
    
    for i=1:m % 计算不同转速下的实际剪切速率和实际剪切应力
        gamma_dot_brreg(i)=1.7023*N(i); % gamma_dot_brreg为实际剪切速率（s^-1）
        tau_brreg(i)=0.50771*theta(i); % tau_brreg为实际剪切应力（Pa）
    end
    
    %% 使用线性回归计算屈服应力tau_y和塑性粘度mu_p
    % 公式参考：吴翊, 李永生, 胡庆军. 应用数理统计[M]. 国防科技大学出版社, 1995.
    x=gamma_dot_brreg; % 自变量x
    y=tau_brreg; % 因变量y
    
    L_xx=0; % x^2求和值
    L_xy=0; % x*y求和值
    for i=1:m
        L_xx=L_xx+x(i)^2; % 对x^2求和
        L_xy=L_xy+x(i)*y(i); % 对x*y求和
    end
    
    slope=L_xy/L_xx; % 计算斜率slope
    
    %% 将线性回归计算得到的斜率赋值给粘度
    mu=slope; % 将斜率赋值给粘度
    
end