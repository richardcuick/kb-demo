function [K,n]=Rheology_PowerLaw(data_rheology)
    % 根据测试温度、压力条件下的幂律流体六速粘度计数据，采用回归法计算得到稠度系数K和流性指数n
    % 计算方法：使用最小二乘法进行回归计算
    % 输入参数：data_rheology为测试温度、压力条件下的幂律流体六速粘度计数据，第一列为转速，第二列为读数
    % 输出参数：K为稠度系数（Pa·s^n），n为流性指数
    
    N=data_rheology(:,1); % N为六速粘度计的转速
    theta=data_rheology(:,2); % theta为六速粘度计的读数
    
    m=length(N); % m为六速粘度计数据的组数，此情况考虑到现场无法输入六速全部读数时，钻井液流变参数仍然可以使用回归方法进行计算
    
    for i=1:m % 计算不同转速下的实际剪切速率和实际剪切应力的自然对数值，log为自然对数
        gamma_dot_brreg(i)=log(1.7023*N(i)); % gamma_dot_brreg为实际剪切速率（s^-1）的自然对数
        tau_brreg(i)=log(0.50771*theta(i)); % tau_brreg为实际剪切应力（Pa）的自然对数
    end
    
    %% 使用线性回归计算稠度系数K和流性指数n
    % 公式参考：吴翊, 李永生, 胡庆军. 应用数理统计[M]. 国防科技大学出版社, 1995.
    x=gamma_dot_brreg; % 自变量x
    y=tau_brreg; % 因变量y
    
    x_sum=0; % x求和值
    y_sum=0; % y求和值
    for i=1:m
        x_sum=x_sum+x(i); % 对自变量x求和
        y_sum=y_sum+y(i); % 对因变量y求和
    end
    
    x_ave=x_sum/m; % 对自变量x求平均值
    y_ave=y_sum/m; % 对因变量y求平均值
    
    L_xx=0; % (x-x_ave)^2求和值，即各个自变量与自变量平均值间差值平方和
    L_yy=0; % (y-y_ave)^2求和值，即各个因变量与因变量平均值间差值平方和
    L_xy=0; % (x-x_ave)*(y-y_ave)求和值
    for i=1:m
        L_xx=L_xx+(x(i)-x_ave)^2; % 各个自变量与自变量平均值间差值平方和
        L_yy=L_yy+(y(i)-y_ave)^2; % 各个因变量与因变量平均值间差值平方和
        L_xy=L_xy+(x(i)-x_ave)*(y(i)-y_ave); % (x-x_ave)*(y-y_ave)求和值
    end
    
    slope=L_xy/L_xx; % 计算斜率slope
    intercept=y_ave-(L_xy/L_xx)*x_ave; % 计算截距intercept
    r=L_xy/((L_xx*L_yy)^0.5); % 计算相关系数r
    R2=r^2; % 决定系数/拟合优度
    
    %% 将线性回归计算得到的截距和斜率分别赋值给稠度系数和流性指数
    K=exp(intercept); % 将e^截距的赋值给稠度系数
    n=slope; % 将斜率赋值给流性指数
end