function [tau_y,mu_p]=Rheology_Bingham(data_rheology)
    % 根据测试温度、压力条件下的宾汉流体六速粘度计数据，采用回归法计算得到剪切应力tau_y和塑性黏度mu_p
    % 计算方法：使用最小二乘法进行回归计算
    % 输入参数：data_rheology为测试温度、压力条件下的宾汉流体六速粘度计数据，第一列为转速，第二列为读数
    % 输出参数：tau_y为屈服应力（Pa），mu_p为塑性黏度（Pa·s）
    
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
    
    %% 将线性回归计算得到的截距和斜率分别赋值给屈服应力和塑性粘度
    tau_y=intercept; % 将截距赋值给屈服应力
    mu_p=slope; % 将斜率赋值给塑性粘度
end