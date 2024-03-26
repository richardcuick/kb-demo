function [tau_y,K,n]=Rheology_HerschelBulkley(data_rheology)
    % 根据测试温度、压力条件下的赫巴流体六速粘度计数据，采用回归法计算得到剪切应力tau_y、稠度系数K和流性指数n
    % 计算方法：使用最小二乘法进行回归计算
    % 输入参数：data_rheology为测试温度、压力条件下的赫巴流体六速粘度计数据，第一列为转速，第二列为读数
    % 输出参数：tau_y为剪切应力（Pa），K为稠度系数（Pa·s^n），n为流性指数

    N=data_rheology(:,1); % N为六速粘度计的转速
    theta=data_rheology(:,2); % theta为六速粘度计的读数

    m=length(N); % m为六速粘度计数据的组数，此情况考虑到现场无法输入六速全部读数时，钻井液流变参数仍然可以使用回归方法进行计算

    %% 试算求解赫巴流变参数
    steps=floor(0.50771*theta(1)/0.001); % steps为试算次数
    for j=1:steps
        tau_y_trial(j)=0.001*j; % 不同试算次数下的屈服应力试算值
        for i=1:m % 计算不同转速下的实际剪切速率和实际剪切应力的自然对数值，log为自然对数
            gamma_dot_hrreg_trial(i)=log(1.7023*N(i)); % gamma_dot_hrreg_trial为实际剪切速率（s^-1）的自然对数
            tau_hrreg_trial(i)=log(0.50771*theta(i)-tau_y_trial(j)); % tau_hrreg_trial为实际剪切应力（Pa）的自然对数
        end

        % 使用线性回归计算屈服应力试算值tau_y_trial对应的稠度系数K_trial和流性指数n_trial
        % 公式参考：吴翊, 李永生, 胡庆军. 应用数理统计[M]. 国防科技大学出版社, 1995.
        x=gamma_dot_hrreg_trial; % 自变量x
        y=tau_hrreg_trial; % 因变量y

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
        
        % 将线性回归计算得到的截距和斜率分别赋值给稠度系数和流性指数
        K_trial(j)=exp(intercept); % 将e^截距的赋值给稠度系数
        n_trial(j)=slope; % 将斜率赋值给流性指数

        % 计算平均相对误差mre_h_trial
        for i=1:m
            gamma_dot_hrreg(i)=1.7023*N(i); % 计算实际剪切速率
            tau_hrreg(i)=0.50771*theta(i); % 计算实际剪切应力
            tau_hrreg_model(i)=tau_y_trial(j)+K_trial(j)*(gamma_dot_hrreg(i))^n_trial(j); % 计算K_trial和n_trial计算得到的屈服应力
        end

        mre_h_trial=0; % 平均相对误差
        for i=1:m
            re_h_trial=abs((tau_hrreg_model(i)-tau_hrreg(i))/tau_hrreg(i)); % 第i个数据组的相对误差，模型预测屈服应力与实际屈服应力之间的相对误差
            mre_h_trial=mre_h_trial+re_h_trial; % 相对误差求和
        end
        mre_h_trial=mre_h_trial/m; % 平均相对误差=相对误差之和/组数

        mre_h_trial_total(j)=mre_h_trial; % 存储每一次试算时的平均相对误差

        if n_trial(j)>1 % 如果流性指数大于1，则跳出本次循环，进行下一次循环
            mre_h_trial_total(j)=1; % 流性指数大于1时，平均相对误差赋值为1，同时跳出本次循环
            continue;
        end
    end

    %% 确定平均相对误差mre_h_trial_total的最小值及其对应的j
    j_min=1; % 存储相对误差最小值的位置，假设第一个数值为平均相对误差最小值
    mre_h_trial_total_min=mre_h_trial_total(1); % 假设第一个数值为平均相对误差最小值
    for i=2:steps % 搜索相对误差最小值
        if mre_h_trial_total(i)<mre_h_trial_total_min
            mre_h_trial_total_min=mre_h_trial_total(i); % 如果某处相对误差小于存储的相对误差最小值，则替换相对误差最小值
            j_min=i; % 存储相对误差最小值的位置
        end
    end

    %% 输出最小平均相对误差时的流变参数为最终的赫巴流变参数
    tau_y=tau_y_trial(j_min);
    K=K_trial(j_min);
    n=n_trial(j_min);
end