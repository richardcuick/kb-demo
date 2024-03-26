function [result] = Function_ParameterSol11(data)
    % 输入：数据data存储不同温度、压力及相应流体参数（共3列），其中第一列为温度（℃），第二列为压力（MPa），最后一列为相应流体参数
    % 输出：result依次存储参考温度（℃）、参考压力（MPa）、常数C、系数ALPHA、系数BETA、系数GAMMA、系数DELTA、系数EPSILON、决定系数、平均相对误差
    
    n=size(data); % 数据（自变量+因变量）共多少行多少列
    n_row=n(1,1); % 样本观测次数
    n_col=n(1,2); % 自变量+因变量个数，因变量只有一个，所以自变量有（n_col-1）个

    %% 以温度最小值且压力最小值作为模型基准参考值
    min_T=data(1,1);
    for i=1:1:n_row
        if min_T>=data(i,1)
            min_T=data(i,1); % 温度最小值，℃
        else
        end
    end

    for i=1:1:n_row
        if min_T==data(i,1)
            min_T_location(i,1)=i; % 记录所有温度最小值所在位置
        else
            min_T_location(i,1)=0;
        end
    end

    min_P=data(1,2);
    for i=1:1:n_row
        if min_P>=data(i,2)
            min_P=data(i,2); % 压力最小值，MPa
        else
        end
    end

    for i=1:1:n_row
        if min_P==data(i,2)
            min_P_location(i,1)=i; % 记录所有压力最小值所在位置
        else
            min_P_location(i,1)=0;
        end
    end

    for i=1:1:n_row
        if min_T_location(i,1)==min_P_location(i,1) && min_T_location(i,1)~=0 && min_P_location(i,1)~=0
            min_TP_location=i; % 温度、压力最小值所在位置
        elseif min_T_location(i,1)~=0
            min_TP_location=i; % 温度最小值所在位置
        end
    end

    T_0=data(min_TP_location,1); % T_0为参考温度，℃
    P_0=data(min_TP_location,2); % P_0为参考压力，MPa

    %% 数据预处理（将原始数据处理成多元线性方程中的x1、x2、···、xp、y）（data→DATA，data一定是个3列的数据，而DATA则是个≥3列的数据）
    for i=1:1:n_row
        DATA(i,1)=data(i,2)-P_0;
        DATA(i,2)=(data(i,2)-P_0)^2;
        DATA(i,3)=data(i,1)-T_0;
        DATA(i,4)=(data(i,1)-T_0)*(data(i,2)-P_0);
        DATA(i,5)=(data(i,1)-T_0)*(data(i,2)-P_0)^2;
        DATA(i,6)=data(i,3);
    end

    %% 计算多元线性回归系数beta及相关系数r
    [beta,r] = Function_MultipleLinearRegression_WithIntercept(DATA);

    %% 将回归系数beta转化为模型中系数
    INTERCEPT=beta(1,1);
    ALPHA=beta(2,1);
    BETA=beta(3,1);
    GAMMA=beta(4,1);
    DELTA=beta(5,1);
    EPSILON=beta(6,1);

    %% 计算相对误差re及平均相对误差mre
    parameter_model=zeros(n_row,1); % 流体参数模型预测值初值

    for i=1:1:n_row
       parameter_model(i,1)=INTERCEPT+ALPHA*(data(i,2)-P_0)+BETA*(data(i,2)-P_0)^2+GAMMA*(data(i,1)-T_0)+DELTA*(data(i,1)-T_0)*(data(i,2)-P_0)+EPSILON*(data(i,1)-T_0)*(data(i,2)-P_0)^2; % 求流体参数模型预测值
    end

    re=zeros(n_row,1); % 相对误差初值
    for i=1:1:n_row
        re(i,1)=100*abs((data(i,3)-parameter_model(i,1))/data(i,3)); % 计算相对误差，%
    end

    mre=0; % 平均相对误差初值
    for i=1:1:n_row
        mre=mre+re(i,1)/n_row; % 计算平均相对误差，%
    end

    %% 计算标准化残差，查找并剔除离群值，构建新基础数据data_new
    % 计算标准化残差
    Q_e=0; % 残差平方和
    for i=1:1:n_row
        Q_e=Q_e+(data(i,3)-parameter_model(i,1))^2; % 计算残差平方和
    end

    ZRESID=zeros(n_row,1); % 标准化残差
    e=zeros(n_row,1); % 未标准化残差
    S_zre=(Q_e/(n_row-2))^0.5; % 残差的标准差
    for i=1:1:n_row
        e(i,1)=data(i,3)-parameter_model(i,1); % 未标准化残差=真实值-预测值
        ZRESID(i,1)=e(i,1)/S_zre; % 计算标准化残差
    end

    % 判断离群值所在位置
    ZRESID_OverTwoTimes=0;
    for i=1:1:n_row
        if abs(ZRESID(i,1))>2
            count(i,1)=i; % 是离群值时，记录所在位置
            ZRESID_OverTwoTimes=ZRESID_OverTwoTimes+1;
        else
            count(i,1)=0; % 不是离群值时，位置为0
        end
    end

    % 将离群值剔除，构建新基础数据
    co=1;
    for i=1:1:n_row
        if i==count(i,1)
        else
            data_new(co,:)=data(i,:); % 新基础数据
            co=co+1;
        end
    end

    % 确定数据几行几列
    n_new=size(data_new); % 数据（自变量+因变量）共多少行多少列
    n_row_new=n_new(1,1); % 样本观测次数
    n_col_new=n_new(1,2); % 自变量+因变量个数，因变量只有一个，所以自变量有（n_col_new-1）个

    if ZRESID_OverTwoTimes<=0 % 如果没有离群值
        %% 输出结果
%         fprintf('T_0(degree Celsius): %8.6f\n',T_0); % 参考温度，℃
%         fprintf('P_0(MPa): %8.6f\n',P_0); % 参考压力，MPa
%         fprintf('INTERCEPT: %8.6f\n',INTERCEPT); % 常数
%         fprintf('ALPHA: %8.6f\n',ALPHA); % 系数
%         fprintf('BETA: %8.6f\n',BETA); % 系数
%         fprintf('GAMMA: %8.6f\n',GAMMA); % 系数
%         fprintf('DELTA: %8.6f\n',DELTA); % 系数
%         fprintf('EPSILON: %8.6f\n',EPSILON); % 系数
%         fprintf('COD: %8.6f\n',r^2); % 决定系数
%         fprintf('MRE(Percent): %8.6f\n',mre); % 平均相对误差
        
        result(1,1)=T_0; % 参考温度，℃
        result(2,1)=P_0; % 参考压力，MPa
        result(3,1)=INTERCEPT; % 常数
        result(4,1)=ALPHA; % 系数
        result(5,1)=BETA; % 系数
        result(6,1)=GAMMA; % 系数
        result(7,1)=DELTA; % 系数
        result(8,1)=EPSILON; % 系数
        result(9,1)=r^2; % 决定系数
        result(10,1)=mre; % 平均相对误差
        
        %% 绘图
%         figure(1)
%         scatter3(data(:,1),data(:,2),data(:,3),'red','filled'); % 真实值散点图
%         hold on;
%         scatter3(data(:,1),data(:,2),parameter_model(:,1),'blue'); % 模型预测值散点图
%         xlabel('Temperature(degree Celsius)','FontName','黑体','FontSize',10);
%         ylabel('Pressure(MPa)','FontName','黑体','FontSize',10);
%         zlabel('parameter','FontName','黑体','FontSize',10);
%         legend('Actual Value','Predicted Value');

    elseif ZRESID_OverTwoTimes>0 % 如果有离群值
        %% 修正模型
        %% 以温度最小值且压力最小值作为模型基准参考值
        min_T_new=data_new(1,1);
        for i=1:1:n_row_new
            if min_T_new>=data_new(i,1)
                min_T_new=data_new(i,1); % 温度最小值，℃
            else
            end
        end

        for i=1:1:n_row_new
            if min_T_new==data_new(i,1)
                min_T_location_new(i,1)=i; % 记录所有温度最小值所在位置
            else
                min_T_location_new(i,1)=0;
            end
        end

        min_P_new=data_new(1,2);
        for i=1:1:n_row_new
            if min_P_new>=data_new(i,2)
                min_P_new=data_new(i,2); % 压力最小值，MPa
            else
            end
        end

        for i=1:1:n_row_new
            if min_P_new==data_new(i,2)
                min_P_location_new(i,1)=i; % 记录所有压力最小值所在位置
            else
                min_P_location_new(i,1)=0;
            end
        end

        for i=1:1:n_row_new
            if min_T_location_new(i,1)==min_P_location_new(i,1) && min_T_location_new(i,1)~=0 && min_P_location_new(i,1)~=0
                min_TP_location_new=i; % 温度、压力最小值所在位置
            elseif min_T_location_new(i,1)~=0
                min_TP_location_new=i; % 温度最小值所在位置
            end
        end

        T_0_new=data_new(min_TP_location_new,1); % T_0_new为参考温度，℃
        P_0_new=data_new(min_TP_location_new,2); % P_0_new为参考压力，MPa

        %% 数据预处理（将原始数据处理成多元线性方程中的x1、x2、···、xp、y）（data_new→DATA_new，data_new一定是个3列的数据，而DATA_new则是个≥3列的数据）
        for i=1:1:n_row_new
            DATA_new(i,1)=data_new(i,2)-P_0_new;
            DATA_new(i,2)=(data_new(i,2)-P_0_new)^2;
            DATA_new(i,3)=data_new(i,1)-T_0_new;
            DATA_new(i,4)=(data_new(i,1)-T_0_new)*(data_new(i,2)-P_0_new);
            DATA_new(i,5)=(data_new(i,1)-T_0_new)*(data_new(i,2)-P_0_new)^2;
            DATA_new(i,6)=data_new(i,3);
        end

        %% 计算多元线性回归系数beta_new及相关系数r_new
        [beta_new,r_new] = Function_MultipleLinearRegression_WithIntercept(DATA_new);

        %% 将回归系数beta_new转化为模型中系数
        INTERCEPT_new=beta_new(1,1);
        ALPHA_new=beta_new(2,1);
        BETA_new=beta_new(3,1);
        GAMMA_new=beta_new(4,1);
        DELTA_new=beta_new(5,1);
        EPSILON_new=beta_new(6,1);

        %% 计算相对误差re_new及平均相对误差mre_new
        parameter_model_new=zeros(n_row_new,1); % 流体参数模型预测值初值

        for i=1:1:n_row_new
            parameter_model_new(i,1)=INTERCEPT_new+ALPHA_new*(data_new(i,2)-P_0_new)+BETA_new*(data_new(i,2)-P_0_new)^2+GAMMA_new*(data_new(i,1)-T_0_new)+DELTA_new*(data_new(i,1)-T_0_new)*(data_new(i,2)-P_0_new)+EPSILON_new*(data_new(i,1)-T_0_new)*(data_new(i,2)-P_0_new)^2; % 求流体参数模型预测值
        end

        re_new=zeros(n_row_new,1); % 相对误差初值
        for i=1:1:n_row_new
            re_new(i,1)=100*abs((data_new(i,3)-parameter_model_new(i,1))/data_new(i,3)); % 计算相对误差，%
        end

        mre_new=0; % 平均相对误差初值
        for i=1:1:n_row_new
            mre_new=mre_new+re_new(i,1)/n_row_new; % 计算平均相对误差，%
        end

        %% 输出结果
%         fprintf('T_0(degree Celsius): %8.6f\n',T_0_new); % 参考温度，℃
%         fprintf('P_0(MPa): %8.6f\n',P_0_new); % 参考压力，MPa
%         fprintf('INTERCEPT: %8.6f\n',INTERCEPT_new); % 常数
%         fprintf('ALPHA: %8.6f\n',ALPHA_new); % 系数
%         fprintf('BETA: %8.6f\n',BETA_new); % 系数
%         fprintf('GAMMA: %8.6f\n',GAMMA_new); % 系数
%         fprintf('DELTA: %8.6f\n',DELTA_new); % 系数
%         fprintf('EPSILON: %8.6f\n',EPSILON_new); % 系数
%         fprintf('COD: %8.6f\n',r_new^2); % 决定系数
%         fprintf('MRE(Percent): %8.6f\n',mre_new); % 平均相对误差
        
        result(1,1)=T_0_new; % 参考温度，℃
        result(2,1)=P_0_new; % 参考压力，MPa
        result(3,1)=INTERCEPT_new; % 常数
        result(4,1)=ALPHA_new; % 系数
        result(5,1)=BETA_new; % 系数
        result(6,1)=GAMMA_new; % 系数
        result(7,1)=DELTA_new; % 系数
        result(8,1)=EPSILON_new; % 系数
        result(9,1)=r_new^2; % 决定系数
        result(10,1)=mre_new; % 平均相对误差

        %% 绘图
%         figure(2)
%         scatter3(data_new(:,1),data_new(:,2),data_new(:,3),'red','filled'); % 真实值散点图
%         hold on;
%         scatter3(data_new(:,1),data_new(:,2),parameter_model_new(:,1),'blue'); % 模型预测值散点图
%         xlabel('Temperature(degree Celsius)','FontName','黑体','FontSize',10);
%         ylabel('Pressure(MPa)','FontName','黑体','FontSize',10);
%         zlabel('parameter','FontName','黑体','FontSize',10);
%         legend('Actual Value','Predicted Value');

    end
end