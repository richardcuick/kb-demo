function [beta,r] = Function_MultipleLinearRegression_WithoutIntercept(DATA)
    %% 说明：本函数为常数项为0的多元线性回归函数
    % 输入数据：DATA为n×(p+1)矩阵，每行数据代表一组观测值，共n行；前p列存储自变量数据、最后一列存储因变量数据，共p+1列
    % 输出数据：beta为自变量前系数，p×1矩阵
    % 输出数据：r为相关系数（r^2为决定系数/拟合优度）
    
    %% 判断样本数、自变量及因变量个数
    n=size(DATA); % 判断数据（自变量+因变量）共多少行多少列
    n_row=n(1,1); % 样本观测次数
    n_col=n(1,2); % 自变量+因变量个数，因变量只有一个，所以自变量有（n_col-1）个
    
    %% 计算多元线性回归系数    
    % 计算L_ij（p×p方阵，i代表行，j代表列）
    L_ij=zeros(n_col-1); % L_ij初值设为0
    for i=1:1:n_col-1
        for j=1:1:n_col-1
            for t=1:1:n_row
                L_ij(i,j)=L_ij(i,j)+DATA(t,i)*DATA(t,j); % 计算L_ij
            end
        end
    end
    
    % 计算L_iy（p×1矩阵，i代表行）
    L_iy=zeros(n_col-1,1); % L_iy初值设为0
    for i=1:1:n_col-1
        for t=1:1:n_row
            L_iy(i,1)=L_iy(i,1)+DATA(t,i)*DATA(t,n_col); % 计算L_iy
        end
    end
    
    % 计算自变量前系数
%     beta=pinv(L_ij)*L_iy; % 自变量前系数（涉及到俩矩阵运算①矩阵求逆；②矩阵乘法）
    beta=Function_MatrixMultiplication(Function_InverseMatrix(L_ij),L_iy); % 自变量前系数（涉及到俩矩阵运算①矩阵求逆；②矩阵乘法）
    
    %% 计算相关系数r
    U=0; % 回归平方和初值
    Q_e=0; % 残差平方和初值  
    
    for i=1:1:n_col-1
        U=U+beta(i,1)*L_iy(i,1); % 计算回归平方和
    end
    
    model=zeros(n_row,1);
    for i=1:1:n_row
        for j=1:1:n_col-1
            model(i,1)=model(i,1)+beta(j,1)*DATA(i,j); % 计算模型预测值
        end
    end
    
    for i=1:1:n_row
        Q_e=Q_e+(DATA(i,n_col)-model(i,1))^2; % 计算残差平方和
    end
    
    r=(U/(Q_e+U))^0.5; % 相关系数
    
end