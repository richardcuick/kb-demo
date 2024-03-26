function [beta,r] = Function_MultipleLinearRegression_WithIntercept(DATA)
    % 输入：DATA为n×(p+1)矩阵，每行数据代表一组观测值，共n行；前p列存储自变量数据、最后一列存储因变量数据，共p+1列
    % 输出：beta为(p+1)×1矩阵，第一行数据是常数，其余为自变量前系数
    % 输出：r为相关系数（r^2为决定系数/拟合优度）
    
    %% 计算多元线性回归系数
    n=size(DATA); % 判断数据（自变量+因变量）共多少行多少列
    n_row=n(1,1); % 样本观测次数
    n_col=n(1,2); % 自变量+因变量个数，因变量只有一个，所以自变量有（n_col-1）个
    
    % 求变量平均值variable_ave（最后一个为因变量平均值，其余为自变量平均值）
    variable_ave=zeros(n_col,1); % 变量平均值初值设为0
    for i=1:1:n_col
        for j=1:1:n_row
            variable_ave(i,1)=variable_ave(i,1)+DATA(j,i)/n_row; % 求变量平均值
        end
    end
    
    % 计算L_ij（p×p方阵，第一个字母代表行，第二个代表列）
    L_ij=zeros(n_col-1); % L_ij初值设为0
    for i=1:1:n_col-1
        for j=1:1:n_col-1
            for t=1:1:n_row
                L_ij(i,j)=L_ij(i,j)+(DATA(t,i)-variable_ave(i,1))*(DATA(t,j)-variable_ave(j,1)); % 计算L_ij
            end
        end
    end
    
    % 计算L_iy（p×1矩阵，第一个字母代表行，第二个代表列）
    L_iy=zeros(n_col-1,1); % L_iy初值设为0
    for i=1:1:n_col-1
        for t=1:1:n_row
            L_iy(i,1)=L_iy(i,1)+(DATA(t,i)-variable_ave(i,1))*(DATA(t,n_col)-variable_ave(n_col,1)); % 计算L_iy
        end
    end
    
    % 计算自变量前系数
%     beta_x=pinv(L_ij)*L_iy; % 自变量前系数（涉及到俩矩阵运算①矩阵求逆；②矩阵乘法）
    beta_x=Function_MatrixMultiplication(Function_InverseMatrix(L_ij),L_iy); % 自变量前系数（涉及到俩矩阵运算①矩阵求逆；②矩阵乘法）
    
    % 计算常数
    beta_0=variable_ave(n_col,1); % 常数初值
    for i=1:1:n_col-1
        beta_0=beta_0-variable_ave(i,1)*beta_x(i,1); % 计算常数
    end
    
    beta=[beta_0;beta_x]; % 多元线性回归系数
    
    %% 计算相关系数r
    L_yy=0; % 总平方和初值
    U=0; % 回归平方和初值
    
    for i=1:1:n_row
        L_yy=L_yy+(DATA(i,n_col)-variable_ave(n_col,1))^2; % 计算总平方和
    end
    
    for i=1:1:n_col-1
        U=U+beta_x(i,1)*L_iy(i,1); % 计算回归平方和
    end
    
    r=(U/L_yy)^0.5; % 相关系数
    
end

