function [InverseMatrix] = Function_InverseMatrix(Matrix)
    %% 说明：本函数利用待定系数法及列主元三角分解法求逆函数
    % 输入数据：Matrix为n×n方阵
    % 输出数据：InverseMatrix为Matrix的逆矩阵，n×n方阵
    
    n=length(Matrix); % 方阵Matrix的维度
    b=eye(n); % n阶单位阵b
    
    % 求LU
    for i=1:1:n
        p(i,1)=i; % 排列阵P，用向量p表示
    end
    LU=Matrix; 
    
    for k=1:1:n-1
        [s,i]=max(abs(LU(k:n,k))); % 查找方阵Matrix的下三角元素的列主元s及其所在行号i（i对应的行号为下三角列元素最大值的位置行号！！）
        ik=i+k-1; % ik为列主元s在方阵Matrix中的行号

        if s==0 % 存在列主元最大值为0，此时矩阵Matrix为奇异矩阵，不可逆，没有逆矩阵
            quit;
        end

        if ik~=k % 如果列主元s的位置ik不在对角线上，此时进行行交换以使得列主元s在对角线位置k处
            m=p(k,1); % 中间变量m
            p(k,1)=p(ik,1); % 将“第k列元素最大值所在行号ik”赋值给排列阵p的第k个元素，进而用排列阵第k个元素的数值表征第k列元素最大值所在行号
            p(ik,1)=m; % 将第k列对角线处元素行号k赋值给p(ik,1)
            
            lk=LU(k,:); % 中间变量lk
            LU(k,:)=LU(ik,:); % 将第k列元素最大值所在的第ik行的元素交换到第k行
            LU(ik,:)=lk; % 将第k行元素交换到第ik行
        end

        LU(k+1:n,k)=LU(k+1:n,k)/LU(k,k);
        LU(k+1:n,k+1:n)=LU(k+1:n,k+1:n)-LU(k+1:n,k)*LU(k,k+1:n);
    end
    
    % 求逆矩阵InverseMatrix
    y=zeros(n,1);
    for ii=1:1:n
        y(1,1)=b(p(1,1),ii);
        for i=2:1:n
            y(i,1)=b(p(i,1),ii)-LU(i,1:i-1)*y(1:i-1,1); % 向前回代求y
        end

        InverseMatrix(n,ii)=y(n,1)/LU(n,n);
        for i=(n-1):-1:1
            InverseMatrix(i,ii)=(y(i,1)-LU(i,i+1:n)*InverseMatrix(i+1:n,ii))/LU(i,i); % 向后回代求InverseMatrix
        end
    end
end