function matrix3 = Function_MatrixMultiplication(matrix1,matrix2)
    %% 说明：本函数为矩阵乘法函数
    % 输入数据：Matrix1为m×s矩阵，Matrix2为s×n矩阵
    % 输出数据：Matrix3为m×n矩阵
    
    [m,s1]=size(matrix1);
    [s2,n]=size(matrix2);
    matrix3=zeros(m,n);
    if(s1~=s2)
        fprintf('Please Check The Input Data!');
    else
        for i=1:m % Matrix1的行
            for j=1:n % Matrix2的列
                for k=1:s1
                    matrix3(i,j)=matrix3(i,j)+matrix1(i,k)*matrix2(k,j);
                end
            end
        end
    end
end