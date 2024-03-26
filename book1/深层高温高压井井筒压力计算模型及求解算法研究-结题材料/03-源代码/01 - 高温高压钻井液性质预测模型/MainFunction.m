close all;clear all;clc; % 关闭图窗、清空变量、清除命令窗口内容
%% 导入基础数据data
data=csvread('data18.csv',1,0); % data存储不同温度、压力及相应流体参数，其中第一列为温度（℃），第二列为压力（MPa），其余为相应流体参数
n=size(data); % 数据（自变量+因变量）共多少行多少列
n_row=n(1,1); % 样本观测次数
n_col=n(1,2); % 自变量+因变量个数（因变量n_col-2个）

%% 模型参数计算
for i=1:1:n_col-2 % i代表列，也代表因变量
    Data=[data(:,1),data(:,2),data(:,i+2)]; % 温度、压力与第i个因变量构成的n_row×3矩阵
    
    % 第1个模型（Sorelle – 1982）
    [result_01] = Function_ParameterSol01(Data); % 第1个模型参数计算结果
    rc_01=size(result_01);
    row_01=rc_01(1,1); % 第1个模型参数计算结果构成的列向量有几行
    for j=1:1:row_01
        RESULT(j,i)=result_01(j,1); % 将第1个模型参数计算结果存储到RESULT中
    end
    
    % 第2个模型（Politte – 1985）
    [result_02] = Function_ParameterSol02(Data); % 第2个模型参数计算结果
    rc_02=size(result_02);
    row_02=rc_02(1,1); % 第2个模型参数计算结果构成的列向量有几行
    for j=1:1:row_02
        RESULT(j+row_01+1,i)=result_02(j,1); % 将第2个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01,i)=result_02(j,1); % 将第2个模型参数计算结果存储到RESULT中
    end
    
    % 第3个模型（Minton – 1988）
    [result_03] = Function_ParameterSol03(Data); % 第3个模型参数计算结果
    rc_03=size(result_03);
    row_03=rc_03(1,1); % 第3个模型参数计算结果构成的列向量有几行
    for j=1:1:row_03
        RESULT(j+row_01+row_02+2,i)=result_03(j,1); % 将第3个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01+row_02,i)=result_03(j,1); % 将第3个模型参数计算结果存储到RESULT中
    end
    
    % 第4个模型（Fisk – 1989）
    [result_04] = Function_ParameterSol04(Data); % 第4个模型参数计算结果
    rc_04=size(result_04);
    row_04=rc_04(1,1); % 第4个模型参数计算结果构成的列向量有几行
    for j=1:1:row_04
        RESULT(j+row_01+row_02+row_03+3,i)=result_04(j,1); % 将第4个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01+row_02+row_03,i)=result_04(j,1); % 将第4个模型参数计算结果存储到RESULT中
    end
    
    % 第5个模型（鄢捷年 – 1990）
    [result_05] = Function_ParameterSol05(Data); % 第5个模型参数计算结果
    rc_05=size(result_05);
    row_05=rc_05(1,1); % 第5个模型参数计算结果构成的列向量有几行
    for j=1:1:row_05
        RESULT(j+row_01+row_02+row_03+row_04+4,i)=result_05(j,1); % 将第5个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01+row_02+row_03+row_04,i)=result_05(j,1); % 将第5个模型参数计算结果存储到RESULT中
    end
    
    % 第6个模型（Eirik – 1998）
    [result_06] = Function_ParameterSol06(Data); % 第6个模型参数计算结果
    rc_06=size(result_06);
    row_06=rc_06(1,1); % 第6个模型参数计算结果构成的列向量有几行
    for j=1:1:row_06
        RESULT(j+row_01+row_02+row_03+row_04+row_05+5,i)=result_06(j,1); % 将第6个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01+row_02+row_03+row_04+row_05,i)=result_06(j,1); % 将第6个模型参数计算结果存储到RESULT中
    end
    
    % 第7个模型（周福建 – 1999）
    [result_07] = Function_ParameterSol07(Data); % 第7个模型参数计算结果
    rc_07=size(result_07);
    row_07=rc_07(1,1); % 第7个模型参数计算结果构成的列向量有几行
    for j=1:1:row_07
        RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+6,i)=result_07(j,1); % 将第7个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06,i)=result_07(j,1); % 将第7个模型参数计算结果存储到RESULT中
    end
    
    % 第8个模型（汪海阁 – 2000）
    [result_08] = Function_ParameterSol08(Data); % 第8个模型参数计算结果
    rc_08=size(result_08);
    row_08=rc_08(1,1); % 第8个模型参数计算结果构成的列向量有几行
    for j=1:1:row_08
        RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07+7,i)=result_08(j,1); % 将第8个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07,i)=result_08(j,1); % 将第8个模型参数计算结果存储到RESULT中
    end
    
    % 第9个模型（张琰 – 2000）
    [result_09] = Function_ParameterSol09(Data); % 第9个模型参数计算结果
    rc_09=size(result_09);
    row_09=rc_09(1,1); % 第9个模型参数计算结果构成的列向量有几行
    for j=1:1:row_09
        RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07+row_08+8,i)=result_09(j,1); % 将第9个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07+row_08,i)=result_09(j,1); % 将第9个模型参数计算结果存储到RESULT中
    end
    
    % 第10个模型（赵胜英 – 2009）
    [result_10] = Function_ParameterSol10(Data); % 第10个模型参数计算结果
    rc_10=size(result_10);
    row_10=rc_10(1,1); % 第10个模型参数计算结果构成的列向量有几行
    for j=1:1:row_10
        RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07+row_08+row_09+9,i)=result_10(j,1); % 将第10个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07+row_08+row_09,i)=result_10(j,1); % 将第10个模型参数计算结果存储到RESULT中
    end
    
    % 第11个模型（Zamora – 2013）
    [result_11] = Function_ParameterSol11(Data); % 第11个模型参数计算结果
    rc_11=size(result_11);
    row_11=rc_11(1,1); % 第11个模型参数计算结果构成的列向量有几行
    for j=1:1:row_11
        RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07+row_08+row_09+row_10+10,i)=result_11(j,1); % 将第11个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07+row_08+row_09+row_10,i)=result_11(j,1); % 将第11个模型参数计算结果存储到RESULT中
    end
    
    % 第12个模型（赵向阳 – 2013）
    [result_12] = Function_ParameterSol12(Data); % 第12个模型参数计算结果
    rc_12=size(result_12);
    row_12=rc_12(1,1); % 第12个模型参数计算结果构成的列向量有几行
    for j=1:1:row_12
        RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07+row_08+row_09+row_10+row_11+11,i)=result_12(j,1); % 将第12个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07+row_08+row_09+row_10+row_11,i)=result_12(j,1); % 将第12个模型参数计算结果存储到RESULT中
    end
    
    % 第13个模型（高禹 – 2019）
    [result_13] = Function_ParameterSol13(Data); % 第13个模型参数计算结果
    rc_13=size(result_13);
    row_13=rc_13(1,1); % 第13个模型参数计算结果构成的列向量有几行
    for j=1:1:row_13
        RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07+row_08+row_09+row_10+row_11+row_12+12,i)=result_13(j,1); % 将第13个模型参数计算结果存储到RESULT中
%         RESULT(j+row_01+row_02+row_03+row_04+row_05+row_06+row_07+row_08+row_09+row_10+row_11+row_12,i)=result_13(j,1); % 将第13个模型参数计算结果存储到RESULT中
    end
    
end

RESULT
