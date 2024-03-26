function rho_o=Density_Oil(T,P,OilType)
    % 输入温度、压力、油相类型，计算油相的密度
    % 输入参数：压力（Pa）、温度（℃）、油相类型
    % 输出参数：密度（kg/m^3）
    
    T=T*1.8+32;             % 温度从℃转化为℉
    P=P*0.00014504;         % 压力从Pa转化为psi
    if OilType==1   % 柴油
        % 2013-Zamora柴油D1模型
        a1=7.1101;b1=2.61E-5;c1=-1.23E-10;a2=-2.90E-3;b2=8.14E-8;c2=-1.37E-12;
    elseif OilType==2   % 矿物油
        % 2013-Zamora矿物油MO1模型
        a1=6.7422;b1=3.32E-5;c1=-3.46E-10;a2=-2.85E-3;b2=6.31E-8;c2=-8.40E-13;
    elseif OilType==3   % 合成有机物
        % 2013-Zamora合成有机物S1模型
        a1=6.6962;b1=2.83E-5;c1=-1.90E-10;a2=-2.72E-3;b2=6.87E-8;c2=-1.00E-12;
    end
    rho_o=(a1+b1*P+c1*P^2)+(a2+b2*P+c2*P^2)*T;  % 2013-Zamora油相密度模型，密度单位为lb/gal
    rho_o=rho_o*119.82643;  % 钻井液密度单位从lb/gal转化为kg/m^3
end