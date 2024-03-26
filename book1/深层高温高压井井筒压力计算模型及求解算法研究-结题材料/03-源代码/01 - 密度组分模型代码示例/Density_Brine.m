function rho_b=Density_Brine(T,P,SaltContent)
    % 输入温度、压力、盐度，计算盐水的密度
    % 输入参数：压力（Pa）、温度（℃）、盐度（wt%）
    % 输出参数：密度（kg/m^3）
    
    % 清水密度计算
    B0=8.63186;B1=-3.31977E-3;B2=2.37170E-5;    % 1982-Sorelle模型中系数
    T=T*1.8+32;             % 温度从℃转化为℉
    P=P*0.00014504;         % 压力从Pa转化为psi
    rho_w=B0+B1*T+B2*P;     % 1982-Sorelle清水密度模型
    rho_w=rho_w*119.82643;  % 钻井液密度单位从lb/gal转化为kg/m^3
    if SaltContent==0 % 清水的密度
        rho_b=rho_w;
    else
        a1=10.0290;b1=1.68E-5;c1=1.11E-10;a2=-3.09E-3;b2=3.43E-8;c2=-6.36E-13;  % 2013-Zamora中19.3 wt% CaCl2盐水模型中的系数
        T=T*1.8+32;             % 温度从℃转化为℉
        P=P*0.00014504;         % 压力从Pa转化为psi
        rho_b1=(a1+b1*P+c1*P^2)+(a2+b2*P+c2*P^2)*T;  % 2013-Zamora中19.3 wt% CaCl2盐水模型
        rho_b1=rho_b1*119.82643;    % 钻井液密度单位从lb/gal转化为kg/m^3
        rho_b=rho_w+(rho_b1-rho_w)/(19.3-0)*(SaltContent-0);    % 对于任意输入盐度，将19.3 wt%和清水模型进行插值求解任意盐度下的密度值
    end
end