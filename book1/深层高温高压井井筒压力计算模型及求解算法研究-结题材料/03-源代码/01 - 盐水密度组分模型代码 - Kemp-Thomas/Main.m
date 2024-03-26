clear all;close all;clc;
%%
T=297; % 流体温度，K
P_0=0.1; % 大气压，MPa
P=0.1; % 流体压力，MPa
delta_P=P-P_0; % 压差，MPa

M_CaCl2=40+35.5*2; % 氯化钙的摩尔质量，g/mol
M_CaBr2=40+80*2; % 氯化钙的摩尔质量，g/mol

v_Ca_1=1; % 钙的化学计量数（总溶液）
v_Cl_1=1; % 氯的化学计量数（总溶液）
v_Br_1=1; % 溴的化学计量数（总溶液）
v_Ca_2=1; % 钙的化学计量数（总溶液中的氯化钙溶液部分）
v_Cl_2=2; % 氯的化学计量数（总溶液中的氯化钙溶液部分）
v_Ca_3=1; % 钙的化学计量数（总溶液中的溴化钙溶液部分）
v_Br_3=2; % 溴的化学计量数（总溶液中的溴化钙溶液部分）

z_Ca=2; % 钙离子电荷数
z_Cl=1; % 氯离子电荷数
z_Br=1; % 溴离子电荷数

rho_w=997.3; % 水的密度，kg/m^3

%%
F_CaCl2=0.24; % 氯化钙的质量分数
F_CaBr2=0.253; % 溴化钙的质量分数
F_H2O=1-F_CaCl2-F_CaBr2; % 水的质量分数

b_CaCl2=(1000*F_CaCl2)/(M_CaCl2*F_H2O); % 氯化钙的质量摩尔浓度，mol/kg
b_CaBr2=(1000*F_CaBr2)/(M_CaBr2*F_H2O); % 溴化钙的质量摩尔浓度，mol/kg

b_Ca=b_CaCl2+b_CaBr2; % 钙的质量摩尔浓度，mol/kg
b_Cl=2*b_CaCl2; % 氯的质量摩尔浓度，mol/kg
b_Br=2*b_CaBr2; % 溴的质量摩尔浓度，mol/kg

I=0.5*(b_Ca*2^2+b_Cl*1^2+b_Br*1^2); % 离子强度，mol/kg

Av_0=-0.7511;
Av_1=5.6658*10^(-6);
Av_2=1.5472*10^(-5);
Av_3=-4.3457*10^(-8);
Av_4=3.2265*10^(-11);
Av=1*10^(-6)*exp(Av_0+Av_1*T*delta_P+Av_2*T^2+Av_3*T^2*delta_P+Av_4*T^2*delta_P^2); % 水的德拜-休克尔极限定律斜率，m^3/mol

DH=((Av*I^0.5)/(2*(1+I^0.5)))*(v_Ca_1*z_Ca^2+v_Cl_1*z_Cl^2+v_Br_1*z_Br^2); % 摩尔体积计算公式中的德拜-休克尔项，m^3/mol

A_CaCl2=[-13.4695*10^(-5),0,48.2352*10^(-10);9.8240*10^(-7),0,-27.4851*10^(-12);-15.7830*10^(-10),1.2488*10^(-12),34.3893*10^(-15)];
PHI_V_0_CaCl2=0;
for i=1:1:3
    for j=1:1:3
        PHI_V_0_CaCl2=PHI_V_0_CaCl2+A_CaCl2(i,j)*T^(i-1)*delta_P^(j-1); % 氯化钙的无限稀释表观摩尔体积，m^3/mol
    end
end

A_CaBr2=[-13.6882*10^(-5),0,52.5124*10^(-10);10.5661*10^(-7),0,-28.5866*10^(-12);-16.5140*10^(-10),1.1752*10^(-12),37.3889*10^(-15)];
PHI_V_0_CaBr2=0;
for i=1:1:3
    for j=1:1:3
        PHI_V_0_CaBr2=PHI_V_0_CaBr2+A_CaBr2(i,j)*T^(i-1)*delta_P^(j-1); % 溴化钙的无限稀释表观摩尔体积，m^3/mol
    end
end

B_CaCl2=[4.9949*10^(-6),-7.1552*10^(-9),2.8028*10^(-11);-22.6474*10^(-9),0,0;3.3668*10^(-11),0,-2.2953*10^(-8)];
B_0_CaCl2=0;
for i=1:1:3
    for j=1:1:3
        if i==3 && j==3
            continue;
        else
            B_0_CaCl2=B_0_CaCl2+B_CaCl2(i,j)*T^(i-1)*delta_P^(j-1);
        end
    end
end

B_CaCl2_Ca=B_0_CaCl2+B_CaCl2(3,3)*(0.5*(v_Ca_2*z_Ca^2+v_Cl_2*z_Cl^2)/v_Ca_2)*b_Ca;
B_CaCl2_Cl=B_0_CaCl2+B_CaCl2(3,3)*(0.5*(v_Ca_2*z_Ca^2+v_Cl_2*z_Cl^2)/v_Cl_2)*b_Cl;

B_CaBr2=[5.8218*10^(-6),-7.4630*10^(-9),2.2843*10^(-11);-25.4804*10^(-9),0,0;3.6141*10^(-11),0,-3.4634*10^(-8)];
B_0_CaBr2=0;
for i=1:1:3
    for j=1:1:3
        if i==3 && j==3
            continue;
        else
            B_0_CaBr2=B_0_CaBr2+B_CaBr2(i,j)*T^(i-1)*delta_P^(j-1);
        end
    end
end

B_CaBr2_Ca=B_0_CaBr2+B_CaBr2(3,3)*(0.5*(v_Ca_3*z_Ca^2+v_Br_3*z_Br^2)/v_Ca_3)*b_Ca;
B_CaBr2_Br=B_0_CaBr2+B_CaBr2(3,3)*(0.5*(v_Ca_3*z_Ca^2+v_Br_3*z_Br^2)/v_Br_3)*b_Br;

N_Ca=(z_Ca^2*b_Ca)/(z_Ca^2*b_Ca);
N_Cl=(z_Cl^2*b_Cl)/(z_Cl^2*b_Cl+z_Br^2*b_Br);
N_Br=(z_Br^2*b_Br)/(z_Cl^2*b_Cl+z_Br^2*b_Br);

X3_CaCl2_n=0.5*(v_Cl_2*N_Cl*B_CaCl2_Cl+v_Br_3*N_Br*B_CaBr2_Br)*b_Ca;
X3_CaCl2_p=0.5*(v_Ca_2*N_Ca*B_CaCl2_Ca)*b_Cl;
X3_CaBr2_n=0.5*(v_Cl_2*N_Cl*B_CaCl2_Cl+v_Br_3*N_Br*B_CaBr2_Br)*b_Ca;
X3_CaBr2_p=0.5*(v_Ca_3*N_Ca*B_CaBr2_Ca)*b_Br;

X3_CaCl2=X3_CaCl2_n+X3_CaCl2_p;
X3_CaBr2=X3_CaBr2_n+X3_CaBr2_p;

PHI_V_CaCl2=PHI_V_0_CaCl2+DH+X3_CaCl2;
PHI_V_CaBr2=PHI_V_0_CaBr2+DH+X3_CaBr2;

rho=(1+b_CaCl2*M_CaCl2/1000+b_CaBr2*M_CaBr2/1000)/(1/rho_w+b_CaCl2*PHI_V_CaCl2+b_CaBr2*PHI_V_CaBr2); % 盐水密度，kg/m^3

