function [K,n] = Rheology_PowerLaw(data_rheology)
    % ������������ģʽ�³��ϵ��K��Pa*s^n��������ָ��n
    % �������ݣ�data_rheology����T_0��P_0������ճ�ȼƲ������ݣ���һ��Ϊת�٣��ڶ���Ϊ����
    % ������ݣ�K�������ϵ����Pa*s^n����n��������ָ��
    
    N=data_rheology(:,1); % ����ճ�ȼ�ת�٣�rpm
    theta=data_rheology(:,2); % ����ճ�ȼƶ���
    m=length(N); % ����ճ�ȼ���������
    
    for i=1:1:m
        gammadot_reg(i)=log(1.7023*N(i)); % ʵ�ʼ������ʣ�s^-1��ȡ��Ȼ����
        tau_reg(i)=log(0.50771*theta(i)); % ʵ�ʼ���Ӧ����Pa��ȡ��Ȼ����
    end
    
    %% ʹ�����Իع������ϵ��K������ָ��n
    % ��ʽ�ο������, ������, �����. Ӧ������ͳ��[M]. �����Ƽ���ѧ������, 1995.
    x=gammadot_reg; % �Ա���x
    y=tau_reg; % �����y
    
    x_sum=0; % x���ֵ
    y_sum=0; % y���ֵ
    for i=1:m
        x_sum=x_sum+x(i); % ���Ա���x���
        y_sum=y_sum+y(i); % �������y���
    end
    
    x_ave=x_sum/m; % ���Ա���x��ƽ��ֵ
    y_ave=y_sum/m; % �������y��ƽ��ֵ
    
    L_xx=0; % (x-x_ave)^2���ֵ���������Ա������Ա���ƽ��ֵ���ֵƽ����
    L_yy=0; % (y-y_ave)^2���ֵ��������������������ƽ��ֵ���ֵƽ����
    L_xy=0; % (x-x_ave)*(y-y_ave)���ֵ
    for i=1:m
        L_xx=L_xx+(x(i)-x_ave)^2; % �����Ա������Ա���ƽ��ֵ���ֵƽ����
        L_yy=L_yy+(y(i)-y_ave)^2; % ����������������ƽ��ֵ���ֵƽ����
        L_xy=L_xy+(x(i)-x_ave)*(y(i)-y_ave); % (x-x_ave)*(y-y_ave)���ֵ
    end
    
    slope=L_xy/L_xx; % ����б��slope
    intercept=y_ave-(L_xy/L_xx)*x_ave; % ����ؾ�intercept
    r=L_xy/((L_xx*L_yy)^0.5); % �������ϵ��r
    R2=r^2; % ����ϵ��/����Ŷ�
    
    %% �����Իع����õ��Ľؾ��б�ʷֱ�ֵ�����ϵ��������ָ��
    K=exp(intercept); % ��e^�ؾ�ĸ�ֵ�����ϵ��
    n=slope; % ��б�ʸ�ֵ������ָ��
end

