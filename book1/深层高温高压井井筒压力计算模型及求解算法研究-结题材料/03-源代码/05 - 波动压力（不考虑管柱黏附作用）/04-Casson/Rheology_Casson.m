function [tau_y_c,eta_infinity]=Rheology_Casson(data_rheology)
    % ���ݲ����¶ȡ�ѹ�������µĿ�ɭ��������ճ�ȼ����ݣ����ûع鷨����õ���ɭ����Ӧ��tau_y_c�Ϳ�ɭ��eta_infinity
    % ���㷽����ʹ����С���˷����лع����
    % ���������data_rheologyΪ�����¶ȡ�ѹ�������µĿ�ɭ��������ճ�ȼ����ݣ���һ��Ϊת�٣��ڶ���Ϊ����
    % ���������tau_y_cΪ��ɭ����Ӧ����Pa����eta_infinityΪ��ɭ�ȣ�Pa��s��
    
    N=data_rheology(:,1); % NΪ����ճ�ȼƵ�ת��
    theta=data_rheology(:,2); % thetaΪ����ճ�ȼƵĶ���
    
    m=length(N); % mΪ����ճ�ȼ����ݵ���������������ǵ��ֳ��޷���������ȫ������ʱ���꾮Һ���������Ȼ����ʹ�ûع鷽�����м���
    
    for i=1:m % ���㲻ͬת���µ�ʵ�ʼ������ʺ�ʵ�ʼ���Ӧ��
        gamma_dot_reg(i)=(1.7023*N(i))^0.5; % gamma_dot_regΪʵ�ʼ������ʣ�s^-1��
        tau_reg(i)=(0.50771*theta(i))^0.5; % tau_regΪʵ�ʼ���Ӧ����Pa��
    end
    
    %% ʹ�����Իع���㿨ɭ����Ӧ��tau_y_c�Ϳ�ɭ��eta_infinity
    % ��ʽ�ο������, ������, �����. Ӧ������ͳ��[M]. �����Ƽ���ѧ������, 1995.
    x=gamma_dot_reg; % �Ա���x
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
    
    %% �����Իع����õ��Ľؾ��б�ʷֱ�ֵ����ɭ����Ӧ���Ϳ�ɭճ��
    tau_y_c=intercept^2; % ���ؾำֵ����ɭ����Ӧ��
    eta_infinity=slope^2; % ��б�ʸ�ֵ����ɭճ��
end