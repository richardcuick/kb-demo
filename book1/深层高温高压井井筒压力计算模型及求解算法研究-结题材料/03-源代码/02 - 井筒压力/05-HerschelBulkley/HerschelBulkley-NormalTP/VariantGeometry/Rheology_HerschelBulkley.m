function [tau_y,K,n]=Rheology_HerschelBulkley(data_rheology)
    % ���ݲ����¶ȡ�ѹ�������µĺհ���������ճ�ȼ����ݣ����ûع鷨����õ�����Ӧ��tau_y�����ϵ��K������ָ��n
    % ���㷽����ʹ����С���˷����лع����
    % ���������data_rheologyΪ�����¶ȡ�ѹ�������µĺհ���������ճ�ȼ����ݣ���һ��Ϊת�٣��ڶ���Ϊ����
    % ���������tau_yΪ����Ӧ����Pa����KΪ���ϵ����Pa��s^n����nΪ����ָ��

    N=data_rheology(:,1); % NΪ����ճ�ȼƵ�ת��
    theta=data_rheology(:,2); % thetaΪ����ճ�ȼƵĶ���

    m=length(N); % mΪ����ճ�ȼ����ݵ���������������ǵ��ֳ��޷���������ȫ������ʱ���꾮Һ���������Ȼ����ʹ�ûع鷽�����м���

    %% �������հ��������
    steps=floor(0.50771*theta(1)/0.001); % stepsΪ�������
    for j=1:steps
        tau_y_trial(j)=0.001*j; % ��ͬ��������µ�����Ӧ������ֵ
        for i=1:m % ���㲻ͬת���µ�ʵ�ʼ������ʺ�ʵ�ʼ���Ӧ������Ȼ����ֵ��logΪ��Ȼ����
            gamma_dot_hrreg_trial(i)=log(1.7023*N(i)); % gamma_dot_hrreg_trialΪʵ�ʼ������ʣ�s^-1������Ȼ����
            tau_hrreg_trial(i)=log(0.50771*theta(i)-tau_y_trial(j)); % tau_hrreg_trialΪʵ�ʼ���Ӧ����Pa������Ȼ����
        end

        % ʹ�����Իع��������Ӧ������ֵtau_y_trial��Ӧ�ĳ��ϵ��K_trial������ָ��n_trial
        % ��ʽ�ο������, ������, �����. Ӧ������ͳ��[M]. �����Ƽ���ѧ������, 1995.
        x=gamma_dot_hrreg_trial; % �Ա���x
        y=tau_hrreg_trial; % �����y

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
        
        % �����Իع����õ��Ľؾ��б�ʷֱ�ֵ�����ϵ��������ָ��
        K_trial(j)=exp(intercept); % ��e^�ؾ�ĸ�ֵ�����ϵ��
        n_trial(j)=slope; % ��б�ʸ�ֵ������ָ��

        % ����ƽ��������mre_h_trial
        for i=1:m
            gamma_dot_hrreg(i)=1.7023*N(i); % ����ʵ�ʼ�������
            tau_hrreg(i)=0.50771*theta(i); % ����ʵ�ʼ���Ӧ��
            tau_hrreg_model(i)=tau_y_trial(j)+K_trial(j)*(gamma_dot_hrreg(i))^n_trial(j); % ����K_trial��n_trial����õ�������Ӧ��
        end

        mre_h_trial=0; % ƽ��������
        for i=1:m
            re_h_trial=abs((tau_hrreg_model(i)-tau_hrreg(i))/tau_hrreg(i)); % ��i��������������ģ��Ԥ������Ӧ����ʵ������Ӧ��֮���������
            mre_h_trial=mre_h_trial+re_h_trial; % ���������
        end
        mre_h_trial=mre_h_trial/m; % ƽ��������=������֮��/����

        mre_h_trial_total(j)=mre_h_trial; % �洢ÿһ������ʱ��ƽ��������

        if n_trial(j)>1 % �������ָ������1������������ѭ����������һ��ѭ��
            mre_h_trial_total(j)=1; % ����ָ������1ʱ��ƽ�������ֵΪ1��ͬʱ��������ѭ��
            continue;
        end
    end

    %% ȷ��ƽ��������mre_h_trial_total����Сֵ�����Ӧ��j
    j_min=1; % �洢��������Сֵ��λ�ã������һ����ֵΪƽ����������Сֵ
    mre_h_trial_total_min=mre_h_trial_total(1); % �����һ����ֵΪƽ����������Сֵ
    for i=2:steps % ������������Сֵ
        if mre_h_trial_total(i)<mre_h_trial_total_min
            mre_h_trial_total_min=mre_h_trial_total(i); % ���ĳ��������С�ڴ洢����������Сֵ�����滻��������Сֵ
            j_min=i; % �洢��������Сֵ��λ��
        end
    end

    %% �����Сƽ��������ʱ���������Ϊ���յĺհ��������
    tau_y=tau_y_trial(j_min);
    K=K_trial(j_min);
    n=n_trial(j_min);
end