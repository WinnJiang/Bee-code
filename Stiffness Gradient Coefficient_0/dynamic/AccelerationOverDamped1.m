clear;clc
%过阻尼
%导入位移数据
load displacement1.txt

N=1;
%zetaspan=0.1:0.1:1;
displacement=displacement1;
timespan=displacement(1:3001,1);
startpoint=1;
for i=1:1
    %分配各阻尼下端点位移
    endpoint=3000+startpoint;
    dis(:,i)=displacement(startpoint:endpoint,2);
    startpoint=endpoint+1;
       
    %定义一个符号变量t
    syms t
    %编辑要拟合的公式，设置变量，设置系数
    ft=fittype('exp(-zeta*w*t)*(a*exp((zeta^2-1)^0.5*w*t)+b*exp(-(zeta^2-1)^0.5*w*t))','independent','t','coefficients',{'a','b','zeta','w'});
    %ft=fittype('exp(-zeta*w*t)*(a*exp((zeta^2-1)^0.5*w*t))','independent','t','coefficients',{'a','zeta','w'});
    options = fitoptions(ft);
    %设定拟合初值
    a0(i)=2.048;%初始振幅初值
    b0(i)=-0.3868;
    %b0(i)=max(dis(:,i));%初始振幅初值
    zeta0(i)=i*0.1+1;%阻尼系数初值
    %T(i)=2*abs(timespan(find(dis(:,i)==max(dis(:,i))))-timespan(find(dis(:,i)==min(dis(:,i)))));%周期
    w0(i)=3.3e5;%固有频率初值
    %phi(i)=pi/2;%相位初值
    options.StartPoint=[a0(i) b0(i) zeta0(i) w0(i)];
    %设定拟合所用系数的上下限
    options.Lower = [0 -100 0 1e3];
    options.Upper=[100 0 10 1e8];
    %进行拟合
    cfun=fit(timespan(N:end),dis(N:end,i),ft,options);
    para(:,i)=[cfun.a,cfun.b,cfun.zeta,cfun.w];
    %显示拟合函数，数据必须为列向量形式
    yi(:,i)=cfun(timespan(N:end));
    figure (1)
    hold on
    plot(timespan(N:end),dis(N:end,i),'r',timespan(N:end),yi(:,i),'b-');%单位为微米
    xlabel('Time (s)');
    ylabel('Displacement (um)');
    title('位移');
    if i==10
        hold off
    else
        hold on
    end
    %两次求导求加速度函数
    syms x
    disfun(x)=1e-6*exp(-cfun.zeta*cfun.w*x)*(cfun.a*exp((cfun.zeta^2-1)^0.5*cfun.w*x)+cfun.b*exp(-(cfun.zeta^2-1)^0.5*cfun.w*x));%位移函数，单位化为m
    %disfun(x)=1e-6*exp(-cfun.zeta*cfun.w*x)*(cfun.a*exp((cfun.zeta^2-1)^0.5*cfun.w*x));
    velfun(x)=diff(disfun(x));%速度函数
    accfun(x)=diff(velfun(x));%加速度函数
    vel(:,i)=eval(velfun(timespan(N:end)));%计算速度
    acc(:,i)=eval(accfun(timespan(N:end)));%计算加速度
    maxacc(i)=max(acc(:,i));
    %绘图加速度
    figure(2)
    hold on
    plot(timespan(N:end), acc(:,i));
    xlabel('Time (s)');
    ylabel('Accleration (m/s2)');
    title('加速度');
    if i==20
        hold off
    else
        hold on
    end
    figure(3)
    hold on
    plot(timespan(N:end), vel(:,i));
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    title('速度');
    if i==20
        hold off
    else
        hold on
    end
end

% figure (1)
% hold on
% for i=1:10
%     plot(t(1:end-1),acc(:,i))
% end
% hold off
% figure(2)
% hold on
% for i=1:10
%     plot(t(1:end),vel(:,i))
% end
% hold off
% figure(2)
% hold on
% for i=1:10
%     plot(t(1:end),vel(:,i))
% end
% hold off


%% %curve fitting using costum function

%%
