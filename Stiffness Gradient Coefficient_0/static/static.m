clear all
clc

%参数设置
L=38.54e-6; %初始长度 米
alpha=pi/2;      %集中力作用角
d = 6.21e-6;%直径 m
A = pi*d^2/4;%面积
I = pi*d^4/64;%惯性矩

E_base = 22.8e6;%杨氏模量 Pa
E_gradient=0;%/微米

h=1e-9;%仿真计算步长 m
Lspan=38.54e-6:-1e-8:38.54e-6;%长度范围 m 
Pspan=0.3193e-5:1e-7:0.3193e-5;%力的范围 N
tspan=0:h:L;%仿真时间范围

for iii=1:length(Lspan)
    L=Lspan(iii);
    t_length(iii)=length(0:h:L); %由于长度改变，仿真总时间范围也会改变
    for ii=1:length(Pspan)      
        z(1)=alpha;z(2)=0;%初始条件
        %P=3.491e-5;
        P=Pspan(ii);%集中力 N
       
        Result=zeros(t_length(iii),2);%定义结果矩阵
        % 进行当前工况的结果
        for i=1:t_length(iii)   
           s=tspan(i);
           z=RK(z,s,h,E_base,E_gradient,L,I,P);
           Result(i,:)=z;   % 将每次迭代的值保存
        end
        %定义转角phi矩阵
        phi=zeros(1,t_length(iii));
        %计算当前工况的转角phi
        phi=Result(1:t_length(iii),1)-Result(t_length(iii),1).*ones(t_length(iii),1);
        %求当前工况下梁的各节点挠度,小到大指固定端到自由端
        N=100;%定义100个节点
        for j=1:N
            %p_start=length(phi(iii,ii,:))/N*(N-j)+1;
            p_start=fix(t_length(iii)/N*(N-j))+1;%向0取整
            S(j,1)=sum(cos(phi(p_start:t_length(iii),1))*h);%X长度
            S(j,2)=sum(sin(phi(p_start:t_length(iii),1))*h);%Y挠度
        end
        %寻找当前工况下X长度最大值，计算偏差，并保存
        X00(iii,ii)=max(S(:,1));
        S=zeros(N,2);
        diff(iii,ii)=abs(X00(iii,ii)-22.54e-6);
        
        fprintf('L=%d, P=%s\n', L, P);%监控计算进度
    end
end
min_value=min(min(diff))%寻找最小值
[L_coordinate,P_coordinate] =find(diff==min_value);%寻找最小值对应的坐标
P_value=Pspan(P_coordinate);%力的大小
P_position=Lspan(L_coordinate);%力的作用位置
X0_value=X00(L_coordinate,P_coordinate);%对应的X长度

%计算寻找出的最优值的数据

%初始条件
z(1)=alpha;z(2)=0;
tspan_O=0:h:P_position;
Result_O=zeros(length(tspan_O),2);%定义结果矩阵
for i=1:length(tspan_O)   % 进行多次迭代
    s_O=tspan_O(i);
    z=RK(z,s_O,h,E_base,E_gradient,P_position,I,P_value);
    Result_O(i,:)=z;   % 将每次迭代的值保存
end

%计算转角phi
phi_O=Result_O(:,1)-Result_O(end,1).*ones(length(tspan_O),1);
%求各节点挠度,小到大指固定端到自由端
N=101;%101个节点
for j=1:N
    p_start_O=fix(length(phi_O)/N*(N-j))+1;
    S_O(j,1)=sum(cos(phi_O(p_start_O:end,1))*h);%X长度
    S_O(j,2)=sum(sin(phi_O(p_start_O:end,1))*h);%Y挠度
end
figure(1)
plot(tspan_O,phi_O)
figure(2)
plot(S_O(:,1),S_O(:,2))

%导出初始变形条件
%定义11个节点共66个变量，x,y,z, dx/ds,dy/ds,dz/ds
Node_IniDefor=zeros(66,1);
Node_IniDefor(1:6:61)=S_O(1:10:101,1);%x
Node_IniDefor(2:6:62)=S_O(1:10:101,2);%y
Node_IniDefor(64:-6:4)=cos(phi_O(1:(length(phi_O)-1)/10:length(phi_O)));%dx/ds
Node_IniDefor(65:-6:5)=sin(phi_O(1:(length(phi_O)-1)/10:length(phi_O)));%dy/ds
%将初始条件保存进Excel
writematrix(Node_IniDefor,'Node_IniDefor.xls','sheet',1,'range','A1:A66');

function z3=Dif(E_base,E_gradient,L,s,I,z2,z1,P)
z3=-(-E_base*E_gradient*exp(E_gradient*(L-s))*I*z2+P*sin(z1))/(E_base*exp(E_gradient*(L-s))*I);
end

function z=RK(z,s,h,E_base,E_gradient,L,I,P)
KK1=z(2);
KK2=z(2)+h/2*KK1;
KK3=z(2)+h/2*KK2;
KK4=z(2)+h*KK3;
z(1)=z(1)+h/6*(KK1+2*KK2+2*KK3+KK4);
K1=Dif(E_base,E_gradient,L,s,I,z(2),z(1),P);
K2=Dif(E_base,E_gradient,L,s,I,z(2)+h/2*K1,z(1)+h/2,P);
K3=Dif(E_base,E_gradient,L,s,I,z(2)+h/2*K2,z(1)+h/2,P);
K4=Dif(E_base,E_gradient,L,s,I,z(2)+h*K3,z(1)+h,P);
z(2)=z(2)+h/6*(K1+2*K2+2*K3+K4);
end

