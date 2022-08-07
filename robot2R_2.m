%建立模型
clc;clear all;
L=1;%连杆长
m=5;
I=[0,0,0;0,0,0;0,0,1/12*m*L^2];%惯性矩
g=9.8;

%D-H参数
L1=Revolute('d',0,'a',L,'alpha',0,'r',[-L/2,0,0],'I',I,'m',5,'qlim',[-2*pi,2*pi]);
L2=Revolute('d',0,'a',L,'alpha',0,'r',[-L/2,0,0],'I',I,'m',5,'qlim',[-2*pi,2*pi]);
robot = SerialLink([L1,L2],'name','2R','gravity',[0,g,0]);%x,y,z，y为竖直方向
robot.display();

%仿真选项（二维）
plotpot = {'workspace',[-5 5 -5 5 -5 5],'nobase','notiles','nowrist','top','fps',60};

%初始条件
theta=[5*pi/9,pi/2];
qz=theta-[pi/2,0];
qd=[-0.15,-0.25];%末端速度
y0=[qz]';
tspan=[0:0.04:10];%模拟时长10s
%tspan=[0:0.05:12];%模拟时长12s,i=234
tau=[0,0];

%解微分方程求轨迹
tic
f=@(t,y) [qd(1)*(cos(y(1)+y(2))/sin(y(2)))+qd(2)*(sin(y(1)+y(2))/sin(y(2)));
          qd(1)*(-(cos(y(1)+y(2))+cos(y(1)))/sin(y(2)))+qd(2)*(-(sin(y(1)+y(2))+sin(y(1)))/sin(y(2)))];
[tlist,ylist]=ode45(f,tspan,y0);
toc

%得到期待的角度相关值
q_d=ylist%角度对theta1、2
qd_d=[];%角速度
for i=1:251
    qd_d(i,:)=[qd(1)*(cos(q_d(i,1)+q_d(i,2))/sin(q_d(i,2)))+qd(2)*(sin(q_d(i,1)+q_d(i,2))/sin(q_d(i,2)));
               qd(1)*(-(cos(q_d(i,1)+q_d(i,2))+cos(q_d(i,1)))/sin(q_d(i,2)))+qd(2)*(-(sin(q_d(i,1)+q_d(i,2))+sin(q_d(i,1)))/sin(q_d(i,2)))]
end
qdd_d=[];%角加速度
qdd_d=diff(qd_d,1,1)/0.01;
qdd_d=[0,0;qdd_d];

q_d=q_d';
qd_d=qd_d';
qdd_d=qdd_d';

%运动轨迹
position_t=robot.fkine(q_d');
P=position_t.t;
px=[];pz=[];
for i=linspace(L,length(position_t),length(position_t))
    P=position_t(i).t;
    px=[px;P(1)];
    pz=[pz;P(2)];
end

%开始演示
figure(1);
robot.plot(ylist(:,L:robot.n),plotpot{:});hold on;%二维
plot(px,pz);xlabel("px/m");ylabel("pz/m");title("笛卡尔坐标空间轨迹图");

%关节力矩、角度、角速度、角加速度
%力矩
tau_d=[];M=[];D=[];C=[];G=[];
M11=(5/3)*m*L^2+m*L^2*cos(q_d(2,i));
M12=(1/3)*m*L^2+(1/2)*m*L^2*cos(q_d(2,i));%=M21
M22=(1/3)*m*L^2;
M0=[M11,M12;M12,M22];
D11=0;%D11=D22
D12=(-1/2)*m*L^2*sin(q_d(2,i));
D21=1/2*m*L^2*sin(q_d(2,i));
D0=[D11,D12;D21,D11];
for i=1:251
    M(i,:)=M0*qdd_d(:,i);
    D(i,:)=D0*qd_d(:,i);
    D(i,:)=D(i,:)*qd_d(:,i);
    C(i,:)=[-m*L^2*sin(q_d(2,i))*qd_d(1,i)*qd_d(2,i),0];
    G(i,:)=[(3/2)*m*g*L*sin(q_d(1,i))+(1/2)*m*g*L*sin(q_d(1,i)+q_d(2,i)),(1/2)*m*g*L*sin(q_d(1,i)+q_d(2,i))];
end
tau_d=M+D+C+G;

%作图
figure(2);
subplot(4,1,1)
plot(tspan,tau_d);title('力矩');legend('关节1','关节2');xlabel("s");ylabel("N.m");
subplot(4,1,2)
plot(tspan,q_d);title('角度');legend('关节1','关节2');xlabel("s");ylabel("rad");
subplot(4,1,3)
plot(tspan,qd_d);title('角速度');legend('关节1','关节2');xlabel("s");ylabel("rad/s");
subplot(4,1,4);
plot(tspan,qdd_d);title('角加速度');legend('关节1','关节2');xlabel("s");ylabel("rad/s^2");

%末端执行器笛卡尔坐标
figure(3);
plot(tspan,px);hold on;
plot(tspan,pz);
xlabel("s");ylabel("m");title("末端执行器笛卡尔坐标-时间轨迹图")
legend("x","z");

%雅可比矩阵计算线速度
xzd=[];
for i=1:251
    xzd(:,i)=[-(sin(q_d(1,i)+q_d(2,i))+sin(q_d(1,i))),-sin(q_d(1,i)+q_d(2,i));
        (cos(q_d(1,i)+q_d(2,i))+cos(q_d(1,i))),cos(q_d(1,i)+q_d(2,i))]*qd_d(:,i);
end
figure(4);
xd=xzd(1,:);
zd=xzd(2,:);
plot(tspan,xd);hold on;
plot(tspan,zd);
xlabel("s");ylabel("m/s");title("末端执行器线速度-时间轨迹图")
legend("x","z");

%奇异位型
% J=[-(sin(q_d(1,i)+q_d(2,i))+sin(q_d(1,i))),-sin(q_d(1,i)+q_d(2,i));
%     (cos(q_d(1,i)+q_d(2,i))+cos(q_d(1,i))),cos(q_d(1,i)+q_d(2,i))];=[J11,J12;J21,J22]

J=L^2*sin(q_d(2,:));%雅可比矩阵行列式值
w=[];%计算得w=sqrt((J11J22)^2+(J12J21)^2-2*J11J12J21J22)
for i=1:251
    w(:,i)=(-(sin(q_d(1,i)+q_d(2,i))+sin(q_d(1,i)))*-(sin(q_d(1,i)+q_d(2,i))+sin(q_d(1,i))))^2
           +(-sin(q_d(1,i)+q_d(2,i))*cos(q_d(1,i)+q_d(2,i))+cos(q_d(1,i)))^2
           -2*((-(sin(q_d(1,i)+q_d(2,i))+sin(q_d(1,i)))*-(sin(q_d(1,i)+q_d(2,i))+sin(q_d(1,i))))*(-sin(q_d(1,i)+q_d(2,i))*cos(q_d(1,i)+q_d(2,i))+cos(q_d(1,i))));
     w(:,i)=sqrt( w(:,i));
end
figure(5);
plot(tspan,w);hold on;
xlabel("t/s");ylabel("w");title("Yoshikawa");
res=find(abs(J)<0.1);%雅可比矩阵行列式值为0时奇异
scatter(tspan(res),w(res));
