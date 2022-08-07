%建立模型
clc;clear all;
L=1;%连杆长
m=5;
I=[0,0,0;0,0,0;0,0,1/12*m*L^2];%惯性矩
g=9.8;
%L1=Revolute('d',0,'a',L,'alpha',0,'r',[-L/2,0,0],'m',5,'qlim',[-2*pi,2*pi]);%'modified'
%L2=Revolute('d',0,'a',L,'alpha',0,'r',[-L/2,0,0],'m',5,'qlim',[-2*pi,2*pi]);没有转动惯量
%D-H参数
L1=Revolute('d',0,'a',L,'alpha',0,'r',[-L/2,0,0],'I',I,'m',5,'qlim',[-2*pi,2*pi]);%'modified'
L2=Revolute('d',0,'a',L,'alpha',0,'r',[-L/2,0,0],'I',I,'m',5,'qlim',[-2*pi,2*pi]);
robot = SerialLink([L1,L2],'name','2R','gravity',[0,g,0]);%x,y,z，y为竖直方向
robot.display();

%仿真选项（二维）
plotpot = {'workspace',[-5 5 -5 5 -5 5],'nobase','notiles','nowrist','top','fps',60};

%%第一问
%初始条件
theta=[pi/6,0];
qz=theta-[pi/2,0];
qd=[0,0];
y0=[qz,qd]';
tspan=[0:0.01:5];%5s
tau=[0,0];%力矩

%解微分方程
tic
[tlist,ylist]=ode45(@(t,y) [y(robot.n+1:end);robot.accel(y(1:robot.n)',y(robot.n+1:end)',tau)],tspan,y0);
toc
%获取笛卡尔坐标
position_t=robot.fkine(ylist(:,L:robot.n));
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
%robot.plot(ylist(:,L:robot.n));hold on;%三维
plot(px,pz);xlabel("px/m");ylabel("pz/m");title("笛卡尔坐标空间轨迹图");

%关节角度、角速度、角加速度
figure(2);
q_t=ylist(:,L:robot.n);
%subplot(1,3,1);
plot(tspan,q_t);title('角度');legend('关节1','关节2');xlabel("s");ylabel("rad");
figure(3);
qd_t=ylist(:,robot.n+1:end);
%subplot(1,3,2);
plot(tspan,qd_t);title('角速度');legend('关节1','关节2');xlabel("s");ylabel("rad/s");
figure(4);
qdd_t=robot.accel(q_t,qd_t,zeros(size(q_t)));
%subplot(1,3,3);
plot(tspan,qdd_t);title('角加速度');legend('关节1','关节2');xlabel("s");ylabel("rad/s^2");

%末端执行器的笛卡尔坐标
figure(5);
subplot(2,1,1);
plot(tspan,px);xlabel("s");ylabel("m");title("末端x坐标");hold on;
subplot(2,1,2);
plot(tspan,pz);xlabel("s");ylabel("m");title("末端z坐标");hold on;

%雅可比矩阵t=3.2s n=320,可查q_t表得
q_tn=[-1.8102 -0.1431]
J0=robot.jacob0(q_tn);

%roblocks
%mdl_robot
