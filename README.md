# ROBOTICS-2R
机器人学大作业- matlab 仿真二连杆机器人的控制


1）设机器人关节空间初始状态为 θ=[30°,0°]，dθ=[0,0] ，关节力矩为[0,0]，利用 Simulink 进行时长 5s 的运动学仿真，并解决如下问题：

画出机器人的关节角度、角速度、角加速度轨迹

画出末端执行器的笛卡尔坐标轨迹

计算仿真时间 t=3.2s 时的雅可比矩阵



2)设机器人关节空间初始状态为 θ=[100°,90°],末端执行器在笛卡尔空间的线速度指令为（常数）v=[-0.15,-0.25]m/s,用计算力矩控制规则实现机器人末端执行器的线速度指令运动，并利用 Simulink 进行时长 10s 的动力学仿真，并解决如下问题：

画出计算力矩控制规则框图，并进行解释说明

画出机器人的关节力矩、角度、角速度、角加速度轨迹

画出末端执行器笛卡尔坐标轨迹

利用雅可比矩阵计算并画出末端执行器的笛卡尔坐标速度轨迹

画出 XZ 平面（即只考虑 XZ 平面内的速度雅可比矩阵）内 Yoshikawa 可操作性度量值的轨迹（w=sqrt(det[J(q)（J(q)）T]) ），分析 w 和机器人位形的关系
