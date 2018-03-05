% 首先定义整个动力学方程的输入
% 有些内容为常量，不过暂时没有，也暂时定义成符号变量。
% 输入为6个关节的角度： q1 - q6 
% 6个关节的角速度： dq1 - dq6
% 6个关节的角加速度 ： ddq1 - ddq6
% 需要计算的中间量 
% 6个连杆的角速度 ： w1- w6
% 6个连杆的角加速度 ： dw1 - dw6 
% % 6个连杆的线速度 ： v1 - v6 
% 6个连杆的线加速度 ： dv1 - dv6
% 6个连杆中心位置的线加速度 ： dvc1 - dvc6
% 6个连杆中心处惯性力 ： F1 - F6
% 6个连杆中心处转矩 ： N1 - N6
% 逆向推导关节转矩的输出量：
% 6个关节上的受力 :  f1 - f6
% 6个关节上的力矩 ： n1 - n6
% 6个关节Z方向的输出转矩 ： torque_1 - torque_6

%计算中间需要的常量 定义 ：
% 6个连杆的质量 ： m1 - m6
% 6个连杆的质心位置向量 ： Pcm_1 - Pcm_6
% 6个连杆的中心点位置的转动惯量 ： I_1 - I_6 ， Ixy_1, Ixx_1,Ixz_1,Iyz_1,Iyy_1,Izz_1
% 6个坐标系之间的变换矩阵 ： R12 - R56 ， 逆向矩阵 R21 - R65
% 6个关节的正弦、余弦值 ： c1 - c6 ,s1 -s6


%定义各变量
syms q1   q2   q3   q4   q5   q6
syms dq1  dq2  dq3  dq4  dq5  dq6
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6

syms w1  w2  w3  w4  w5  w6
syms dw1 dw2 dw3 dw4 dw5 dw6
syms v1   v2   v3   v4   v5   v6
syms dv1  dv2  dv3  dv4  dv5  dv6
syms dvc1 dvc2 dvc3 dvc4 dvc5 dvc6

syms F1 F2 F3 F4 F5 F6
syms N1 N2 N3 N4 N5 N6

syms f1 f2 f3 f4 f5 f6
syms n1 n2 n3 n4 n5 n6
syms torque_1 torque_2 torque_3 torque_4 torque_5 torque_6

%定义需要输入的常量
syms m1 m2 m3 m4 m5 m6
syms Pcm_1 Pcm_2 Pcm_3 Pcm_4 Pcm_5 Pcm_6
syms cm_x_1 cm_x_2 cm_x_3 cm_x_4 cm_x_5 cm_x_6
syms cm_y_1 cm_y_2 cm_y_3 cm_y_4 cm_y_5 cm_y_6
syms cm_z_1 cm_z_2 cm_z_3 cm_z_4 cm_z_5 cm_z_6
syms P1   P2   P3   P4   P5   P6
syms L1 L2 L3 L4 L5 L6
% syms Px_1 Py_1 Pz_1
% syms Px_2 Py_2 Pz_2
% syms Px_3 Py_3 Pz_3
% syms Px_4 Py_4 Pz_4
% syms Px_5 Py_5 Pz_5
% syms Px_6 Py_6 Pz_6

P1 = [0;0;0];
P2 = [L1;0;0];
P3 = [L2;0;L5];
P4 = [L3;0;L4];
P5 = [0;0;0];
P6 = [0;-L6;0];

Pcm_1 = [cm_x_1;cm_y_1;cm_z_1];
Pcm_2 = [cm_x_2;cm_y_2;cm_z_2];
Pcm_3 = [cm_x_3;cm_y_3;cm_z_3];
Pcm_4 = [cm_x_4;cm_y_4;cm_z_4];
Pcm_5 = [cm_x_5;cm_y_5;cm_z_5];
Pcm_6 = [cm_x_6;cm_y_6;cm_z_6];

syms I_1 I_2 I_3 I_4 I_5 I_6
syms Ixy_1 Ixx_1 Ixz_1 Iyz_1 Iyy_1 Izz_1
syms Ixy_2 Ixx_2 Ixz_2 Iyz_2 Iyy_2 Izz_2
syms Ixy_3 Ixx_3 Ixz_3 Iyz_3 Iyy_3 Izz_3
syms Ixy_4 Ixx_4 Ixz_4 Iyz_4 Iyy_4 Izz_4
syms Ixy_5 Ixx_5 Ixz_5 Iyz_5 Iyy_5 Izz_5
syms Ixy_6 Ixx_6 Ixz_6 Iyz_6 Iyy_6 Izz_6

I_1=[Ixx_1 -Ixy_1 -Ixz_1;
     -Ixy_1 Iyy_1 -Iyz_1;
     -Ixz_1 -Iyz_1 Izz_1];
 
I_2=[Ixx_2 -Ixy_2 -Ixz_2;
     -Ixy_2 Iyy_2 -Iyz_2;
     -Ixz_2 -Iyz_2 Izz_2];
 
I_3=[ Ixx_3 -Ixy_3 -Ixz_3;
     -Ixy_3  Iyy_3 -Iyz_3;
     -Ixz_3 -Iyz_3  Izz_3];

I_4=[ Ixx_4 -Ixy_4 -Ixz_4;
     -Ixy_4  Iyy_4 -Iyz_4;
     -Ixz_4 -Iyz_4  Izz_4];
 
I_5=[ Ixx_5 -Ixy_5 -Ixz_5;
     -Ixy_5  Iyy_5 -Iyz_5;
     -Ixz_5 -Iyz_5  Izz_5];
 
I_6=[ Ixx_6 -Ixy_6 -Ixz_6;
     -Ixy_6  Iyy_6 -Iyz_6;
     -Ixz_6 -Iyz_6  Izz_6];
 
%定义旋转矩阵
syms R12 R23 R34 R45 R56
syms R21 R32 R43 R54 R65
syms c1 c2 c3 c4 c5 c6 
syms s1 s2 s3 s4 s5 s6

R12 = [c2 -s2 0;0,0,1;-s2 -c2 0;];
R23 = [c3 -s3 0;s3 c3 0;0 0 1];
R34 = [c4 -s4 0;0 0 1;-s4 -c4 0];
R45 = [c5 -s5 0;0 0 1;s5 c5 0];
R56 = [c6 -s6 0;0 0 1;-s6 -c6 0];

R21 = inv(R12)
R32 = inv(R23)
R43 = inv(R34)
R54 = inv(R45)
R65 = inv(R56)

s1 = sin(q1);
s2 = sin(q2);
s3 = sin(q3);
s4 = sin(q4);
s5 = sin(q5);
s6 = sin(q6);

c1 = cos(q1);
c2 = cos(q2);
c3 = cos(q3);
c4 = cos(q4);
c5 = cos(q5);
c6 = cos(q6);

w1  = dq1*[0;0;1]
dw1 = ddq1*[0;0;1]
dv1 = [0;0;0]
dvc1 = [0;0;0]
F1 = m1*dvc1 
N1 = I_1*dw1 + cross(w1,I_1*w1)

w2 = R12*w1 + dq2*[0;0;1]
dw2 = R12*dw1 + cross(R12*w1,dq2*[0;0;1]) + ddq2*[0;0;1]
dv2 = R12*(cross(dw1,P2) + cross(w1,cross(w1,P2)) + dv1)
tempCross = cross(w2,Pcm_2)
dvc2 = cross(dw2,Pcm_2) + cross(w2,tempCross) + dv2
F2 = m2* dvc2
N2 = I_2*dw2 + cross(w2,I_2*w2)
% 
w3 = R23*w2 + dq3*[0;0;1]
dw3 = R23*dw2 + cross(R23*w2,dq3*[0;0;1]) + ddq3*[0;0;1]
dv3 = R23*(cross(dw2,P3) + cross(w2,cross(w2,P3)) + dv2)
tempCross = cross(w3,Pcm_3)
dvc3 = cross(dw3,Pcm_3) + cross(w3,tempCross) + dv3
F3 = m3* dvc3
N3 = I_3*dw3 + cross(w3,I_3*w3)
% 
% 
w4 = R34*w3 + dq4*[0;0;1]
dw4 = R34*dw3 + cross(R34*w3,dq4*[0;0;1]) + ddq4*[0;0;1]
dv4 = R34*(cross(dw3,P4) + cross(w3,cross(w3,P4)) + dv3)
tempCross = cross(w4,Pcm_4)
dvc4 = cross(dw4,Pcm_4) + cross(w4,tempCross) + dv4
F4 = m4* dvc4
N4 = I_4*dw4 + cross(w4,I_4*w4)
% 
% 
w5 = R45*w4 + dq5*[0;0;1]
dw5 = R45*dw4 + cross(R45*w4,dq5*[0;0;1]) + ddq5*[0;0;1]
dv5 = R45*(cross(dw4,P5) + cross(w4,cross(w4,P5)) + dv4)
tempCross = cross(w5,Pcm_5)
dvc5 = cross(dw5,Pcm_5) + cross(w5,tempCross) + dv5
F5 = m5* dvc5
N5 = I_5*dw5 + cross(w5,I_5*w5)
% 
% 
w6 = R56*w5 + dq6*[0;0;1]
dw6 = R56*dw5 + cross(R56*w6,dq6*[0;0;1]) + ddq6*[0;0;1]
dv6 = R56*(cross(dw5,P6) + cross(w5,cross(w5,P6)) + dv5)
tempCross = cross(w6,Pcm_6)
dvc6 = cross(dw6,Pcm_6) + cross(w6,tempCross) + dv6
F6 = m6* dvc6
N6 = I_6*dw6 + cross(w6,I_6*w6)
% 
f6 = F6
n6 = N6 + cross(Pcm_6,F6)  
torque_6 = n6(3,1)

f5 = R65*f6 + F5
n5 = N5 + R65*n6 + cross(Pcm_5,F5) + cross(P6,R65*f6)
torque_5 = n5(3,1)
% 
% 
f4 = R54*f5 + F4
n4 = N4 + R54*n5 + cross(Pcm_4,F4) + cross(P5,R54*f5)
torque_4 = n4(3,1)

% 
% f3 = R43*f4 + F3
% n3 = N3 + R43*n4 + cross(Pcm_3,F3) + PcrossRdotF(P4,R43,f4)
% torque_3 = n3*[0 0 1]
% 
% 
% f2 = R32*f3 + F2
% n2 = N2 + R32*n3 + cross(Pcm_2,F2) + PcrossRdotF(P3,R32,f3)
% torque_2 = n2*[0 0 1]
% 
% 
% f1 = R21*f2 + F1
% n1 = N1 + R21*n2 + cross(Pcm_1,F1) + PcrossRdotF(P2,R21,f2)
% torque_1 = n1*[0 0 1]


