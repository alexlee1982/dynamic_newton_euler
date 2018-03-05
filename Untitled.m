% syms w1 wx1 wy1 wz1
% syms I_1 I_2 I_3 I_4 I_5 I_6
% syms Ixy_1 Ixx_1 Ixz_1 Iyz_1 Iyy_1 Izz_1
% 
% I_1=[Ixx_1 -Ixy_1 -Ixz_1;
%      -Ixy_1 Iyy_1 -Iyz_1;
%      -Ixz_1 -Iyz_1 Izz_1];
% 
% w1 = [wx1;wy1;wz1]
% 
% w = WCrossIdotW(w1,I_1)
% w = expand(w)
syms Px Py Pz Ixx r11 r12 r13 r21 r22 r23 r31 r32 r33 fx fy fz
syms p R1 R2 R3 f Rv PRf

p = [Px;Py;Pz];
f = [fx;fy;fz];
% 
% R1 = [r11 r12 r13];
% R2 = [r21 r22 r23];
% R3 = [r31 r32 r33];

Rv = [R1;R2;R3];
PRf = cross(p,Rv)


R1 = [r11 r12 r13];
R2 = [r21 r22 r23];
R3 = [r31 r32 r33];

Rv = eval(Rv)
PRf = eval(PRf)
PRf = PRf*f
