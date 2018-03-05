function [ PRF ] = PcrossRdotF( P,R,F )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
syms Px Py Pz Ixx r11 r12 r13 r21 r22 r23 r31 r32 r33 fx fy fz
syms p R1 R2 R3 f

Px = P(1,1);
Py = P(2,1);
Pz = P(3,1);

r11 = R(1,1);
r12 = R(1,2);
r13 = R(1,3);

r21 = R(2,1);
r22 = R(2,2);
r23 = R(2,3);

r31 = R(3,1);
r32 = R(3,2);
r33 = R(3,3);

PRF=[fx*(Py*r31 - Pz*r21) + fy*(Py*r32 - Pz*r22) + fz*(Py*r33 - Pz*r23);
    - fx*(Px*r31 - Pz*r11) - fy*(Px*r32 - Pz*r12) - fz*(Px*r33 - Pz*r13);
    fx*(Px*r21 - Py*r11) + fy*(Px*r22 - Py*r12) + fz*(Px*r23 - Py*r13)
    ];
end

