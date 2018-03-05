function [ wIw ] = WCrossIdotW( w ,I )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
% 计算w x I .w
% w为
syms Wx Wy Wz Ixx Iyy Izz Ixy Ixz Iyz
syms I1 I2 I3 Iv

Wx = w(1,1);
Wy = w(2,1);
Wz = w(3,1);

Ixx = I(1,1);
Iyy = I(2,2);
Izz = I(3,3);

Ixy = -I(1,2);
Ixz = -I(1,3); 
Iyz = -I(3,2);


% w = [Wx;Wy;Wz];
% Iv = [I1;I2;I3]
% 
% 
% wIw = cross(w,Iv)

% I1 = [Ixx -Ixy -Ixz];
% I2 = [-Ixy Iyy -Iyz];
% I3 = [-Ixz -Iyz Izz];
% 
% Iv = eval(Iv)
% wIw = eval(wIw)
wIw =  [Wz*(Iyz*Wz + Izz*Wy) - Wy*(Iyz*Wy + Iyy*Wz) - Wx*(Ixz*Wy - Ixy*Wz);
 Wx*(Ixz*Wx + Ixx*Wz) + Wy*(Iyz*Wx - Ixy*Wz) - Wz*(Ixz*Wz + Izz*Wx);
 Wy*(Ixy*Wy + Iyy*Wx) - Wx*(Ixy*Wx + Ixx*Wy) + Wz*(Ixz*Wy - Iyz*Wx)];


end

