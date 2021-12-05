clc
clear
close all



wy = 0.1;
wz = 0.1;


Ix= 0.07;
Iy= 0.0504;
Iz = 0.0109;

%the following omegas are in body frame
wx_vector = [0.2 : 0.01 : 2*pi]';
wy_vector = ones(length(wx_vector),1)*wy;
wz_vector = ones(length(wx_vector),1)*wz;

%the following h is in body frame
h = [Ix*wx_vector  Iy*wy_vector Iz*wz_vector];

h_in_target = [1 0 0]';









