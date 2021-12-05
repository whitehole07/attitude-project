%% CASE 2 REPEATING GROUNDTRACK COMPUTATION, AND OPTIONAL PLOT AND ANIMATION

clc 
clear
close all

%compute repeating groundtrack

%% Orbit characterization
% a unknown because it has to be calculated with the period constraint,
J2 = 0; % it means: no=0/yes=1 J2 perturbation
we = deg2rad(15.04 /3600);
e = 0;     %[-]
i = deg2rad(98);    %[rad]
OM = deg2rad(0);  %[rad]
om = deg2rad(40);   %[rad]
f0 = deg2rad(0);  %[rad]  
muP =  astroConstants(13);      %[km^3/s^2]
k = 15; % satellite revolutions
m = 1; % Earth revolutions
theta_g0 = 0;
t0  = 0;
param = 10000;


[a_sol, a_no_J2] = animations_repeated_groundTrack (e, i, om, OM, f0, J2, m, k, theta_g0, muP, we, t0, param);


