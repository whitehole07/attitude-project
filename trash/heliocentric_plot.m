clear all
close all
clc







% data
a = 2.4501e+4; %[km]
e = 0.6665; %[-]
i  =30; %[deg]
OM = 0; %[rad]
om = 0; %[rad]
th = 0; %[rad]
muE = astroConstants(13);



% stars northern hemisphere
RAAN_star1 = 	[20 41 25.915];
declination_star1 =	[45 16 49.22]; 

RAAN_star2 = [19 44 58.5]; 
declination_star2 = [45 7 51]; 

RAAN_star3 = [20 22 13.7]; 
declination_star3 = [40 15 24.05]; 



% stars austral hemisphere
% 
% %markeb
% RAAN_star1 =[9 22 0.6];
% declination_star1 = [-55 0 39];
% %aspidiske
% RAAN_star2 =[9 17 5.4];
% declination_star2 = [-59 16 31];
% % Avior
% RAAN_star3 =[8 22 30.8];
% declination_star3 = [ -59 30 34];




epsilon = deg2rad(23.5);
% orbit keplerian parameters
kepEI = [a, e, deg2rad(i), deg2rad(OM), deg2rad(om), deg2rad(th)];
[r0, v0] = kep2car(kepEI(1), kepEI(2), kepEI(3), kepEI(4),kepEI(5),kepEI(6), muE);

AU = astroConstants(2); %km
r_sun = [1, 0, 0]'*AU;
kepEI_sun = [norm(r_sun), 0, epsilon, 0, 0, 0];




options=odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
T = 2*pi*sqrt( a^3/muE );  
tspan = linspace( 0, 2*T, 10000);
% Perform the integration without perturbation
[time, state0]  = ode113 (@(t,y) tbp_ode_J2(t, y, muE,0), tspan, [r0; v0], options );

%compute h, useful direction to plot
r_state = state0( :, 1:3);
v_state = state0 ( :, 4:6);
h = cross(r_state, v_state);
h_norm_vect = vecnorm(h, 2, 2);

%plot orbit for sc orbit and sun apparent motion around Earth
t0 = 0;
f0 = 0;
 [f, n_periods] = Kepler (800*T, kepEI_sun(2), kepEI_sun(1), muE, f0, t0)
deltaTh = f-f0;
stepTh = deltaTh/100;
[X_sun, Y_sun, Z_sun] = plotOrbit(kepEI_sun, muE, deltaTh, stepTh);

%% plot

%plot Earth, orbit, sun apparent (circular) motion
k = 1e-3;
Terra_3D;
plot3 (state0(:, 1), state0(:, 2), state0(:, 3));
hold on 
%plot3 (X_sun, Y_sun, Z_sun, 'Marker', 'o', 'markerSize', 10);
quiver3(0,0,0, k*X_sun(1), k*Y_sun(1), k*Z_sun(1));
quiver3(0,0,0, k*X_sun(end), k*Y_sun(end), k*Z_sun(end));

quiver3(state0(1, 1), state0(1, 2), state0(1, 3), h(1,1),h(1,2), h(1,3), 'autoscale', 'off');



% plot direction of star sensor
[star1_in, star2_in, star3_in] = stars(RAAN_star1, declination_star1,RAAN_star2, declination_star2, RAAN_star3, declination_star3)
k = 1e6;
plot3(k*star1_in(1), k*star1_in(2),  k*star1_in(3), 'marker', 'diamond');
plot3(k*star2_in(1), k*star2_in(2),  k*star2_in(3), 'marker', 'pentagram');
plot3(k*star3_in(1), k*star3_in(2),  k*star3_in(3), 'marker', '*');
 
 
k = 1e6;
quiver3 (state0(1, 1), state0(1, 2), state0(1, 3), k*star3_in(1), k*star3_in(2), k*star3_in(3), 'autoscale', 'off');
 
h_vect =h(1, :);
angle = rad2deg(dot (h_vect/norm(h_vect), star3_in)) % angle between h norm and direction of stars

x_axis_body = cross (h_vect, state0 (1, 1:3)) ;
k = 1e-3;
quiver3 (state0(1, 1), state0(1, 2), state0(1, 3), k*x_axis_body(1),  k*x_axis_body(2),  k*x_axis_body(3), 'autoscale', 'off');
new_direction_of_sensor = h_vect 


legend ('Earth', 'orbit', 'sun initial', 'sun final', 'h', 'star1', 'star2', 'star3', 'direction of stars', 'x_axis');

%%
%  angle = rad2deg(mod(dot (star3_in, h(1, :)), 2*pi));
%  angle1=rad2deg(mod(dot (star3_in, h(1, :)), 2*pi));
%  angle2=rad2deg(mod(dot (star3_in, z_sensor_B),2*pi));

Rot_1 = [cos(angle1), 0, -sin(angle1); 0 1 0; sin(angle1), 0, cos(angle1)];
Rot_2 =[1 0 0;cos(angle2), 0, sin(angle2);  -sin(angle2), 0, cos(angle2)];
z_sensor_B = Rot_2*Rot_1*h(1,3);
quiver3 (0,0,0, k*z_sensor_B(1), k*z_sensor_B(2), k*z_sensor_B(3), 'autoscale', 'off');
alfa = rad2deg(acos(star3_in(3)/norm(star3_in)));
beta = rad2deg(atan (star3_in(1) / (star3_in(2))));
