%% CODE TO MAKE ANY SIMULATION 
clc
clear all
close all





%% euler angles initial condition to have coherent attitude rotation matrices at initial instant

%we put different angles in phi and theta to represent the same initial
%attitude for the two representations
Euler0312 = [pi/4, 0, 0]';
Euler0131 = [0 pi/4 0]';

%% omegas in body reference frame

%values for initial omegas
wx =0.2;
wy = 0;
wz = 0;
% wy = 0.1;
% wz = 0.25;
w0 = [wx, wy, wz];


%% spaecraft inertia properties
[I_x_paneldep, I_y_paneldep, I_z_paneldep, I_x_body, I_y_body, I_z_body] = satellite_inertia();

flag_inertia = 1; % put 1 if the spacecraft configuration is with solar panels fully deployed, 0 if the spacecraft is compact, and without solar panels extended.
if flag_inertia == 1
I =  [I_x_paneldep, I_y_paneldep, I_z_paneldep];

else
I  = [I_x_body, I_y_body, I_z_body];

end


J = diag (I);
J_inv = inv(J);


%% orbit characterization

%define various options for the ODE solver
muP= astroConstants(13);
options=odeset('RelTol', 1e-13, 'AbsTol', 1e-14);


% define constants about heliocentric system
n_sun = 2*pi/(365*24*3600);
R_sun = 149600000;
epsilon = 23.45*pi/180;
theta_g0 = 0;



% orbit keplerian parameters assigned
a = 2.4501e+4; %[km]
e = 0.6665; %[-]
i_0  =30; %[deg]
OM_0 = 0; %[rad]
om_0 = 0; %[rad]
true_anomaly0 = 0; %[rad]
muE = astroConstants(13);
[r0,v0] = kep2car(a, e, i_0, OM_0, om_0, true_anomaly0, muE);



% Set time span
% a = 1/( 2/norm(r0) - dot(v0,v0)/muP );      % Semi-major axis [km]
T = 2*pi*sqrt( a^3/muP );                        % Orbital period [1/s]
n = 2*pi/T; % checked: mean angular velocity of elliptical orbit
tspan = linspace( 0, 200*T, 10000);

[time, state]  = ode113 (@(t,y) tbp_ode_J2(t, y, muP,0), tspan, [r0; v0], options );
r = state(:, 1:3);
v = state(:, 4:6);
h = cross (r0, v0);


r_norm= vecnorm (r, 2,2);
r_versor = r ./ (r_norm);

v_norm= vecnorm (v, 2,2);
v_versor = v ./ (v_norm);

h_norm = norm(h);
h_versor = h ./ (h_norm);


Ts = 0.1; %s, sample time of sensors
f = 1/ Ts;
%%  stuff for magnetic torque

%clear we t0 theta_g0 init_mjd2000

we = 7.2916e-05;  % earth rotation speed
t0 = 0;           % init time
theta_g0 = 0;     % theta greenwich at t0
init_mjd2000 = 0; % init mjd2000

igrfSg = fopen('igrfSg.txt');
igrfSh = fopen('igrfSh.txt');

g_mat = cell2mat(textscan(igrfSg, '%f %f %f %f'));
h_mat = cell2mat(textscan(igrfSh, '%f %f %f %f'));


%[Mat_g, Mat_h, Mat_gsv, Mat_hsv]= ghnorm();


%% star tracker data

% stars_ right ascension declination 
matrix_stars = [(18 + 42/60 + 47/3600)*15, 59 + 37/60 + 50/3600;...
                (22 + 28/60          )*15, 57 + 41/60 + 45/3600; ...
                (21 + 08/60 + 52/3600)*15, 38 + 56/60 + 51/3600]; 
            
[Sx, Sy, Sz] = sph2cart(matrix_stars(:,1)*pi/180, matrix_stars(:,2)*pi/180, ones(3,1)); 
star1_in = [Sx(1); Sy(1); Sz(1)]; 
star2_in = [Sx(2); Sy(2); Sz(2)]; 
star3_in =  [Sx(3); Sy(3); Sz(3)]; 

A_SB = [1 0 0; 0 1 0; 0 0 1]; 
% Matrix of rotation between star sensor frame
% and body frame is Identity matrix in such a way that the sensor has as axis of view the same z_body axis


%% perform simulation

out = sim('euler_equation_dave.slx');

%% squeeze outputs relevant 

t = squeeze (out.tout);
omin = (squeeze(out.ominvector))';
ombody = squeeze(out.omega);
Eulerangles = squeeze (out.Eulerangles);


%% plot stl model of spacecraft
TR_new = plot3(0, 0, 0);
data = stlread('satellite_def.stl');
trisurf ( data );
axis equal
axis([-10 10 -10 10 -10 10 ]*1000)


for i = 1 : length(t)

data_new = data.Points * out.A_final(:, :, i);
delete (TR_new);

TR_new = triangulation(data.ConnectivityList, data_new(:, 1),  data_new(:, 2),  data_new(:, 3));
TR_new = trisurf (TR_new);

SAT.EdgeColor = 'k';
SAT.FaceAlpha = 0.5;
colormap parula
  
drawnow

end


%% model geometry of spacecraft and solar panels

%spacecraft main body
n1=[1;0;0]; n4=-n1;
n2=[0;1;0]; n5=-n2;
n3=[0;0;1]; n6=-n3;
%solar panels
n7=[1;0;0]; n8=-n7;
n9=[1;0;0]; n10=-n9; %back


a_body = 1.1; %m along z
b_body = 0.7; %m along y 
c_body = 0.7; %m along x

% body references
% x = c
% y = b
% z = a


Sat_dim = [c_body b_body a_body]; %in this order to be consistent with srp block
Normal = [n1 n2 n3 n4 n5 n6];




%% control torque
Thrust = 1; %N ?? value invented, to be confirmed once we select a datasheet for a constant thrust
arm = [0.35, 0, 0; 
            -0.35, 0, 0; 
             0.35, 0, 0; 
            -0.35, 0, 0; 
             0, 0.35, -0.55;
             0, 0.35, 0.55]; 

force_dir = [0, 0, -1;
                  0, 0, -1;
                  -1, 0, 0;
                  -1, 0, 0;
                  0, -1, 0;
                  0, -1, 0;];


% m ?? values taken as half of the dimensions of the spacecraft 
MIB = 0.1; %s minimum impulse bit, time of actuator  



%% linear observer


A_sys = [eye(3,3), zeros(3,3)]; 
B_sys = [eye(3,3), zeros(3,3)]; 
C_sys = [eye(3,3), zeros(3,3)]; 
D_sys = [eye(3,3), zeros(3,3)];
Q_weight = diag ([1 2 3]);
R_weight  = diag ([1 2 3]);
[P_LQR,G_LQR,L_LQR] = icare(A_sys, B_sys, Q_weight, R_weight, zeros(6,3), eye(6,6), zeros(6,6));

