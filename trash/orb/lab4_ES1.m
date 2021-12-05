%% CASE 1


clc 
clear
close all


% Orbit characterization
a = 7171.010; %[km]
e = 0;     %[-]
i = deg2rad(30);    %[rad]
OM = deg2rad(0);  %[rad]
om = deg2rad(0);   %[rad]
f0 = deg2rad(0);  %[rad]  
muP =  astroConstants(13);      %[km^3/s^2]
T = 2*pi*sqrt( a^3/muP );        % Orbital period [s]

[r0,v0] = kep2car(a, e, i, OM, om, f0, muP);



% ground track function parameters
theta_g0 = 0;
we = deg2rad(15.04 /3600);
J2 = 1; % it means no J2 perturbation
t0  = 0;
tf = 5*T;
param = 1000;
%perform ground track calculus
[alfa, delta, lon, lat] = groundTrack (r0, v0, theta_g0, muP, we, t0, tf,J2, param);


%% ODE solving
% put n = 0 for perform integration without J2 effect, put n = 1 for J2 perturbation
tspan = linspace( t0, tf, param); 
options=odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
n= J2;
[time, state]  = ode113 (@(t,y) tbp_ode_J2(t, y, muP,n), tspan, [r0; v0], options );
% define r and v
r = state (:, 1 : 3);
v = state (:, 4 : 6);

X= r(:, 1);
Y = r(:, 2);
Z = r(:, 3);




%% plot of longitude and latitude
figure (1)
plot (wrapTo180(rad2deg(lon)) , rad2deg(lat), '.');
hold on
xlim([-180 180])
ylim([-90 90])
% axis tight

% Gonzalone = 'gonzalo.jpg';
Earth_image = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
a = imread( Earth_image);
h = image(xlim, -ylim, a); 
uistack(h,'bottom')

yticks([-90: 30 : 90]);
xticks([-180: 30 : 180]); 
axis([-180 180, -90 90])

%% animation of spacecraft along groudtrack
h = plot(nan, nan, 'or');
for j = 1 : length (X) 
delete (h) ;
h = plot (wrapTo180(rad2deg(lon(j))) , rad2deg(lat(j)), 'or');
drawnow
end





%% animation of orbit

Terra_3D;
plot3(X, Y, Z, '--');
h=plot3(nan,nan,nan,'or');
    
for j = 1:length (X)
    animatedline(X, Y, Z, 'Color', [0 0.4470 0.7410], 'linewidth', 2);
    set(h,'XData',X(j),'YData',Y(j),'ZData',Z(j));
    drawnow
end
