clear
close all
clc

%% Constants
mu = astroConstants(13);

r0 = [1500 0 0];
v0 = [0 15 0];

%% Orbit1
odeOptions = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[t, y] = ode113(@(t,y) odeTwoBp(t, y, mu, 0, 0), [0 100], [r0 v0], odeOptions);
r = y(:, 1:3);
v = y(:, 4:end);

%% Orbit2
simout = sim("model.slx");
rs = simout.rr.data;
vs = simout.vv.data;