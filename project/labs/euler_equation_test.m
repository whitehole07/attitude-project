clear
close all
clc

Jx = 0.0700; %kgm^2
Jy = 0.0504; %kgm^2
Jz = 0.0109; %kgm^2

wx0 = 0.45; %rad/s 
wy0 = 0.52; %rad/s
wz0 = 0.55; %rad/s

w0 = [wx0 wy0 wz0];

J = diag([Jx Jy Jz]);
inv_J = inv(J);

%% 1st validation
lambda = ((Jz - Jx) * wz0)/Jx;

model = sim("euler_equation.slx");
t = model.tout;
wx = model.w.Data(:, 1);
wy = model.w.Data(:, 2);
wz = model.w.Data(:, 3);

wx_a = wx0 * cos(lambda*t) - wy0 * sin(lambda*t);
wy_a = wx0 * sin(lambda*t) + wy0 * cos(lambda*t);
wz_a = ones(1, size(t, 1)) * wz0;

plot(t, wx_a, t, wy_a, t, wz_a)
hold on
plot(t, wx, '--', t, wy, '--', t, wz, '--')
grid on

%% 2nd validation
clear
close all
clc

Jx = 0.0100; %kgm^2
Jy = 0.0500; %kgm^2
Jz = 0.0700; %kgm^2

wx0 = 0; %rad/s 
wy0 = 2*pi; %rad/s
wz0 = 0; %rad/s

w0 = [wx0 wy0 wz0];

J = diag([Jx Jy Jz]);
inv_J = inv(J);

model = sim("euler_equation.slx");
t = model.tout;
wx = model.w.Data(:, 1);
wy = model.w.Data(:, 2);
wz = model.w.Data(:, 3);

plot(t, wx, t, wy, t, wz)
grid on

