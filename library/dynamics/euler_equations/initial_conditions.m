clear Jx Jy Jz wx0 wy0 wz0 w0 J inv_J

Jx = 0.0700; %kgm^2
Jy = 0.0504; %kgm^2
Jz = 0.0109; %kgm^2

wx0 = 0.2; %rad/s 
wy0 = 0.1; %rad/s
wz0 = 0.0; %rad/s

w0 = [wx0 wy0 wz0];

J = diag([Jx Jy Jz]);
inv_J = inv(J);