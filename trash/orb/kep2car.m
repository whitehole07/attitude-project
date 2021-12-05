function [r,v] = kep2car(a, e, i, OM, om, th, mu)
% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE:
% [r, v] = kep2car(a, e, i, OM, om, th, mu)
%
% DESCRIPTION:
% Conversion from Keplerian elements to Cartesian coordinates. Angles in
% radians.
%
% INPUT:
% a  [1x1] Semi-major axis         [km]
% e  [1x1] Eccentricity            [-]
% i  [1x1] Inclination             [rad]
% OM [1x1] RAAN                    [rad]
% om [1x1] Pericentre anomaly      [rad]
% th [1x1] True anomaly            [rad]
% mu [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% r  [3x1] Position vector         [km]
% v  [3x1] Velocity vector         [km/s]


p=a*(1-e^2);
r_mod=p/(1+e*cos(th));

rPF=r_mod*[cos(th); sin(th); 0];
   
vPF=(sqrt((mu)/p))*[-sin(th); e+cos(th); 0];
               
R3_OM=[cos(OM),  sin(OM), 0;...
       -sin(OM), cos(OM), 0;...
       0,        0,       1];
   
R1_i=[1,       0,      0;...
      0,  cos(i), sin(i);...
      0, -sin(i), cos(i)];
   
R3_om=[ cos(om), sin(om),  0;...
       -sin(om), cos(om),  0;...
              0,       0,  1];

Teci_PF=(R3_om)*(R1_i)*(R3_OM);
TPF_eci=Teci_PF';

r=(TPF_eci)*(rPF);
v=(TPF_eci)*(vPF);


end

