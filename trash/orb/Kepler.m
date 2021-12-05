function [f, n_periods] = Kepler (t, e, a, muP, f0, t0)


% Kepler equation solver for the two-body problem
%
% PROTOTYPE
% [f, n_periods] = Kepler (t, e, a, muP, f0, t0)
%
% INPUT:
% t[1] Delta Time of motion [T]
% e[1] Eccentricity of orbit [-]
% a[1] Semimajor axis of the orbit [L]
% muP[1] Gravitational parameter of the primary [L^3/T^2]
% f0 [1] Initial True anomaly [rad]
% t0 [1] Initial Time [T]
%
% OUTPUT:
% F [1] True anomaly at time t [rad]
%n_periods [1] number of full periods run by the body in time t [-]


% calculate conditions for eccentric anomaly and mean anomaly, at inital
% and final state
E0 = 2 * atan ( sqrt( (1-e)/ (1+e)) * tan (f0 / 2) );
M0 = E0 - e*sin(E0);
M = M0 + sqrt (muP/(a^3) ) * (t - t0);

%split M in number of full orbits (n_periods) and Mrest, which is 0 <
%Mrest < 2pi
Mrest  = wrapTo2Pi(M);
n_periods = floor (M/(2*pi));


% solve non linear equation of Keplerian motion
Eguess = Mrest + (e*sin(Mrest)) / (1 - sin(e+Mrest) + sin(Mrest) ) ; 
fun = @(E) E - e*sin(E) - Mrest;
options=optimset('TolFun', 1e-09);
Esol = fzero(fun, Eguess, options); 

%calculate final true anomaly, paired with found Esol
f = 2* atan  ( sqrt( (1+e) / (1-e)) * tan (Esol / 2) );   
f  = wrapTo2Pi(f);
frad = rad2deg(f);

end







