function [alfa, delta, lon, lat] = groundTrack (r0, v0, theta_g0, muP, we, t0, tf, J2, param)
% groundTrack for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% [alfa, delta, lon, lat] = groundTrack (r0, v0, theta_g0, muP, we, t0, tf, n)
%
% INPUT:
% r0                [3x1] Initial position in Cartesian coordinates             [km]
% v0                [3x1] Initial velocity in Cartesian coordinates             [km/s]
% theta_g0       [1] initial true anomaly of Greenwich meridian           [rad]
% muP             [1] Gravitational parameter of the primary                 [km^3/s^2]
% we                [1] Angular velocity of Earth along its axis                 [rad/s]
% t0                 [1] Initial time to evaluate groundtrack                     [s]
% tf                  [1] Final time to evaluate groundtrack                      [s]
% J2                  [1] parameter to define type of ODE solving: with J2 effect or without J2 effect [-]
% param           [1] parameter to choose the number of time for the ODE
%                           solving, the animation of the orbit and the
%                           groundtrack. Desirable: upper than 1e+3.                                                [-]
% OUTPUT:
% alfa              [1] right ascension of the spacecraft with respect to Earth Center                  [rad]
% delta            [1] declination of the spacecraft with respect to Earth Center                        [rad]
% lon               [n_timex1] Longitude vector                                                                      [rad]
% lat               [n_timex1]  Latitude vector                                                                         [rad]

%PUT J2= 1 IF YOU WANT TO consider the effect of PERTURBATION of 2BP

% NOTE: If the orbit in the main script is characterized in Keplerian
% coordinates, before calling groundTrack.m USE kep2car.m to convert the
% parameters from Keplerian [a, e, i, OM, om, f] to Cartesian [r0, v0]






% 
% % controllare dentro o fuori?
% 
% [r,v] = kep2car(a, e, i, OM, om, th, muP);


% ODE solving to get state vector
% put n = 0 for perform integration without J2 effect, put n = 1 for J2 perturbation
tspan = linspace( t0, tf, param); 
options=odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
n = J2;
[time, state]  = ode113 (@(t,y) tbp_ode_J2(t, y, muP,n), tspan, [r0; v0], options );



% define r and v
r = state (:, 1 : 3);
v = state (:, 4 : 6);

x= r(:, 1);
y = r(:, 2);
z = r(:, 3);

r_norm  = vecnorm (r, 2,2);




% declination

delta = asin (z ./ r_norm);

% right ascension
alfa = acos ((x./r_norm)./(cos(delta))) .* sign (y./r_norm) + 2*pi*(y./r_norm < 0);


% % right ascension con for
% for i = 1 :  c
% 
% if (y./r_norm > 0)
% 
% 
% 
% % else
% %      alfa = 2*pi - acos ((x./r_norm)./(cos(delta)));
% % 
% % end

%define theta Greenwich
t = linspace (t0, tf, param); % occhio che l'incremento di tempo per ora è arbitrario e fisso
theta_g = [];

for i = 1 : length (t)

theta_g_i = theta_g0 + we * (t(i) - t0);
theta_g = [theta_g; theta_g_i];

end


% longitude

lambda  = alfa - theta_g; % controllare che le dimensioni abbiano senso (che siano o entrambi riga o entrambi colonna)
lon = lambda;

%latitude

phi = delta;
lat = phi;


end



