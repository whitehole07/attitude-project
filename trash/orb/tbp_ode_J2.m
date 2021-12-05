function [dy] = tbp_ode_J2 ( ~ , y, muP, n)
% ode_2bp ODE system for the two-body problem corrected with second zonal harmonical J2 (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% muP[1] Gravitational parameter of the primary [L^3/T^2]
% n[1] Parameter /if equal to 1, consider Earth’s oblateness – Effect of J2[-]
% 
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]

%PUT N= 1 IF YOU WANT TO consider the effect of PERTURBATION of 2BP

%Position and velocity
r = y(1:3);
v = y(4:6);

% Distance from the primary
r_norm= norm(r);


if n == 1       %consider the effect of perturbation of 2BP
J = astroConstants(9); %J2 factor for Earth
R = astroConstants(23); %[km] Earth radius


a = 3/2/(r_norm)^4 * (J*muP*(R)^2) * [ r(1)/r_norm*(5* (r(3)/r_norm)^2 -1) 
                                                                    r(2)/r_norm*(5* (r(3)/r_norm)^2 -1)
                                                                    r(3)/r_norm*(5* (r(3)/r_norm)^2 -3)];
    
    dy = [v(1); v(2); v(3);
        -muP/r_norm^3*r(1)+a(1);
       -muP/r_norm^3*r(2)+a(2);
        -muP/r_norm^3*r(3)+a(3)];
    
else
    
     dy = [v(1); v(2); v(3);
        -muP/r_norm^3*r(1);
        -muP/r_norm^3*r(2);
        -muP/r_norm^3*r(3)];
    
end
