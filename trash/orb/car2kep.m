function [a, mod_e, i, OM, om, th] = car2kep(r, v, mu)
% car2kep.m - Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE:
% [a, mod_e, i, OM, om, th] = car2kep(r, v, mu)
%
% DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements. Angles in
% radians.
%
% INPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]
% mu [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% a [1x1] Semi-major axis [km]
% mod_e [1x1] Eccentricity [-]
% i [1x1] Inclination [rad]
% OM [1x1] RAAN [rad]
% om [1x1] Pericentre anomaly [rad]
% th [1x1] True anomaly [rad]




%% modulo r e v
    mod_r=norm(r);
    mod_v=norm(v);
   
%% h
    h=cross(r,v);
    mod_h=norm(h);
%% i
    
    k=[0;0;1];
    i=acos((h(3,1))/mod_h);
    
%% eccenticita'

    e=(1/mu)*((mod_v^2-mu/mod_r)*r-dot(r,v)*v);
    mod_e=norm(e);
    
%% energia

    eps=(0.5*mod_v^2)-(mu/mod_r);
    a=-mu/(2*eps);
    
%% linea dei nodi
    
    N=cross(k,h);
    mod_N=norm(N);
    
%% Ascensione Retta del Nodo Ascendente
    
    if(N(2,1)>=0)
        OM=acos(N(1,1)/mod_N);
    else
        OM=2*pi-acos(N(1,1)/mod_N);
    end  
%% anomalia pericentro

    if (e(3,1)>=0)
        om=acos(dot(N,e)/(mod_N*mod_e));
    else
        om=2*pi-acos(dot(N,e)/(mod_N*mod_e));
    end
    
%% velocità radiale 

    v_r=dot(r,v)/mod_r;

%% anomalia vera 

    if(v_r>=0)
        th=acos(dot(e,r)/(mod_e*mod_r));
    else
        th= 2*pi - acos(dot(e,r)/(mod_e*mod_r));
    end
     
    
  


end

