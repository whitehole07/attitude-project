% SOLAR RADIATION PRESSURE TORQUE CALCULATION
%
% computes the solar radiation pressure torque applied on a cubic
% spacecraft based on its dimension, orientation, and time of the
% year to compute solar position

% /!\ this is only the first version, aka experimental draft
% TO-DO : We should account for the planet shadowing on the spacecraft for
% the next version, as this model only represents the sun and the satellite
% in the system.

function T = srp(Sat_dim,N,A_BN,t)
    % Sat_dim is the cubic satellite dimentions [x y z]' in m
    % N is the set of faces' orientations (normals) of the satellite
    % given as [n1 n2 n3 n4 n5 n6] with ni = [xi yi zi]'
    % A_BN is the rotation matrix between the satellite and the inertial frame
    % t is the time of the year in s
    
    % reflection coefficients definition
    rho_s = 0.5; % [no unit] specular refelection coefficient
    rho_d = 0.1; % [no unit] diffuse refelection coefficient

    %faces surface definition
    surface = [Sat_dim(2)*Sat_dim(3) Sat_dim(1)*Sat_dim(3) Sat_dim(1)*Sat_dim(2)];
    
    %pressure definition
    c = 299792.458; % [m/s] light speed
    Fe = 1358; % [W/m2] power per unit surface
    P = Fe/c; % [Pa] Solar Radiation Pressure
    
    %sun orientation definition
    T = 365.25*24*3600; % [s] one year's time span
    n_sun = 2*pi/T; % [rad/s] sun's rotation period around earth
    Rs = 150e9; % [m] earth/sun distance
    epsilon = 23.45*pi/180; % [rad]
    S_B = A_BN*[Rs*cos(n_sun*t) Rs*sin(n_sun*t)*cos(epsilon) Rs*sin(n_sun*t)*sin(epsilon)]';
    
    %computation
    T = zeros([3 1]); % SRP-induced torque
    for n=1:6
        if dot(N(:,n),S_B) > 0
            A = surface(mod(n-1,3)+1);
            F = -P*A*dot(S_B,N(:,n))*((1-rho_s)*S_B+(2*rho_s*dot(S_B,N(:,n))+(2/3)*rho_d)*N(:,n));
            T = T + cross((Sat_dim(mod(n-1,3)+1)/2)*N(:,n),F);
        end
    end
end