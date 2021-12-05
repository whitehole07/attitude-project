function [X, Y, Z] = plotOrbit(kepEI, mu, deltaTh, stepTh)



% plotOrbit.m  
%
% PROTOTYPE:
% [X, Y, Z] = plotOrbit(kepEI, mu, deltaTh, stepTh)
%
% DESCRIPTION:
% The function saves the positions of the satellite on its
% orbit, given a range of true anomalies and a step to analyze this
% interval.
%
% INPUT:
% kepEI     [1x6]  vector of Keplerian parameters to define the satellite and the orbit      [km, - , rad, rad, rad, rad] 
% mu        [1x1] planetary constant          [km^3*s^-2]
% deltaTh  [1x1]  Angle of orbit's positions to evaluate  [rad]
% stepTh   [1x1]  Angular step to evaluate the orbit   [rad]
% 
% OUTPUT:
% X    [1xstepTh+1]  vector of X coordinates of the satellite during the range of orbit    [km]
% Y    [1xstepTh+1]  vector of Y coordinates of the satellite during the range of orbit    [km]
% Z    [1xstepTh+1]  vector of Z coordinates of the satellite during the range of orbit    [km]


[r, v] = kep2car(kepEI(1), kepEI(2), kepEI(3), kepEI(4), kepEI(5), 0, mu);
X=[r(1)];
Y=[r(2)];
Z=[r(3)];

for j=(0+stepTh):stepTh:deltaTh
    [r, v] = kep2car(kepEI(1), kepEI(2), kepEI(3), kepEI(4), kepEI(5), j, mu);
    X=[X r(1)];
    Y=[Y r(2)];
    Z=[Z r(3)];
end
end