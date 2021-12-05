function [a_sol, a_no_J2] = animations_repeated_groundTrack (e, i, om, OM, f0, J2, m, k, theta_g0, muP, we, t0,param)


%no=0/yes=1 J2 perturbation
%compute and plots repeating groundtracks, considering or not J2 effect

%% constants

Te = 86164; % Earth period     %[km^3/s^2]
f = m / k ; % rate 


%% computing a without J2 effect


T = f * Te; % satellite period
tf = k*T;
a_no_J2 = (muP*(T/(2*pi))^2) ^ (1/3);
a = a_no_J2;

%% computing a with J2 effect
if J2 == 1


fun = @(x) (we - (-(3/2*(sqrt(muP)*astroConstants(9)*astroConstants(23)^2)/((1-e^2)^2*x^(7/2)))*cos(i))) / (1/(sqrt( x^3/muP )) + (-(3/2*(sqrt(muP)*astroConstants(9)*astroConstants(23)^2)/((1-e^2)^2*x^(7/2)))*(5/2*sin(i)*sin(i)-2)) + (-(3/2*(sqrt(muP)*astroConstants(9)*astroConstants(23)^2)/((1-e^2)^(3/2)*x^(7/2)))*(1-3/2*sin(i)*sin(i)))) - f*(1+0*x);
options=optimset('TolFun', 1e-09);
a_guess = a_no_J2;
a_sol = fzero(fun, a_guess, options); 
a = a_sol;
T = 2*pi * sqrt(a^3/muP); 
tf = k*T;
flag = 1
else
    a_sol = a_no_J2;

end

[r0,v0] = kep2car(a, e, i, OM, om, f0, muP);

%% performing groundtrack concerning the repetition limits?
% ground track function parameters


[alfa, delta, lon, lat] = groundTrack (r0, v0, theta_g0, muP, we, t0, tf, J2, param);

%% ODE solving
% put n = 0 for perform integration without J2 effect, put n = 1 for J2 perturbation
tspan = linspace( t0, tf, param); 
options=odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
n= J2;
[time, state]  = ode113 (@(t,y) tbp_ode_J2(t, y, muP,n), tspan, [r0; v0], options );
% define r and v
r = state (:, 1 : 3);
v = state (:, 4 : 6);

X= r(:, 1);
Y = r(:, 2);
Z = r(:, 3);


%% plot of longitude and latitude
figure (1)
plot (wrapTo180(rad2deg(lon)) , rad2deg(lat), '.');
hold on
plot (wrapTo180(rad2deg(lon(1))) , rad2deg(lat(1)),'Color', 'g', 'Marker', 'o');
plot (wrapTo180(rad2deg(lon(end))) , rad2deg(lat(end)),'Color',  'g', 'Marker', 's');
xlim([-180 180])
ylim([-90 90])
% axis tight

% Gonzalone = 'gonzalo.jpg';
Earth_image = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
a = imread( Earth_image);
h = image(xlim, -ylim, a); 
uistack(h,'bottom')
xlabel ('Longitude [deg]')
ylabel ('Latitude [deg]')
title('Groundtrack of orbit')
legend ('groundtrack', 'start', 'end')
yticks([-90: 30 : 90]);
xticks([-180: 30 : 180]); 
axis([-180 180, -90 90])

 %% animation of spacecraft along groudtrack
 % COMMENTED BECAUSE IT IS REALLY HEAVY AND MAKES ALL THE CODE SLOWER
% h = plot(nan, nan, 'or');
% for j = 1 : length (X) 
% delete (h) ;
% h = plot (wrapTo180(rad2deg(lon(j))) , rad2deg(lat(j)), 'or');
% drawnow
% end

%% animation of orbit
% figure(2);
% Terra_3D;
% plot3(X, Y, Z, '--');
% h=plot3(nan,nan,nan,'or');
%     
% for j = 1:length (X)
%     animatedline(X, Y, Z, 'Color', [0 0.4470 0.7410], 'linewidth', 3);
%     set(h,'XData',X(j),'YData',Y(j),'ZData',Z(j));
%     drawnow
% end



end

