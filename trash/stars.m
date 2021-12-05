function [star1_in, star2_in, star3_in] = stars(RAAN_star1, declination_star1,RAAN_star2, declination_star2, RAAN_star3, declination_star3)


% Function Description:
%
% According to our orbit, the normal to the orbital plane (the angular momentum vector) passes more or less 
% near the north pole, inclined of 30° wrt the vertical. To choose stars,
% we looked at the night sky on stellarium near north pole (Cape North,
% Norway) and looked for stars around our normal (~60° of azimuth, with the
% highest magnitude possible and "near" one each other, in order to avoi
% problems (eg. exiting for error from the sensor's FOV)
%The identifiied stars are Deneb, Delta Cygni and Sadr.
%
%
% Sources for the data:
% --> (https://stellarium-web.org/) to identify stars
% --> wikipedia to find RAAN, declination and distance of the 3 stars
%
% Express RAAN and declination of the stars as follows:
% 
% RAAN = [hours, minutes, seconds]
% declination = [degrees, primes, seconds]
% 
% For our chosen stars use:
% 
% RAAN_star1 = 	[20 41 25.915];
% declination_star1 =	[45 16 49.22]; 
% 
% RAAN_star2 = [19 44 58.5]; 
% declination_star2 = [45 7 51]; 
% 
% RAAN_star3 = [20 22 13.7]; 
% declination_star3 = [40 15 24.05]; 




% Star 1 (Deneb)

TH_star1 = deg2rad( RAAN_star1(1) + RAAN_star1(2)/60 + RAAN_star1(3)/3600)*15; %[rad]
PHI_star1 = deg2rad(declination_star1(1) + declination_star1(2)/60 + declination_star1(3)/3600) ; %[rad]


% Star 2 (Delta Cygni)

TH_star2 = deg2rad( RAAN_star2(1) + RAAN_star2(2)/60 + RAAN_star2(3)/3600)*15; %[rad]
PHI_star2 = deg2rad(declination_star2(1) + declination_star2(2)/60 + declination_star2(3)/3600) ; %[rad]

% Star 3 (Sadr)

TH_star3 = deg2rad( RAAN_star3(1) + RAAN_star3(2)/60 + RAAN_star3(3)/3600)*15; %[rad]
PHI_star3 = deg2rad(declination_star3(1) + declination_star3(2)/60 + declination_star3(3)/3600) ; %[rad]

% compute star directions in cartesian coordinates

matrix_stars = [TH_star1, PHI_star1;...
                TH_star2, PHI_star2; ...
                TH_star3, PHI_star3]; 
            
[Sx, Sy, Sz] = sph2cart(matrix_stars(:,1), matrix_stars(:,2), ones(3,1)); 
star1_in = [Sx(1); Sy(1); Sz(1)]; 
star2_in = [Sx(2); Sy(2); Sz(2)]; 
star3_in =  [Sx(3); Sy(3); Sz(3)]; 

end



