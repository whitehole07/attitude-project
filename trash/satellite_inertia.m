function [I_x_paneldep, I_y_paneldep, I_z_paneldep, I_x_body, I_y_body, I_z_body] = satellite_inertia()

% this function aims to define satellite and solar panels inertia matrix in order 
% to be able to perform SRP disturbance torque (and probably other) calculations 

%   the reference satellite is a Starlink prototype (Feb 2018, according to
%   wikipedia), with 2 solar panels

%   References:

%   1) https://docs.rs-online.com/dd20/0900766b8171c04b.pdf 
%   --> to approximate mass and z-dimension of the solar panel. I used the
%       z-dimension proposed dividing it by 2.
%
%   2)http://spaceflight101.com/falcon-9-paz/microsat-2a-2b/
%   --> for body and solar panels x-y dimensions
%
%   3) https://web.archive.org/web/20201117160214/https://space.skyrocket.de/doc_sdat/starlink-v1-0.htm
%   --> for starlink v1.0 's mass

% total mass of a single starlink satellite is approximately 260kg

% PROBLEMS

%    1)   In the end we have mass for the latest (and really used) version of
%         starlink satellite and the design of one of the first designs of the
%         satellite (eg. final version, v1.0, has just 1 solar panel while we have
%         2 in our model. 
%
%           -->   I don't think this is such a big deal cause we don't want
%                 to model a real spacecraft but we just needed some realistic and credible
%                 data to run ou simulation
%
%   2)     We approximately know the overall mass of a starlink satellite
%          (~ 260kg) but we don't know how much is from the 'body' and how much is
%          from the panels. I suggest the following:
%
%          From reference 1) we can obtain the density of the solar panel
%          rho = 9.3/(1.082*0.796*0.035) = 308.5138 kg/m^3
%
%          Then, we can multiply the density obtain per our volume
%          V  = 8*2*0.015 = 0.2400  m^3
%          m_panels = 2*V*rho = 148.0866 kg
%
%          Therefore, we have:
%          m_body = 260 - m_panel = 111.9134 kg

% body
% x = c
% y = b
% z = a

a_body = 1.1; %m
b_body = 0.7; %m
c_body = 0.7; %m
m_body = 111.9134; %kg

%inertia properties of the spacecraft without panel deployment
I_x_body = (m_body/12) * (a_body^2 + b_body^2);
I_y_body = (m_body/12) * (a_body^2 + b_body^2);
I_z_body = (m_body/12) * (b_body^2 + c_body^2);
I_body = [I_x_body,I_y_body,I_z_body];
I_body = diag(I_body);

% solar panels (x1)
%to compare with notation used in class we have
% x = c
% y = b
% z = a

a_panel = 2; %m
b_panel = 8; %m
c_panel = 0.015; %m
m_panel = 148.0866/2; %kg

%inertia properties of the spacecraft of the only solar panels
I_x_panel = (m_panel/12) * (a_panel^2 + b_panel^2);
I_y_panel = (m_panel/12) * (a_panel^2 + b_panel^2);
I_z_panel = (m_panel/12) * (b_panel^2 + c_panel^2);

%inertia properties of the spacecraft with panel fully deployed
I_x_paneldep = I_x_panel + I_x_body;
I_y_paneldep = I_y_panel + I_y_body;
I_z_paneldep = I_z_panel + I_z_body;





%% old code, probably strange definitions of x,y,z and Ix Iy Iz
% % % % body
% % % x_body = 1.1; %m
% % % y_body = 0.7; %m
% % % z_body = 0.7; %m
% % % m_body = 111.9134; %kg
% % % 
% % % I_x_body = (m_body/12) * (x_body^2 + y_body^2);
% % % I_y_body = (m_body/12) * (x_body^2 + y_body^2);
% % % I_z_body = (m_body/12) * (y_body^2 + z_body^2);
% % % 
% % % I_body = [I_x_body,I_y_body,I_z_body];
% % % I_body = diag(I_body);
% % % 
% % % % solar panels (x1)
% % % 
% % % x_panel = 8; %m
% % % y_panel = 2; %m
% % % z_panel = 0.015; %m
% % % m_panel = 148.0866/2; %kg
% % % 
% % % I_x_panel = (m_panel/12) * (x_panel^2 + y_panel^2);
% % % I_y_panel = (m_panel/12) * (x_panel^2 + y_panel^2);
% % % I_z_panel = (m_panel/12) * (y_panel^2 + z_panel^2);
% % % 
% % % I_panel = [I_x_panel,I_y_panel,I_z_panel];
% % % I_panel = diag(I_panel);
% % % 
% % % 

end

