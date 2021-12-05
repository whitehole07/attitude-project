clear we t0 theta_g0 init_mjd2000

we = 7.2916e-05;  % earth rotation speed
t0 = 0;           % init time
theta_g0 = 0;     % theta greenwich at t0
init_mjd2000 = 0; % init mjd2000

igrfSg = fopen('igrfSg.txt');
igrfSh = fopen('igrfSh.txt');

g_mat = cell2mat(textscan(igrfSg, '%f %f %f %f'));
h_mat = cell2mat(textscan(igrfSh, '%f %f %f %f'));