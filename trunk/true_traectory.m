load Ephem.mat

fprintf('True traectory calculation:\n');
fprintf('- read and interpolation...');
% t_eph = Ephem(:, 4)*24*60*60 + Ephem(:, 5)*60*60 + Ephem(:, 6)*60 + + Ephem(:, 7)*60;
% t_eph = t_eph - t_eph(1);
toe_eph = Ephem(:, 19); toe_eph = toe_eph - toe_eph(1);
Xist.Crs = interp1(toe_eph, Ephem(:, 12), tmod, 'pchip');
Xist.dn = interp1(toe_eph, Ephem(:, 13), tmod, 'pchip');
for i = 2:length(toe_eph)
    while abs(Ephem(i, 14) - Ephem(i-1, 14)) > pi
        if (Ephem(i, 14) - Ephem(i-1, 14)) > pi
            Ephem(i, 14) = Ephem(i, 14) - 2*pi;
        elseif (Ephem(i, 14) - Ephem(i-1, 14)) < -pi
            Ephem(i, 14) = Ephem(i, 14) + 2*pi;
        end
    end
end
Xist.M0 = interp1(toe_eph, Ephem(:, 14), tmod, 'pchip');
Xist.Cuc = interp1(toe_eph, Ephem(:, 15), tmod, 'pchip');
Xist.e = interp1(toe_eph, Ephem(:, 16), tmod, 'pchip');
Xist.Cus = interp1(toe_eph, Ephem(:, 17), tmod, 'pchip');
Xist.sqrtA = interp1(toe_eph, Ephem(:, 18), tmod, 'pchip');
N_r =  26 / (1+6.5);
Xist.A = Xist.sqrtA.^2 / N_r;
% toe = interp1(toe_eph, Ephem(:, 19), tmod, 'linear');
Xist.Cic = interp1(toe_eph, Ephem(:, 20), tmod, 'pchip');
Xist.Omega = interp1(toe_eph, Ephem(:, 21), tmod, 'pchip');
Xist.Cis = interp1(toe_eph, Ephem(:, 22), tmod, 'pchip');
Xist.i0 = interp1(toe_eph, Ephem(:, 23), tmod, 'pchip');
Xist.Crc = interp1(toe_eph, Ephem(:, 24), tmod, 'pchip');
Xist.omega = interp1(toe_eph, Ephem(:, 25), tmod, 'pchip');
Xist.Omega_dot = interp1(toe_eph, Ephem(:, 26), tmod, 'pchip');
Xist.i_dot = interp1(toe_eph, Ephem(:, 27), tmod, 'pchip');
fprintf('complete\n');

fprintf('- gradients and coordinates...');
Xist.E = zeros(1, Nmod);
Xist.x0 = nan(1, Nmod); Xist.y0 = nan(1, Nmod); Xist.z0 = nan(1, Nmod);

Xist = get_orbit_XYZ(Xist, 1, Nmod, pi);

dTmod = dTmod / (N_r^(3/2));

Xist.d_omega = diff(Xist.omega) / dTmod;
Xist.d_omega(end+1) = Xist.d_omega(end);
Xist.dd_omega = diff(Xist.d_omega) / dTmod;
Xist.dd_omega(end+1) = Xist.dd_omega(end);

Xist.d_Omega = diff(Xist.Omega) / dTmod;
Xist.d_Omega(end+1) = Xist.d_Omega(end);
Xist.dd_Omega = diff(Xist.d_Omega) / dTmod;
Xist.dd_Omega(end+1) = Xist.dd_Omega(end);
Xist.d_lambda = diff(Xist.lambda) / dTmod;
Xist.d_lambda(end+1) = Xist.d_lambda(end);

Xist.d_i = diff(Xist.i) / dTmod;
Xist.d_i(end+1) = Xist.d_i(end);
Xist.dd_i = diff(Xist.d_i) / dTmod;
Xist.dd_i(end+1) = Xist.dd_i(end);
Xist.i0 = Xist.i - Xist.Cic.*cos(2*Xist.u) - Xist.Cis.*sin(2*Xist.u);
Xist.i_dot = diff(Xist.i0) / dTmod;
Xist.i_dot(end+1) = Xist.i_dot(end);

Xist.d_theta = diff(Xist.theta) / dTmod;
Xist.d_theta(end+1) = Xist.d_theta(end);
Xist.dd_theta = diff(Xist.d_theta) / dTmod;
Xist.dd_theta(end+1) = Xist.dd_theta(end);

Xist.d_u = diff(Xist.u) / dTmod;
Xist.d_u(end+1) = Xist.d_u(end);
Xist.dd_u = diff(Xist.d_u) / dTmod;
Xist.dd_u(end+1) = Xist.dd_u(end);

Xist.d_r = diff(Xist.r) / dTmod;
Xist.d_r(end+1) = Xist.d_r(end);
Xist.dd_r = diff(Xist.d_r) / dTmod;
Xist.dd_r(end+1) = Xist.dd_r(end);

Xist.d_A = diff(Xist.A) / dTmod;
Xist.d_A(end+1) = Xist.d_A(end);
Xist.dd_A = diff(Xist.d_A) / dTmod;
Xist.dd_A(end+1) = Xist.dd_A(end);

Xist.d_e = diff(Xist.e) / dTmod;
Xist.d_e(end+1) = Xist.d_e(end);
Xist.dd_e = diff(Xist.d_e) / dTmod;
Xist.dd_e(end+1) = Xist.dd_e(end);

Xist.p = Xist.A.*(1 - Xist.e.^2);
Xist.d_p = diff(Xist.p) / dTmod;
Xist.d_p(end+1) = Xist.d_p(end);
Xist.dd_p = diff(Xist.d_p) / dTmod;
Xist.dd_p(end+1) = Xist.dd_p(end);

[Xist.d_x0 Xist.d_y0 Xist.d_z0] = ...
    get_vector_VxVyVz(Xist.r, Xist.d_r, Xist.lambda, ...
                      Xist.i, Xist.u, Xist.d_u, ...
                      Xist.x0, Xist.y0, Xist.z0);
% 
% Xist.x0 = Xist.x0 / N_r;
% Xist.y0 = Xist.y0 / N_r;
% Xist.z0 = Xist.z0 / N_r;

% Xist.d_x0 = diff(Xist.x0) / dTmod;
% Xist.d_x0(end+1) = Xist.d_x0(end);
% Xist.d_y0 = diff(Xist.y0) / dTmod;
% Xist.d_y0(end+1) = Xist.d_y0(end);
% Xist.d_z0 = diff(Xist.z0) / dTmod;
% Xist.d_z0(end+1) = Xist.d_z0(end);

tmod = tmod / (N_r^(3/2));


fprintf('complete\n');
