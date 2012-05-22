load Ephem.mat

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
Xist.A = Xist.sqrtA.^2;
% toe = interp1(toe_eph, Ephem(:, 19), tmod, 'linear');
Xist.Cic = interp1(toe_eph, Ephem(:, 20), tmod, 'pchip');
Xist.Omega = interp1(toe_eph, Ephem(:, 21), tmod, 'pchip');
Xist.Cis = interp1(toe_eph, Ephem(:, 22), tmod, 'pchip');
Xist.i0 = interp1(toe_eph, Ephem(:, 23), tmod, 'pchip');
Xist.Crc = interp1(toe_eph, Ephem(:, 24), tmod, 'pchip');
Xist.omega = interp1(toe_eph, Ephem(:, 25), tmod, 'pchip');
Xist.Omega_dot = interp1(toe_eph, Ephem(:, 26), tmod, 'pchip');
Xist.i_dot = interp1(toe_eph, Ephem(:, 27), tmod, 'pchip');

Xist.E = zeros(1, Nmod);
Xist.x0 = nan(1, Nmod); Xist.y0 = nan(1, Nmod); Xist.z0 = nan(1, Nmod);

Xist = get_orbit_XYZ(Xist, 1, Nmod, pi);

Xist.d_lambda = diff(Xist.lambda) / dTmod;
Xist.d_lambda(end+1) = Xist.d_lambda(end);
Xist.d_i = diff(Xist.i) / dTmod;
Xist.d_i(end+1) = Xist.d_i(end);
Xist.d_u = diff(Xist.u) / dTmod;
Xist.d_u(end+1) = Xist.d_u(end);
Xist.d_r = diff(Xist.r) / dTmod;
Xist.d_r(end+1) = Xist.d_r(end);

[Xist.d_x0 Xist.d_y0 Xist.d_z0] = ...
    get_vector_VxVyVz(Xist.r, Xist.d_r, Xist.lambda, ...
                      Xist.i, Xist.u, Xist.d_u, ...
                      Xist.x0, Xist.y0, Xist.z0);