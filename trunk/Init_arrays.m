SV_GLO_List = [];
hF_cont = 0;

% Xist - true values
Xist.Crs = nan(1, Nmod);
Xist.dn = nan(1, Nmod);
Xist.M0 = nan(1, Nmod);
Xist.Cuc = nan(1, Nmod);
Xist.e = nan(1, Nmod);
Xist.Cus = nan(1, Nmod);
Xist.sqrtA = nan(1, Nmod);
Xist.A = nan(1, Nmod);
% toe = interp1(toe_eph, Ephem(:, 19), tmod, 'linear');
Xist.Cic = nan(1, Nmod);
Xist.Omega = nan(1, Nmod);
Xist.Cis = nan(1, Nmod);
Xist.i0 = nan(1, Nmod);
Xist.Crc = nan(1, Nmod);
Xist.omega = nan(1, Nmod);
Xist.Omega_dot = nan(1, Nmod);
Xist.i_dot = nan(1, Nmod);
Xist.d_lambda = nan(1, Nmod);
Xist.d_i = nan(1, Nmod);
Xist.d_u = nan(1, Nmod);
Xist.d_r = nan(1, Nmod);
Xist.E = zeros(1, Nmod);
Xist.x0 = nan(1, Nmod); 
Xist.y0 = nan(1, Nmod); 
Xist.z0 = nan(1, Nmod);

% Xextr - extrapolation 
Xextr = Xist;

% Xest - estimation
Xest = Xist;
