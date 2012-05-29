function Ysolve = fsolve_e_p_theta(Xs, Vr_izm, Vu_izm, r_izm)

global mu_earth

% Xs = e, p, theta
e = Xs(1);
p = Xs(2);
theta = Xs(3);

r = p / (1+e*cos(theta));
Vr = sqrt(mu_earth/p)*e*sin(theta);
Vu = sqrt(mu_earth/p)*(1 + e*cos(theta));

Ysolve = [Vr_izm - Vr, Vu_izm - Vu, r_izm - r];

end

