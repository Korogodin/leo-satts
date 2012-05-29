function Ysolve = fsolve_Kepler(Xs, Xizm, Yizm, Zizm, VXizm, VYizm, VZizm)

global mu_earth

%%%%Xsolve = r, u, Omega, i, d_r, d_u
% Xsolve = theta, omega_p, Omega, i, e, p
theta = Xs(1);
omega_p = Xs(2);
Omega = Xs(3);
i = Xs(4);
e = Xs(5);
p = Xs(6);

Xsolve(1) = p / (1 + e*cos(theta));
Xsolve(2) = theta + omega_p;
Xsolve(3) = Omega;
Xsolve(4) = i;
Xsolve(5) = sqrt(mu_earth/p) * e * sin(theta);
if Xsolve(1) ~= 0
    Xsolve(6) = sqrt(mu_earth/p) * (1 + e*cos(theta)) / Xsolve(1);
else
    Xsolve(6) = 0;
end

[x y z] = ...
        get_vector_XYZ(Xsolve(1), Xsolve(3), Xsolve(4), Xsolve(2));

[Vx, Vy, Vz] = get_vector_VxVyVz( Xsolve(1), Xsolve(5), Xsolve(3), Xsolve(4), Xsolve(2), Xsolve(6), x, y, z);

ErrX = Xizm - x;
ErrY = Yizm - y;
ErrZ = Zizm - z;

ErrVx = VXizm - Vx;
ErrVy = VYizm - Vy;
ErrVz = VZizm - Vz;

Ysolve = [ErrX, ErrY, ErrZ, ErrVx, ErrVy, ErrVz];

end

