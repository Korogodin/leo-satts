function Ysolve = fsolve_XandV(Xsolve, Xizm, Yizm, Zizm, VXizm, VYizm, VZizm)

%Xsolve = r, u, Omega, i, d_r, d_u

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

