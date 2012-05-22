function [Vx, Vy, Vz] = get_vector_VxVyVz( r, d_r, lambda, i, u, d_u, x, y, z)
%GET_VECTOR_XYZ Get vector Vx, Vy, Vz

Vu = r.*d_u;
Vr = d_r;

Vx = Vr.*x./r - Vu.*(sin(u).*cos(lambda) + cos(u).*sin(lambda).*cos(i));
Vy = Vr.*y./r - Vu.*(sin(u).*sin(lambda) - cos(u).*cos(lambda).*cos(i));
Vz = Vr.*z./r + Vu.*cos(u).*sin(i);


end

