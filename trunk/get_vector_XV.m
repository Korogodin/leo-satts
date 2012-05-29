function [x, y, z, Vx, Vy, Vz] = get_vector_XV( e, p, theta, omega, Omega, i)
%GET_VECTOR_XV Get vectors x, y, z and Vx, Vy, Vz

global mu_earth

munapi = sqrt(mu_earth / p);
Vr = munapi*e*sin(theta); 
Vu = munapi*(1+e*cos(theta));
r = p / (1+e*cos(theta));
u = theta + omega; 
lambda = Omega;

xyz = U3(-Omega)*U1(-i)*U3(-(theta+omega))*[r; 0; 0];

x = xyz(1);
y = xyz(2);
z = xyz(3);

Vx = Vr.*x./r - Vu.*(sin(u).*cos(lambda) + cos(u).*sin(lambda).*cos(i));
Vy = Vr.*y./r - Vu.*(sin(u).*sin(lambda) - cos(u).*cos(lambda).*cos(i));
Vz = Vr.*z./r + Vu.*cos(u).*sin(i);


end

