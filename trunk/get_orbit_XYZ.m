function X = get_orbit_XYZ(X, start_i, stop_i, E0)

global options_solve Nmod

if ~isnan(E0)
    X.E(start_i) = fsolve(@(E)(E-X.e(start_i)*sin(E) - X.M0(start_i)), E0, options_solve);
else
    X.E(start_i) = fsolve(@(E)(E-X.e(start_i)*sin(E) - X.M0(start_i)), X.E(start_i-1), options_solve);
end

for i = (start_i+1):stop_i
    X.E(i) = fsolve(@(E)(E-X.e(i)*sin(E)  - X.M0(i)), X.E(i-1), options_solve);
end

% X.theta1 = atan( (sqrt(1-X.e.^2).*sin(X.E)) ./ (cos(X.E) - X.e) );
X.theta = 2*atan(sqrt((1+X.e)./(1-X.e)).*tan(X.E/2));
for i = 2:Nmod
    while abs(X.theta(i) - X.theta(i-1)) > pi
        if (X.theta(i) - X.theta(i-1)) > pi
            X.theta(i) = X.theta(i) - 2*pi;
        elseif (X.theta(i) - X.theta(i-1)) < -pi
            X.theta(i) = X.theta(i) + 2*pi;
        end
    end
end

phi = X.omega + X.theta;
X.u = phi + X.Cus .* sin(2*phi) + X.Cuc .* cos(2*phi);
% X.u = phi;

X.r = X.A.*(1 - X.e.*cos(X.E)) + X.Crs .* sin(2*phi) + X.Crc .* cos(2*phi);
% X.r = X.A.*(1 - X.e.*cos(X.E));

X.i = X.i0 + X.Cis .* sin(2*phi) + X.Cic .* cos(2*phi);
% X.i = X.i0;

X.lambda = X.Omega;

for i = 1:Nmod
    xyz = U3(-X.lambda(i))*U1(-X.i(i))*U3(-X.u(i))*[X.r(i); 0; 0];
    X.x0(i) = xyz(1); X.y0(i) = xyz(2); X.z0(i) = xyz(3);
end

end

