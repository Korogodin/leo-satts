classdef SV < handle
    %SV Space Vichle class
    
    properties (SetAccess = private)
        p, % Focal parameter
        e, % Eccentricity
        theta, % True anomaly
        Omega, % Longitude of the ascending node
        omega, % Argument of periapsis
        i, % Inclination
        r, % Radius
        x0, % Coordinates 
        y0,
        z0,
        Vx0, % Velocity
        Vy0,
        Vz0,
        Vr0, 
        Vu, 
        ax0, % Acceleration
        ay0,
        az0,
        name, % Name label
        var_name, % Variable name
        SVOrbit
    end
    
    methods
        
        function sv = SV(p, e, theta, Omega, omega, i, Name, var_name)
            globals;
            sv.p = p;
            sv.e = e;
            sv.theta = theta;
            sv.Omega = Omega;
            sv.omega = omega;
            sv.i = i;
            sv.name = Name;
            sv.var_name = var_name;
            sv.CalcXYZ();
            sv.SVOrbit = Orbit(sv);
            SV_GLO_List{length(SV_GLO_List) + 1} = var_name;
        end
        
        function CalcXYZ(sv)
            globals;
            u = sv.theta + sv.omega;

            sv.r = sv.p ./ (1 + sv.e.*cos(u)); % Height

            sv.x0 = sv.r .* (cos(u).*cos(sv.Omega) - sin(u).*sin(sv.Omega).*cos(sv.i));
            sv.y0 = sv.r .* (cos(u).*sin(sv.Omega) - sin(u).*cos(sv.Omega).*cos(sv.i));
            sv.z0 = sv.r.*sin(u).*sin(sv.i);            
            
            sv.Vr0 = sqrt(mu_earth/sv.p)*sv.e*sin(sv.theta);
            sv.Vu = sqrt(mu_earth/sv.p)*(1 + sv.e*cos(sv.theta));
            
            sv.Vx0 = sv.Vr0*sv.x0/sv.r - sv.Vu*(sin(u)*cos(sv.Omega) + cos(u)*sin(sv.Omega)*cos(sv.i));
            sv.Vy0 = sv.Vr0*sv.y0/sv.r - sv.Vu*(sin(u)*sin(sv.Omega) + cos(u)*cos(sv.Omega)*cos(sv.i));
            sv.Vz0 = sv.Vr0*sv.z0/sv.r - sv.Vu*cos(u)*sin(sv.i);
            
        end
        
    end
    
end

