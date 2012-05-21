classdef Orbit
    
    properties (SetAccess = private)
        X0,
        Y0,
        Z0,
        r,
        theta,
        E,
        T,
        NT,
        dT,
        t
    end
    
    methods
        function orb = Orbit(sv)
            orb.T = 12*60*60;
            orb.NT = 10000;
            orb.dT = orb.T / orb.NT;
            orb.t = orb.dT:orb.dT:orb.T;
            orb.E = orb.t/orb.T*2*pi;
            orb.theta = 2*atan(sqrt((1+sv.e)./(1-sv.e)).*tan(orb.E/2));
            orb.theta = orb.theta + 2*pi*(orb.theta<0);
            orb = orb.GetOrbit(sv);
        end
        
        function orb = GetOrbit(orb, sv)

            u = orb.theta + sv.omega;

            orb.r = sv.p ./ (1 + sv.e.*cos(u)); % Height

            orb.X0 = sv.r .* (cos(u).*cos(sv.Omega) - sin(u).*sin(sv.Omega).*cos(sv.i));
            orb.Y0 = sv.r .* (cos(u).*sin(sv.Omega) - sin(u).*cos(sv.Omega).*cos(sv.i));
            orb.Z0 = sv.r.*sin(u).*sin(sv.i);  
        end
    end
    
end

