classdef State_X < handle
    %State_X state vector
    
    properties (SetAccess = public)
        e,
        p,
        theta,
        d_theta,
        omega,
        Omega,
        d_Omega,
        i,
        d_i,
        X, 
        lambda
    end
    
    methods
        
        function XX = State_X()
            XX.X = nan(9, 1);
        end
        
        function set_scalar(XX)
            XX.e = XX.X(1);
            XX.p = XX.X(2);
            XX.theta = XX.X(3);
            XX.d_theta = XX.X(4);
            XX.omega = XX.X(5);
            XX.Omega = XX.X(6);
            XX.d_Omega = XX.X(7);
            XX.i = XX.X(8);
            XX.d_i = XX.X(9);
        end

        function set_vector(XX)
            XX.X(1) = XX.e;
            XX.X(2) = XX.p;
            XX.X(3) = XX.theta;
            XX.X(4) = XX.d_theta;
            XX.X(5) = XX.omega;
            XX.X(6) = XX.Omega;
            XX.X(7) = XX.d_Omega;
            XX.X(8) = XX.i;
            XX.X(9) = XX.d_i;
        end
        
    end
    
end

