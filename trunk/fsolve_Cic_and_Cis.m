function Ys = fsolve_Cic_and_Cis( Xs, i, d_i, dd_i, u, d_u, dd_u )

%Xs = Cic, Cis, i0, idot
Cic = Xs(1);
Cis = Xs(2);
i0 = Xs(3);
i_dot = Xs(4);

s_i = i0 + Cic*cos(2*u) + Cis*sin(2*u);
s_d_i = i_dot - 2*Cic*sin(2*u)*d_u + 2*Cis*cos(2*u)*d_u;
s_dd_i =  - 4*Cic*cos(2*u)*(d_u)^2 - 2*Cic*sin(2*u)*dd_u ... 
            - 4*Cis*sin(2*u)*(d_u)^2 + 2*Cis*cos(2*u)*dd_u;

Err_i = i - s_i;
Err_d_i = d_i - s_d_i;
Err_dd_i = dd_i - s_dd_i;

Ys = [Err_i, Err_d_i, Err_dd_i];

end

