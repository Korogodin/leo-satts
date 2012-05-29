function F = lsq_i( Xs, start_i, stop_i )
globals;

F = nan(3, stop_i - start_i + 1);

Cic = Xs(1);
Cis = Xs(2);

i0 = Xs(3);
i_dot = Xs(4);
% i0 = Xist.i(1);
% i_dot = (Xist.i(end)-Xist.i(1)) / tmod(end);

u = Xist.u(start_i:stop_i);
d_u = Xist.d_u(start_i:stop_i);
dd_u = Xist.dd_u(start_i:stop_i);
t = tmod(start_i:stop_i) - tmod(start_i);

s_i = i0 + i_dot*t + Cic.*cos(2*u) + Cis.*sin(2*u);
s_d_i = i_dot - 2*Cic*sin(2*u).*d_u + 2*Cis*cos(2*u).*d_u;
s_dd_i =  - 0*4*Cic*cos(2*u).*(d_u).^2 - 2*Cic*sin(2*u).*dd_u ... 
            - 0*4*Cis*sin(2*u).*(d_u).^2 + 2*Cis*cos(2*u).*dd_u;

Err_i = Xist.i(start_i:stop_i) - s_i;
Err_d_i = Xist.d_i(start_i:stop_i) - s_d_i;
Err_dd_i = Xist.dd_i(start_i:stop_i) - s_dd_i;

F(1, :) =  Err_i;
F(2, :) =  Err_d_i;
F(3, :) =  Err_dd_i;

end

