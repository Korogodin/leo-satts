function F = lsq_i( Xs )
globals;

F = nan(3, Nmod);

Cic = Xs(1);
Cis = Xs(2);
i0 = Xs(3);
i_dot = Xs(4);

s_i = i0 + i_dot*tmod + Cic.*cos(2*Xist.u) + Cis.*sin(2*Xist.u);
s_d_i = i_dot - 2*Cic*sin(2*Xist.u).*Xist.d_u + 2*Cis*cos(2*Xist.u).*Xist.d_u;
s_dd_i =  - 4*Cic*cos(2*Xist.u).*(Xist.d_u).^2 - 2*Cic*sin(2*Xist.u).*Xist.dd_u ... 
            - 4*Cis*sin(2*Xist.u).*(Xist.d_u).^2 + 2*Cis*cos(2*Xist.u).*Xist.dd_u;

Err_i = Xist.i - s_i;
Err_d_i = Xist.d_i - s_d_i;
Err_dd_i = Xist.dd_i - s_dd_i;

F(1, :) =  Err_i;
F(2, :) =  Err_d_i;
F(3, :) =  Err_dd_i;

end

