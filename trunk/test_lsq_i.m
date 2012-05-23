opti = optimset('Display','iter','TolFun',1e-26, 'TolX',1e-36);
Xs = lsqnonlin(@(x)(lsq_i(x)), [0 0 0 0], ...
                [-1e-5 -1e-5 -pi -1e-8], [+1e-5 +1e-5 +pi +1e-8], ...
                opti);

figure; 
plot(tmod, Xist.i, ...
    tmod, Xs(3) + tmod*Xs(4) + Xs(1)*cos(2*Xist.u) + Xs(2)*sin(2*Xist.u), ...
    tmod, Xist.i0 + Xist.Cic.*cos(2*Xist.u) + Xist.Cis.*sin(2*Xist.u))
