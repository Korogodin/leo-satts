opti = optimset('Display','off','TolFun',1e-15, 'TolX',1e-15);

short = 0;
if ~short
    lb = [-1e-5 -1e-5 -pi -1e-8];
    ub = [+1e-5 +1e-5 +pi +1e-8];
    X0 = [0 0 0 0];
else
    lb = [-1e-5 -1e-5];
    ub = [+1e-5 +1e-5];
    X0 = [0 0];
end
Xs = lsqnonlin(@(x)(lsq_i(x, 1, 2)), X0, lb, ub, opti);


Cic = Xs(1)*ones(1, Nmod);            
Cis = Xs(2)*ones(1, Nmod);
if ~short
    i0 = Xs(3)*ones(1, Nmod);            
    i_dot = Xs(4)*ones(1, Nmod);
else
    i0 = Xist.i(1)*ones(1, Nmod);
    i_dot = (Xist.i(end)-Xist.i(1)) / tmod(end)*ones(1, Nmod);    
end
   

BW_i = 0.001;
Ki = get_K2(BW_i, dTmod);
for i = 2:Nmod
    if i < Nmod-120
        Xs = lsqnonlin(@(x)(lsq_i(x, i, i+60)), Xs, lb, ub, opti);
    else
        Xs = lsqnonlin(@(x)(lsq_i(x, i-1, i)), Xs, lb, ub, opti);
    end
                
    Cic(i) = Xs(1);
    Cis(i) = Xs(2);
    if ~short
        i0(i) = Xs(3);
        i_dot(i) = Xs(4);
    end
    
    if ~mod(i, fix(Nmod/10))
        fprintf('Done %.0f %% \n', i/Nmod*100);
    end
end

hF = 0;
hF = figure(hF+1);
plot(tmod, Xist.i, ...
    tmod, i0 + Cic.*cos(2*Xist.u) + Cis.*sin(2*Xist.u), ...
    tmod, Xist.i0 + Xist.Cic.*cos(2*Xist.u) + Xist.Cis.*sin(2*Xist.u));
%     tmod, i0 + tmod.*i_dot + Cic.*cos(2*Xist.u) + Cis.*sin(2*Xist.u), ...

hF = figure(hF+1);
subplot(2,2,1); plot(tmod, Xist.Cic, tmod, Cic)
ylabel('Cic, rad');
subplot(2,2,2); plot(tmod, Xist.Cis, tmod, Cis)
ylabel('Cis, rad');
subplot(2,2,3); plot(tmod, Xist.Cic - Cic)
ylabel('\delta Cic, rad');
subplot(2,2,4); plot(tmod, Xist.Cis - Cis)
ylabel('\delta Cis, rad');

hF = figure(hF+1);
subplot(2,2,1); plot(tmod, Xist.i0, tmod, i0)
ylabel('i0, rad');
subplot(2,2,2); plot(tmod, Xist.i_dot, tmod, i_dot)
ylabel('i dot, rad/s');
subplot(2,2,3); plot(tmod, Xist.i0 - i0)
ylabel('\delta i0, rad');
subplot(2,2,4); plot(tmod, Xist.i_dot - i_dot)
ylabel('\delta i dot, rad/s');