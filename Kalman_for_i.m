
F = [1 dTmod 0 0;
     0 1     0 0;
     0 0     1 0;
     0 0     0 1];

XXextr = [Xist.i0(1); Xist.i_dot(1); Xist.Cic(1); Xist.Cis(1)];

std_dot = 1e-12 / 15*dTmod * 0.5 ;
std_ic = 5e-10 / 15*dTmod * 0.1;
std_is = 5e-10 / 15*dTmod * 0.1;
G = [0; std_dot; std_ic; std_is];

Dest = [std_dot^2   0           0           0;
        0           std_dot^2*100   0           0;
        0           0           std_ic^2    0;
        0           0           0           std_is^2];

std_izm = 1e-7;    
R = std_izm^2;    
for i = 1:Nmod
    
    H = [1 dTmod cos(2*Xist.u(i)) sin(2*Xist.u(i))];
    Dextr = F*Dest*F' + G*G';
    Dest = inv(inv(Dextr) + H'*inv(R)*H);
    
    K = Dest*H'*inv(R);
    
    XXest = XXextr + K*(Xist.i(i) + randn(1,1)*std_izm - H*XXextr);
    XXextr = F*XXest;
    
    i0(i) = XXest(1);
    i_dot(i) = XXest(2);
    Cic(i) = XXest(3);
    Cis(i) = XXest(4);
    
    [Xest.x0(i) Xest.y0(i) Xest.z0(i)] = ...
        get_vector_XYZ(Xist.r(i), Xist.lambda(i), i0(i), Xist.u(i));    
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

% Xest.Cic = Cic;
% Xest.Cis = Cis;
% Xest.i_dot = i_dot;
% Xest.i0 = i0;
hF = figure(hF+1);
plot(tmod, sqrt((Xest.x0 - Xist.x0).^2 + (Xest.y0 - Xist.y0).^2 + (Xest.z0 - Xist.z0).^2))
ylabel('\delta xyz, m');
