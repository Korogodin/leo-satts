
std_X = 10 / sqrt(dTmod) * 1;
std_V = 0.01 / sqrt(dTmod) * 1;
Xizm = Xist.x0 + randn(1, Nmod)*std_X;
Yizm = Xist.y0 + randn(1, Nmod)*std_X;
Zizm = Xist.z0 + randn(1, Nmod)*std_X;
VXizm = Xist.d_x0 + randn(1, Nmod)*std_V;
VYizm = Xist.d_y0 + randn(1, Nmod)*std_V;
VZizm = Xist.d_z0 + randn(1, Nmod)*std_V;

Xextr.theta(1) = Xist.theta(1);
Xextr.omega(1) = Xist.omega(1);
Xextr.Omega(1) = Xist.Omega(1);
Xextr.i(1) = Xist.i(1);
Xextr.e(1) = Xist.e(1);
Xextr.p(1) = Xist.p(1);
Xextr.d_theta(1) = Xist.d_theta(1);
Xextr.d_omega(1) = Xist.d_omega(1);
Xextr.d_Omega(1) = Xist.d_Omega(1);
Xextr.d_i(1) = Xist.d_i(1);
Xextr.d_e(1) = 0;
Xextr.d_p(1) = 0;

K_e = get_K2(0.0001, dTmod);
K_p = get_K2(0.0001, dTmod);
K_Omega = get_K2(0.025, dTmod);
K_omega = get_K2(0.0001, dTmod);
K_i = get_K2(0.02, dTmod);
K_theta = get_K2(0.02, dTmod);
K_e(1) = 1;
K_p(1) = 1;
K_Omega(1) = 1;
K_omega(1) = 1;
K_i(1) = 1;
K_theta(1) = 1;

for i = 1:Nmod
    
    Xs(1) = Xextr.theta(i);
    Xs(2) = Xextr.omega(i);
    Xs(3) = Xextr.Omega(i);
    Xs(4) = Xextr.i(i);
    Xs(5) = Xextr.e(i);
    Xs(6) = Xextr.p(i);
    
    Xs = fsolve(@(Xfs)(fsolve_Kepler(Xfs, Xizm(i), Yizm(i), Zizm(i),...
         VXizm(i), VYizm(i), VZizm(i))), Xs, options_solve);
    
    Xest.theta(i) = Xextr.theta(i) + K_theta(1)*(Xs(1) - Xextr.theta(i));
    Xest.d_theta(i) = Xextr.d_theta(i) + K_theta(2)*(Xs(1) - Xextr.theta(i));
    Xextr.theta(i+1) = Xest.theta(i) + Xest.d_theta(i)*dTmod;
    Xextr.d_theta(i+1) = Xextr.d_theta(i);

    Xest.omega(i) = Xextr.omega(i) + K_omega(1)*(Xs(2) - Xextr.omega(i));
    Xest.d_omega(i) = Xextr.d_omega(i) + K_omega(2)*(Xs(2) - Xextr.omega(i));
    Xextr.omega(i+1) = Xest.omega(i) + Xest.d_omega(i)*dTmod;
    Xextr.d_omega(i+1) = Xextr.d_omega(i);

    Xest.Omega(i) = Xextr.Omega(i) + K_Omega(1)*(Xs(3) - Xextr.Omega(i));
    Xest.d_Omega(i) = Xextr.d_Omega(i) + K_Omega(2)*(Xs(3) - Xextr.Omega(i));
    Xextr.Omega(i+1) = Xest.Omega(i) + Xest.d_Omega(i)*dTmod;
    Xextr.d_Omega(i+1) = Xextr.d_Omega(i);

    Xest.i(i) = Xextr.i(i) + K_i(1)*(Xs(4) - Xextr.i(i));
    Xest.d_i(i) = Xextr.d_i(i) + K_i(2)*(Xs(4) - Xextr.i(i));
    Xextr.i(i+1) = Xest.i(i) + Xest.d_i(i)*dTmod;
    Xextr.d_i(i+1) = Xextr.d_i(i);

    Xest.e(i) = Xextr.e(i) + K_e(1)*(Xs(5) - Xextr.e(i));
    Xest.d_e(i) = Xextr.d_e(i) + K_e(2)*(Xs(5) - Xextr.e(i));
    Xextr.e(i+1) = Xest.e(i) + Xest.d_e(i)*dTmod;
    Xextr.d_e(i+1) = Xextr.d_e(i);    
    
    Xest.p(i) = Xextr.p(i) + K_p(1)*(Xs(6) - Xextr.p(i));
    Xest.d_p(i) = Xextr.d_p(i) + K_p(2)*(Xs(6) - Xextr.p(i));
    Xextr.p(i+1) = Xest.p(i) + Xest.d_p(i)*dTmod;
    Xextr.d_p(i+1) = Xextr.d_p(i);
    
    Xest.r(i) = Xest.p(i) ./ (1 + Xest.e(i)*cos(Xest.theta(i)));
    [Xest.x0(i) Xest.y0(i) Xest.z0(i)] = get_vector_XYZ( Xest.r(i),...
                Xest.Omega(i), Xest.i(i), Xest.theta(i)+Xest.omega(i) );

    if ~mod(i, fix(Nmod/10))
        fprintf('Done %.0f %% \n', i/Nmod*100);
    end
end

% Xest.theta = (Xest.e<0).*(Xest.theta - pi) + (Xest.e>=0).*Xest.theta;
% Xest.omega = (Xest.e<0).*(Xest.omega - pi) + (Xest.e>=0).*Xest.omega;
% Xest.e = (Xest.e<0).*Xest.e*(-1) + (Xest.e>=0).*Xest.e;

hF = 0;
hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xist.i, tmod, Xest.i)
ylabel('i, rad');
subplot(2,1,2); plot(tmod, Xist.i - Xest.i)
ylabel('\delta i, rad');

hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xist.omega, tmod, Xest.omega)
ylabel('\omega, rad');
subplot(2,1,2); plot(tmod, Xist.omega - Xest.omega)
ylabel('\delta \omega, rad');

hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xist.Omega, tmod, Xest.Omega)
ylabel('\Omega, rad');
subplot(2,1,2); plot(tmod, Xist.Omega - Xest.Omega)
ylabel('\delta \Omega, rad');

hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xist.theta, tmod, Xest.theta)
ylabel('\theta, rad');
subplot(2,1,2); plot(tmod, Xist.theta - Xest.theta)
ylabel('\delta \theta, rad');

hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xist.e, tmod, Xest.e)
ylabel('e');
subplot(2,1,2); plot(tmod, Xist.e - Xest.e)
ylabel('\delta e');

hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xist.p, tmod, Xest.p)
ylabel('p, m');
subplot(2,1,2); plot(tmod, Xist.p - Xest.p)
ylabel('\delta p, m');

hF = figure(hF+1);
plot(tmod, sqrt((Xest.x0 - Xist.x0).^2 + (Xest.y0 - Xist.y0).^2 + (Xest.z0 - Xist.z0).^2))
ylabel('\delta xyz, m');
