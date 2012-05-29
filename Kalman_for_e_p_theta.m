
F = [1 dTmod 0 0        0   0;       
     0 1     0 0        0   0;       
     0 0     1 dTmod    0   0;       
     0 0     0 1        0   0;       
     0 0     0 0        1   dTmod;   
     0 0     0 0        0   1];     

XXextr = [Xist.e(1)*0; Xist.d_e(1)*0; Xist.p(1)*0 + 20e6; Xist.d_p(1)*0; Xist.theta(1)*0.0; Xist.d_theta(1)*0];
Xs =[0; 20e6; 0];

std_e = 1e-6 / 15*dTmod;
std_p = 1e-5 / 15*dTmod;
std_theta = 1e-3 / 15*dTmod;

Dest = [std_e^2*1e4     0           0           0            0           0;
            0           std_e^2*1e4 0           0            0           0;
            0           0           std_p^2*1e4     0            0           0;
            0           0           0           std_p^2*1e4      0           0;
            0           0           0           0            std_theta^2*1e4 0;
            0           0           0           0            0           std_theta^2*1e4];

G = [0; std_e; 0; std_p; 0; std_theta];

std_Vr = 0.01;
std_Vu = 0.01;
std_r = 5;
Dn = [std_Vr^2 0            0;
      0        std_Vu^2     0
      0        0            std_r^2];
 
c = [1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0];
 
 
for i = 1:Nmod
    
    Vr_izm = Xest.d_r(i);
    Vu_izm = Xest.r(i)*Xest.d_u(i);
    r_izm = Xest.r(i);
%     XXextr(3) = Xist.p(i);

    e = XXextr(1);
    p = XXextr(3);
    theta = XXextr(5);
    munapi = sqrt(mu_earth / p);
    
    
    S0 = [munapi*sin(theta)                          munapi*cos(theta)                                  -p/(1+e*cos(theta))^2*cos(theta);
          -0.5*sqrt(mu_earth)*e*sin(theta)*p^(-3/2)  -0.5*sqrt(mu_earth)*(1+e*cos(theta))*p^(-3/2)      1/(1+e*cos(theta));
          munapi*e*cos(theta)                        -munapi*e*sin(theta)                               p/(1+e*cos(theta))^2*e*sin(theta)]';
    S = (c'*S0')';
        
    Vr_extr = munapi*e*sin(theta); 
    Vu_extr = munapi*(1+e*cos(theta));
    r_extr = p / (1+e*cos(theta));
    
    dY = [Vr_izm; Vu_izm; r_izm] - [Vr_extr; Vu_extr; r_extr];
    
    Dextr = F*Dest*F' + G*G';
%     t1 = S'/Dn*S;
%     t2 = inv(Dextr);
%     Dest = inv(t1 + t2); 

   
    XXest = XXextr + Dest*S'/Dn*dY;
    XXextr = F*XXest;
    
    Xest2.e(i) = XXest(1);
    Xest2.p(i) = XXest(3);
    Xest2.theta(i) = XXest(5);
    Xest2.omega(i) = Xest.u(i) - Xest2.theta(i);

    Xs = fsolve(@(Xfs)(fsolve_e_p_theta(Xfs, Vr_izm, Vu_izm, r_izm)), Xs, options_solve);
    Xest3.e(i) = Xs(1);
    Xest3.p(i) = Xs(2);
    Xest3.theta(i) = Xs(3);

    Xest3.theta(i) = (Xest3.e(i)<0).*(Xest3.theta(i) + pi) + (Xest3.e(i)>=0).*Xest3.theta(i);
    Xest3.omega(i) = Xest.u(i) - Xest3.theta(i);
    Xest3.e(i) = (Xest3.e(i)<0).*Xest3.e(i)*(-1) + (Xest3.e(i)>=0).*Xest3.e(i);    
    
    [Xest2.x0(i) Xest2.y0(i) Xest2.z0(i) Xest2.Vx(i) Xest2.Vy(i) Xest2.Vz(i)] = ...
        get_vector_XV( Xest2.e(i), Xest2.p(i), Xest2.theta(i), Xest2.omega(i), Xest.Omega(i), Xest.i(i));

    [Xest3.x0(i) Xest3.y0(i) Xest3.z0(i) Xest3.Vx(i) Xest3.Vy(i) Xest3.Vz(i)] = ...
        get_vector_XV( Xest3.e(i), Xest3.p(i), Xest3.theta(i), Xest3.omega(i), Xest.Omega(i), Xest.i(i));
    
    if ~mod(i, fix(Nmod/10))
        fprintf('Done %.0f %% \n', i/Nmod*100);
    end
end 
 

hF = 0;
hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xist.e, tmod, Xest.e, tmod, Xest2.e, tmod, Xest3.e)
ylabel('e');
subplot(2,1,2); plot(tmod, Xest.e - Xist.e, tmod, Xest2.e - Xist.e)
ylabel('d e');


hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xist.p, tmod, Xest.p, tmod, Xest2.p, tmod, Xest3.p)
ylabel('p');
subplot(2,1,2); plot(tmod, Xest.p - Xist.p, tmod, Xest2.p - Xist.p)
ylabel('d p');

hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xist.theta, tmod, Xest.theta, tmod, Xest2.theta, tmod, Xest3.theta)
ylabel('\theta');
subplot(2,1,2); plot(tmod, Xest.theta - Xist.theta, tmod, Xest2.theta - Xist.theta)
ylabel('d theta');

hF = figure(hF+1);
subplot(3,1,1); 
plot(tmod, sqrt((Xest.x0 - Xist.x0).^2 + (Xest.y0 - Xist.y0).^2 + (Xest.z0 - Xist.z0).^2))
ylabel('\delta xyz 1, m');
subplot(3,1,2); 
plot(tmod, sqrt((Xest2.x0 - Xist.x0).^2 + (Xest2.y0 - Xist.y0).^2 + (Xest2.z0 - Xist.z0).^2))
ylabel('\delta xyz 2, m');
subplot(3,1,3); 
plot(tmod, sqrt((Xest3.x0 - Xist.x0).^2 + (Xest3.y0 - Xist.y0).^2 + (Xest3.z0 - Xist.z0).^2))
ylabel('\delta xyz 3, m');
