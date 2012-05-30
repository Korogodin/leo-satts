
T = dTmod;
F = [1  0   0   0   0   0   0   0   0;
     0  1   0   0   0   0   0   0   0;
     0  0   1   T   0   0   0   0   0;
     0  0   0   1   0   0   0   0   0;
     0  0   0   0   1   0   0   0   0;
     0  0   0   0   0   1   T   0   0;
     0  0   0   0   0   0   1   0   0;
     0  0   0   0   0   0   0   1   T;
     0  0   0   0   0   0   0   0   1];

Xextr4 = Xest;
Xest4 = Xest;
 
% Xextr.X = [e; p; theta; theta'; omega; Omega; Omega'; i; i'];
Xextr4.X = [0; 25e6; Xest.theta(1); 0; Xest.omega(1); Xest.Omega(1); 0; Xest.i(1); 0];
Xest4.X = Xextr4.X;
% Xs =[0; 20e6; 0];

std_e = 1e-6 / 15*dTmod;
std_p = 100 / 15*dTmod;
std_theta = 1e-3 / 15*dTmod;
std_omega = 1e-3 / 15*dTmod;
std_Omega = 1e-3 / 15*dTmod;
std_i = 1e-3 / 15*dTmod;

Dest = [std_e^2*1e1     0           0                   0               0               0               0               0           0
            0           std_p^2*1e2 0                   0               0               0               0               0           0
            0           0           std_theta^2/1e9     0               0               0               0               0           0
            0           0           0                   std_theta^2/1e9 0               0               0               0           0
            0           0           0                   0               std_omega^2/1e9 0               0               0           0
            0           0           0                   0               0               std_Omega^2/1e9 0               0           0
            0           0           0                   0               0               0               std_Omega^2/1e9 0           0
            0           0           0                   0               0               0               0               std_i^2/1e9 0
            0           0           0                   0               0               0               0               0           std_i^2/1e9];

G =  [1 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 0 0 0 0;
      0 0 1 0 0 0;
      0 0 0 1 0 0;
      0 0 0 0 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 0;
      0 0 0 0 0 1];

Dg = [std_e^2 0       0           0           0           0
      0       std_p^2 0           0           0           0 
      0       0       std_theta^2 0           0           0
      0       0       0           std_omega^2 0           0
      0       0       0           0           std_Omega^2 0
      0       0       0           0           0           std_i^2];

GDgG = G*Dg*G';                        

std_x = 10 / sqrt(dTmod) /10;
std_V = 0.01 / sqrt(dTmod);

Dn = [std_x^2  0         0          0       0       0
      0        std_x^2   0          0       0       0
      0        0         std_x^2    0       0       0
      0        0         0          std_V^2 0       0
      0        0         0          0       std_V^2 0
      0        0         0          0       0       std_V^2];
 
c = [1 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 1 0 0 0; 
     0 0 0 0 0 0 0 1 0];
 
for i = 1:Nmod
    
    x0_izm = Xist.x0(i) + std_x * randn(1,1);
    y0_izm = Xist.y0(i) + std_x * randn(1,1);
    z0_izm = Xist.z0(i) + std_x * randn(1,1);
    Vx_izm = Xist.d_x0(i) + std_V * randn(1,1);    
    Vy_izm = Xist.d_y0(i) + std_V * randn(1,1);
    Vz_izm = Xist.d_z0(i) + std_V * randn(1,1);

    
%     Xextr4.X(1) = Xest.e(i);
%     Xextr4.X(2) = Xest.p(i);
    Xextr4.X(3) = Xest.theta(i);
    Xextr4.X(5) = Xest.omega(i);
    Xextr4.X(6) = Xest.Omega(i);
    Xextr4.X(8) = Xest.i(i);

    e = Xextr4.X(1);
    p = Xextr4.X(2);
    theta = Xextr4.X(3);
    omega = Xextr4.X(5);
    Omega = Xextr4.X(6);
    i0 = Xextr4.X(8);
    munapi = sqrt(mu_earth / p);
    u = theta + omega;
    
    oec = 1 + e*cos(theta);
    es = e*sin(theta);
    
    Sc_x = cos(u)*cos(Omega) - sin(u)*sin(Omega)*cos(i0);
    Sc_y = cos(u)*sin(Omega) + sin(u)*cos(Omega)*cos(i0);
    Sc_z = sin(u)*sin(i0);
    
    Sc_Vx = sin(u)*cos(Omega) + cos(u)*sin(Omega)*cos(i0);
    Sc_Vy = sin(u)*sin(Omega) - cos(u)*cos(Omega)*cos(i0);
    Sc_Vz = cos(u)*sin(i0);

    [x0_extr y0_extr z0_extr Vx_extr Vy_extr Vz_extr] = ...
        get_vector_XV( e, p, theta, omega, Omega, i0);
    
    dY = [x0_izm; y0_izm; z0_izm; Vx_izm; Vy_izm; Vz_izm] - ...
            [x0_extr; y0_extr; z0_extr; Vx_extr; Vy_extr; Vz_extr];
    dY_ist = [Xist.x0(i); Xist.y0(i); Xist.z0(i); Xist.d_x0(i); Xist.d_y0(i); Xist.d_z0(i)] - ...
            [x0_extr; y0_extr; z0_extr; Vx_extr; Vy_extr; Vz_extr];
        
    S0 = [NaN   NaN   1e-19     1e-19     1e-19     1e-19;       
          NaN   NaN   1e-19     1e-19     1e-19     1e-19; 
          NaN   NaN   1e-19     1e-19     1e-19     1e-19;
          NaN   NaN   1e-19     1e-19     1e-19     1e-19;   
          NaN   NaN   1e-19     1e-19     1e-19     1e-19;
          NaN   NaN   1e-19     1e-19     1e-19     1e-19];
      
    % Discriminator for e
    S0(1, 1) = - p * Sc_x / (1+e*cos(theta))^2 * cos(theta);  
    S0(2, 1) = - p * Sc_y / (1+e*cos(theta))^2 * cos(theta);
    S0(3, 1) = - p * Sc_z / (1+e*cos(theta))^2 * cos(theta);
    S0(4, 1) = munapi * (sin(theta)*Sc_x - cos(theta)*Sc_Vx);
    S0(5, 1) = munapi * (sin(theta)*Sc_y - cos(theta)*Sc_Vy);
    S0(6, 1) = munapi * (sin(theta)*Sc_z + cos(theta)*Sc_Vz);
    
    % Discriminator for p  
    S0(1, 2) = Sc_x/oec;
    S0(2, 2) = Sc_y/oec;
    S0(3, 2) = Sc_z/oec;
    S0(4, 2) = - 0.5 * 1/p * munapi * (es *Sc_x - oec*Sc_Vx);
    S0(5, 2) = - 0.5 * 1/p * munapi * (es *Sc_y - oec*Sc_Vy);
    S0(6, 2) = - 0.5 * 1/p * munapi * (es *Sc_z + oec*Sc_Vz);

    % Discriminator for theta  
    Sc_x_theta = - sin(u)*cos(Omega) - cos(u)*sin(Omega)*cos(i0);
    Sc_y_theta = - sin(u)*sin(Omega) + cos(u)*cos(Omega)*cos(i0);
    Sc_z_theta = cos(u)*sin(i0);
    
    S0(1, 3) = Sc_x*p / oec^2 * sin(theta) + p / oec * Sc_x_theta;
    S0(2, 3) = Sc_y*p / oec^2 * sin(theta) + p / oec * Sc_y_theta;
    S0(3, 3) = Sc_z*p / oec^2 * sin(theta) + p / oec * Sc_z_theta;
    S0(4, 3) = 0;
    S0(5, 3) = 0;
    S0(6, 3) = 0;
    
    S01(i) = S0(1,3);
    S02(i) = S0(2,3);
    S03(i) = S0(3,3);
    S04(i) = S0(4,3);
    S05(i) = S0(5,3);
    S06(i) = S0(6,3);
    
%     S0(1, 3) = dY_ist(1) / (Xest.theta(i) - Xextr.theta(i));
%     S0(2, 3) = dY_ist(2) / (Xest.theta(i) - Xextr.theta(i));
%     S0(3, 3) = dY_ist(3) / (Xest.theta(i) - Xextr.theta(i));
%     S0(4, 3) = dY_ist(4) / (Xest.theta(i) - Xextr.theta(i));
%     S0(5, 3) = dY_ist(5) / (Xest.theta(i) - Xextr.theta(i));
%     S0(6, 3) = dY_ist(6) / (Xest.theta(i) - Xextr.theta(i));

    S01_ist(i) = S0(1,3);
    S02_ist(i) = S0(2,3);
    S03_ist(i) = S0(3,3);
    S04_ist(i) = S0(4,3);
    S05_ist(i) = S0(5,3);
    S06_ist(i) = S0(6,3);

    S = (c'*S0')';
    
    Dextr = F*Dest*F' + GDgG;
    t1 = S'/Dn*S;
    t2 = inv(Dextr);
    if sum(sum(isnan(t1 + t2)))
        a 
    end
    Dest = inv(t1 + t2); 

    
    if sum(sum(isnan(Dest)))
        a 
    end
    
    Xest4.X = Xextr4.X + Dest*S'/Dn*dY;
    
    if sum(sum(isnan(Xest4.X)))
        a 
    end
    Xextr4.X = F*Xest4.X;
    
    Xest4.e(i) = Xest4.X(1);
    Xest4.p(i) = Xest4.X(2);
    Xest4.theta(i) = Xest4.X(3);    
%     Xest4.omega(i) = Xest4.X(5);
%     Xest4.Omega(i) = Xest4.X(6);
%     Xest4.i(i) = Xest4.X(8);

    [Xest4.x0(i) Xest4.y0(i) Xest4.z0(i) Xest4.Vx(i) Xest4.Vy(i) Xest4.Vz(i)] = ...
        get_vector_XV( Xest4.e(i), Xest4.p(i), Xest4.theta(i), Xest4.omega(i), Xest4.Omega(i), Xest4.i(i));
    
    if ~mod(i, fix(Nmod/10))
        fprintf('Done %.0f %% \n', i/Nmod*100);
    end
end 
 

hF = 0;
hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xest.e, tmod, Xest4.e, tmod, Xist.e)
ylabel('e');
subplot(2,1,2); plot(tmod, Xest4.e - Xest.e)
ylabel('d e');


hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xest.p, tmod, Xest4.p, tmod, Xist.p)
ylabel('p');
subplot(2,1,2); plot(tmod, Xest4.p - Xest.p)
ylabel('d p');

hF = figure(hF+1);
subplot(2,1,1); plot(tmod, Xest4.theta, tmod, Xest.theta, tmod, Xist.theta)
ylabel('\theta');
subplot(2,1,2); plot(tmod, Xest4.theta - Xest.theta)
ylabel('d theta');

hF = figure(hF+1);
subplot(2,1,1); 
plot(tmod, sqrt((Xest.x0 - Xist.x0).^2 + (Xest.y0 - Xist.y0).^2 + (Xest.z0 - Xist.z0).^2))
ylabel('\delta xyz 1, m');
subplot(2,1,2); 
plot(tmod, sqrt((Xest4.x0 - Xist.x0).^2 + (Xest4.y0 - Xist.y0).^2 + (Xest4.z0 - Xist.z0).^2))
ylabel('\delta xyz 4, m');

hF = figure(hF+1);
subplot(3,2,1); plot(tmod, S01, tmod, S01_ist)
ylabel('1');
subplot(3,2,3); plot(tmod, S02, tmod, S02_ist)
ylabel('2');
subplot(3,2,5); plot(tmod, S03, tmod, S03_ist)
ylabel('3');
subplot(3,2,2); plot(tmod, S04, tmod, S04_ist)
ylabel('4');
subplot(3,2,4); plot(tmod, S05, tmod, S05_ist)
ylabel('5');
subplot(3,2,6); plot(tmod, S06, tmod, S06_ist)
ylabel('6');

