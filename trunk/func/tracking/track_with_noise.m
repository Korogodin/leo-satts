%> ======================================================================
%> @brief Tracking alghorithm
%> @param handles Main handles struct
%> ======================================================================
function track_with_noise( handles, std_x, std_V )

set(handles.pb_Track, 'Enable', 'off'); 
pause(0.01);

globals;

% Filter step
T = dTmod;

% Evolutionary matrix
F = [1   T   0   0   0   0   0   0   0   0   0   0;
     0   1   0   0   0   0   0   0   0   0   0   0;
     0   0   1   T   0   0   0   0   0   0   0   0;
     0   0   0   1   0   0   0   0   0   0   0   0;
     0   0   0   0   1   T   0   0   0   0   0   0;
     0   0   0   0   0   1   0   0   0   0   0   0;
     0   0   0   0   0   0   1   T   0   0   0   0;
     0   0   0   0   0   0   0   1   0   0   0   0;
     0   0   0   0   0   0   0   0   1   T   0   0;
     0   0   0   0   0   0   0   0   0   1   0   0;
     0   0   0   0   0   0   0   0   0   0   1   T;
     0   0   0   0   0   0   0   0   0   0   0   1];

p_mult = 5e7; % To simplify matrix calculations - reducing the dynamic range

% Xextr.X =  [Xest_won.e(1); Xest_won.p(1)/p_mult; Xest_won.theta(1); 0; Xest_won.omega(1); 0; Xest_won.Omega(1); 0; Xest_won.i(1); 0];
% Xextr.X =  [0.01; 0; Xest_won.p(1)/p_mult; 0; 1; Xest_won.d_theta(1)*1.1; 0; 0; 0; Xest_won.d_Omega(1)*0.9; Xest_won.i(1)*0.9; Xest_won.d_i(1)];
Xextr.X =  [0.01; 0; 7e6/p_mult; 0; 1; Xest_won.d_theta(1)*1.1; 0; 0; 0; Xest_won.d_Omega(1)*0.9; Xest_won.i(1)*0.9; Xest_won.d_i(1)];
Xest.X = Xextr.X;

% RMS of shaping noises
std_e = 5e-6 / 15*dTmod;
std_p = 1 / 15*dTmod / p_mult; % [m]
std_theta = 1e-5 / 15*dTmod; % [rad]
std_omega = 1e-6 / 15*dTmod; % [rad]
std_Omega = 1e-8 / 15*dTmod; % [rad]
std_i = 1e-8 / 15*dTmod; % [rad]

Dest = [std_e^2*1e1     0           0           0               0               0               0               0               0               0               0           0
            0           std_e^2*1e2 0           0               0               0               0               0               0               0               0           0
            0           0           std_p^2*1e3 0               0               0               0               0               0               0               0           0
            0           0           0           std_p^2*1e4     0               0               0               0               0               0               0           0
            0           0           0           0               std_theta^2*1e2 0               0               0               0               0               0           0
            0           0           0           0               0               std_theta^2*1e0 0               0               0               0               0           0
            0           0           0           0               0               0               std_omega^2*1e2 0               0               0               0           0
            0           0           0           0               0               0               0               std_omega^2*1e2 0               0               0           0
            0           0           0           0               0               0               0               0               std_Omega^2*1e2 0               0           0
            0           0           0           0               0               0               0               0               0               std_Omega^2*1e2 0           0
            0           0           0           0               0               0               0               0               0               0               std_i^2*1e2 0
            0           0           0           0               0               0               0               0               0               0               0           std_i^2*1e2];

G =  [0 0 0 0 0 0;
      1 0 0 0 0 0;
      0 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 0 0 0 0;
      0 0 1 0 0 0;
      0 0 0 0 0 0;
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

% std_x = 10 / sqrt(dTmod)  ;
% std_V = 0.01 / sqrt(dTmod);

Dn = [std_x^2  0         0          0       0       0
      0        std_x^2   0          0       0       0
      0        0         std_x^2    0       0       0
      0        0         0          std_V^2 0       0
      0        0         0          0       std_V^2 0
      0        0         0          0       0       std_V^2];
 
c = [1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0 0 0; 
     0 0 0 0 0 0 0 0 0 0 1 0];
 
for i = 1:Nmod
    
    x0_izm = Xist.x0(i) + std_x * randn(1,1);
    y0_izm = Xist.y0(i) + std_x * randn(1,1);
    z0_izm = Xist.z0(i) + std_x * randn(1,1);
    Vx_izm = Xist.d_x0(i) + std_V * randn(1,1);    
    Vy_izm = Xist.d_y0(i) + std_V * randn(1,1);
    Vz_izm = Xist.d_z0(i) + std_V * randn(1,1);

    % For testing and debug    
%     Xextr.X(1) = Xest_won.e(i);
%     Xextr.X(3) = Xest_won.p(i)/p_mult;
%     Xextr.X(5) = Xest_won.theta(i);
%     Xextr.X(7) = Xest_won.omega(i);
%     Xextr.X(9) = Xest_won.Omega(i);
%     Xextr.X(11) = Xest_won.i(i);

    e = Xextr.X(1);
    p = Xextr.X(3) * p_mult;
    theta = Xextr.X(5);
    omega = Xextr.X(7);
    Omega = Xextr.X(9);
    i0 = Xextr.X(11);

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
        
    S0 = nan(6,6);

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
    S0(:,2) = S0(:,2) * p_mult;
    
    % Discriminator for theta  
    Sc_x_theta = - sin(u)*cos(Omega) - cos(u)*sin(Omega)*cos(i0);
    Sc_y_theta = - sin(u)*sin(Omega) + cos(u)*cos(Omega)*cos(i0);
    Sc_z_theta = cos(u)*sin(i0);
    Sc_Vx_theta = cos(u)*cos(Omega) - sin(u)*sin(Omega)*cos(i0);
    Sc_Vy_theta = cos(u)*sin(Omega) + sin(u)*cos(Omega)*cos(i0);
    Sc_Vz_theta = -sin(u)*sin(i0);
    
    S0(1, 3) = e*Sc_x*p / oec^2 * sin(theta) + p / oec * Sc_x_theta;
    S0(2, 3) = e*Sc_y*p / oec^2 * sin(theta) + p / oec * Sc_y_theta;
    S0(3, 3) = e*Sc_z*p / oec^2 * sin(theta) + p / oec * Sc_z_theta;
    S0(4, 3) = -munapi*Sc_Vx_theta;
    S0(5, 3) = -munapi*Sc_Vy_theta;
    S0(6, 3) = munapi*Sc_Vz_theta;

%     dx = 0.001;
%     [x0_extr2 y0_extr2 z0_extr2 Vx_extr2 Vy_extr2 Vz_extr2] = ...
%         get_vector_XV( e, p, theta+dx, omega, Omega, i0);
%     dS = [x0_extr2; y0_extr2; z0_extr2; Vx_extr2; Vy_extr2; Vz_extr2] - ...
%             [x0_extr; y0_extr; z0_extr; Vx_extr; Vy_extr; Vz_extr];
%         
%     S0(1, 3) = dS(1) / dx;
%     S0(2, 3) = dS(2) / dx;
%     S0(3, 3) = dS(3) / dx;
%     S0(4, 3) = dS(4) / dx;
%     S0(5, 3) = dS(5) / dx;
%     S0(6, 3) = dS(6) / dx;
    
    % Discriminator for omega
    Sc_x_omega = Sc_x_theta;
    Sc_y_omega = Sc_y_theta;
    Sc_z_omega = Sc_z_theta;
    Sc_Vx_omega = Sc_Vx_theta;
    Sc_Vy_omega = Sc_Vy_theta;
    Sc_Vz_omega = Sc_Vz_theta;

    S0(1,4) = p / oec * Sc_x_omega;
    S0(2,4) = p / oec * Sc_y_omega;
    S0(3,4) = p / oec * Sc_z_omega;
    S0(4,4) = munapi * (e*sin(omega)*Sc_x_omega - oec*Sc_Vx_omega);
    S0(5,4) = munapi * (e*sin(omega)*Sc_y_omega - oec*Sc_Vy_omega);
    S0(6,4) = munapi * (e*sin(omega)*Sc_z_omega + oec*Sc_Vz_omega);

    S01(i) = S0(1,4);
    S02(i) = S0(2,4);
    S03(i) = S0(3,4);
    S04(i) = S0(4,4);
    S05(i) = S0(5,4);
    S06(i) = S0(6,4);
    
    dx = 1e-6;
    [x0_extr2 y0_extr2 z0_extr2 Vx_extr2 Vy_extr2 Vz_extr2] = ...
        get_vector_XV( e, p, theta, omega+dx, Omega, i0);
    dS = [x0_extr2; y0_extr2; z0_extr2; Vx_extr2; Vy_extr2; Vz_extr2] - ...
            [x0_extr; y0_extr; z0_extr; Vx_extr; Vy_extr; Vz_extr];
        
    S0(1, 4) = dS(1) / dx;
    S0(2, 4) = dS(2) / dx;
    S0(3, 4) = dS(3) / dx;
    S0(4, 4) = dS(4) / dx;
    S0(5, 4) = dS(5) / dx;
    S0(6, 4) = dS(6) / dx;    
    
    S01_ist(i) = S0(1,4);
    S02_ist(i) = S0(2,4);
    S03_ist(i) = S0(3,4);
    S04_ist(i) = S0(4,4);
    S05_ist(i) = S0(5,4);
    S06_ist(i) = S0(6,4); 
 
    
    % Discriminator for Omega
    Sc_x_Omega = -cos(u)*sin(Omega) - sin(u)*cos(Omega)*cos(i0);
    Sc_y_Omega = cos(u)*cos(Omega) - sin(u)*sin(Omega)*cos(i0);
    Sc_z_Omega = 0;
    Sc_Vx_Omega = -sin(u)*sin(Omega) + cos(u)*cos(Omega)*cos(i0);
    Sc_Vy_Omega = sin(u)*cos(Omega) + cos(u)*sin(Omega)*cos(i0);
    Sc_Vz_Omega = 0;
    
    S0(1,5) = p/oec*Sc_x_Omega;
    S0(2,5) = p/oec*Sc_y_Omega;
    S0(3,5) = 0;
    S0(4,5) = munapi * (e*sin(omega)*Sc_x_Omega - oec*Sc_Vx_Omega);
    S0(5,5) = munapi * (e*sin(omega)*Sc_y_Omega - oec*Sc_Vy_Omega);
    S0(6,5) = 0;
    
%     dx = 0.001;
%     [x0_extr2 y0_extr2 z0_extr2 Vx_extr2 Vy_extr2 Vz_extr2] = ...
%         get_vector_XV( e, p, theta, omega, Omega+dx, i0);
%     dS = [x0_extr2; y0_extr2; z0_extr2; Vx_extr2; Vy_extr2; Vz_extr2] - ...
%             [x0_extr; y0_extr; z0_extr; Vx_extr; Vy_extr; Vz_extr];
%         
%     S0(1, 5) = dS(1) / dx;
%     S0(2, 5) = dS(2) / dx;
%     S0(3, 5) = dS(3) / dx;
%     S0(4, 5) = dS(4) / dx;
%     S0(5, 5) = dS(5) / dx;
%     S0(6, 5) = dS(6) / dx;  

    % Discriminator for i
    Sc_x_i = sin(u)*sin(Omega)*sin(i0);
    Sc_y_i = -sin(u)*cos(Omega)*sin(i0);
    Sc_z_i = sin(u)*cos(i0);
    Sc_Vx_i = -cos(u)*sin(Omega)*sin(i0);
    Sc_Vy_i = cos(u)*cos(Omega)*sin(i0);
    Sc_Vz_i = cos(u)*cos(i0);
    
    S0(1, 6) = p/oec * Sc_x_i;
    S0(2, 6) = p/oec * Sc_y_i;
    S0(3, 6) = p/oec * Sc_z_i;
    S0(4, 6) = munapi * (e*sin(theta)*Sc_x_i - oec*Sc_Vx_i);
    S0(5, 6) = munapi * (e*sin(theta)*Sc_y_i - oec*Sc_Vy_i);
    S0(6, 6) = munapi * (e*sin(theta)*Sc_z_i + oec*Sc_Vz_i);
    
%     dx = 0.001;
%     [x0_extr2 y0_extr2 z0_extr2 Vx_extr2 Vy_extr2 Vz_extr2] = ...
%         get_vector_XV( e, p, theta, omega, Omega, i0+dx);
%     dS = [x0_extr2; y0_extr2; z0_extr2; Vx_extr2; Vy_extr2; Vz_extr2] - ...
%             [x0_extr; y0_extr; z0_extr; Vx_extr; Vy_extr; Vz_extr];
%         
%     S0(1, 6) = dS(1) / dx;
%     S0(2, 6) = dS(2) / dx;
%     S0(3, 6) = dS(3) / dx;
%     S0(4, 6) = dS(4) / dx;
%     S0(5, 6) = dS(5) / dx;
%     S0(6, 6) = dS(6) / dx;  
    
    S = (c'*S0')';
   
    Dextr = F*Dest*F' + GDgG;
    t1 = S'/Dn*S;
    t2 = inv(Dextr);
    Dest = inv(t1 + t2); 
   
    Xest.X = Xextr.X + Dest*S'/Dn*dY;
    Xextr.X = F*Xest.X;
    
    if Xest.X(1) > 0;
        Xest.e(i) = Xest.X(1);
        Xest.theta(i) = mod_pm_pi(Xest.X(5));    
        Xest.omega(i) = mod_pm_pi(Xest.X(7));
    else
        Xest.e(i) = -Xest.X(1);
        Xest.theta(i) = mod_pm_pi(-Xest.X(5));    
        Xest.omega(i) = mod_pm_pi(-Xest.X(7));
    end    
    Xest.p(i) = Xest.X(3)*p_mult;
    Xest.Omega(i) = mod_pm_pi(Xest.X(9));
    Xest.i(i) = mod_pm_pi(Xest.X(11));
    
    [Xest.x0(i) Xest.y0(i) Xest.z0(i) Xest.d_x0(i) Xest.d_y0(i) Xest.d_z0(i)] = ...
        get_vector_XV( Xest.e(i), Xest.p(i), Xest.theta(i), Xest.omega(i), Xest.Omega(i), Xest.i(i) );
    [Xest.x(i) Xest.y(i) Xest.z(i) tmp1 tmp2 tmp3] = ...
        get_vector_XV( Xest.e(i), Xest.p(i), Xest.theta(i), Xest.omega(i), Xest.Omega(i) - omega_e*tmod(i), Xest.i(i) );
    
    Xest.d_e(i) = Xest.X(2);
    Xest.d_p(i) = Xest.X(4)*p_mult;
    Xest.d_theta(i) = Xest.X(6);
    Xest.d_omega(i) = Xest.X(8);
    Xest.d_Omega(i) = Xest.X(10);
    Xest.d_i(i) = Xest.X(12);
    
    if ~mod(i, fix(Nmod/100))
        set(handles.txt_Track, 'String', sprintf('%.0f %%', i/Nmod*100));
    end
    drawnow
    pause(0.01);    
end 
set(handles.txt_Track, 'String', sprintf('%.0f %%', 100));

set(handles.pb_Track, 'Enable', 'on'); 

set(handles.pan_Orbital_Elements, 'Visible', 'off');
set(handles.pan_Tracking, 'Visible', 'on');

% Output results to form
for i = 1:5
    for j = 1:5
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'DrawMode', 'fast');
        plot_axes_Track(handles, 0, [num2str(i) num2str(j)]);
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'Tag', ['axes_Track_' num2str(i) num2str(j)]);
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'ButtonDownFcn', str2func('@(hObject,eventdata)fig_main(''axes_Track_ButtonDownFcn'',hObject,eventdata,guidata(hObject))'));
    end
end
plot_axes_Earth(handles, 0);
draw_Errors(handles);

end

