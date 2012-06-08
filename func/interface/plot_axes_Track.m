%> ======================================================================
%> @brief plot Tracking graphs
%> @param handles Main handles struct
%> @param hF Flag (handle) for separate window
%> @param ij Indexies of axes (2-el string)
% ======================================================================
function plot_axes_Track(handles, hF, ij)
globals;

if hF 
    figure(hF); 
    hA = gca;
    XLab = 't, s';
else
    hA = eval(['handles.axes_Track_' ij]);
    set(handles.fig_main,'CurrentAxes', hA)
    XLab = '';
    YLab = '';
end

titl = '';
X1 = tmod; X2 = tmod; X3 = tmod; X4 = tmod;
Y1 = nan(1, Nmod); Y2 = nan(1, Nmod); Y3 = nan(1, Nmod); Y4 = nan(1, Nmod);
switch ij
    case '11'
        Y1 = Xest.e;
        Y2 = Xest_won.e;
        YLab = 'e';
        titl = 'Eccentricity';
    case '12'
        Y1 = Xest.d_e;
        Y2 = Xest_won.d_e;
        YLab = 'e';
        titl = 'Eccentricity rate';
    case '13'
        Y1 = Xest.p;
        Y2 = Xest_won.p;
        YLab = 'p, m';
        titl = 'Focal parameter';
    case '14'
        Y1 = Xest.d_p;
        Y2 = Xest_won.d_p;
        YLab = 'p, m';
        titl = 'Focal parameter rate';
    case '21'
        Y1 = Xest.theta;
        Y2 = Xest_won.theta;
        YLab = '\theta, rad';
        titl = 'True anomaly';
    case '22'
        Y1 = Xest.d_theta;
        Y2 = Xest_won.d_theta;
        YLab = '\theta'', rad/s';
        titl = 'True anomaly rate';
    case '23'
        Y1 = Xest.omega;
        Y2 = Xest_won.omega;
        YLab = '\omega, rad';
        titl = 'Argument of periapsis';
    case '24'
        Y1 = Xest.d_omega;
        Y2 = Xest_won.d_omega;
        YLab = '\omega'', rad/s';
        titl = 'Argument of periapsis rate';
    case '31'
        Y1 = Xest.Omega;
        Y2 = Xest_won.Omega;
        YLab = '\Omega, rad';
        titl = 'Longitude of the ascending node';
    case '32'
        Y1 = Xest.d_Omega;
        Y2 = Xest_won.d_Omega;
        YLab = '\Omega'', rad/s';
        titl = 'Longitude of the ascending node rate';
    case '33'
        Y1 = Xest.i;
        Y2 = Xest_won.i;
        YLab = 'i, rad';
        titl = 'Inclination';
    case '34'
        Y1 = Xest.d_i;
        Y2 = Xest_won.d_i;
        YLab = 'i'', rad/s';
        titl = 'Inclination rate';        

        
end

if isnan(Y4(1))
    if isnan(Y3(1))
        if isnan(Y2(1))
            if isnan(Y1(1))
                set(hA, 'Visible', 'off');
            else
                plot(hA, X1, Y1);
                set(hA, 'Visible', 'on');
            end
        else
            plot(hA, X1, Y1, X2, Y2);
            set(hA, 'Visible', 'on');
        end
    else
        plot(hA, X1, Y1, X2, Y2, X3, Y3);
        set(hA, 'Visible', 'on');
    end
else
    plot(hA, X1, Y1, X2, Y2, X3, Y3, X4, Y4);
    set(hA, 'Visible', 'on');
end

if hF
    title(hA, titl);
end

ly = ylabel(YLab);
lx = xlabel(XLab);
grid(hA, 'on');

if ~hF
    set(hA, 'XTick', []);
    set(hA, 'YTick', []);
    set(ly, 'FontSize', Font_Size);
end
