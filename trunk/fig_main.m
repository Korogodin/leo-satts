function varargout = fig_main(varargin)
% FIG_MAIN M-file for fig_main.fig
%      FIG_MAIN, by itself, creates a new FIG_MAIN or raises the existing
%      singleton*.
%
%      H = FIG_MAIN returns the handle to a new FIG_MAIN or the handle to
%      the existing singleton*.
%
%      FIG_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIG_MAIN.M with the given input arguments.
%
%      FIG_MAIN('Property','Value',...) creates a new FIG_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fig_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fig_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fig_main

% Last Modified by GUIDE v2.5 07-Jun-2012 11:47:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fig_main_OpeningFcn, ...
                   'gui_OutputFcn',  @fig_main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before fig_main is made visible.
function fig_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fig_main (see VARARGIN)

% Choose default command line output for fig_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.fig_main, 'Units', 'pixels');
set(0, 'Units', 'pixels');

% Locate form
ScreenSize = get(0,'ScreenSize');
if ((ScreenSize(3) < 1280)||(ScreenSize(4) < 720))
    msgbox('Sceeen size too small!');
    close(handles.fig_main);
end
FigWeigth = 1250; 
FigWeigth_border = fix((ScreenSize(3) - FigWeigth)/2);
FigHeigth = 670;
FigHeigth_border = fix((ScreenSize(4) - FigHeigth)/2);
set(handles.fig_main, 'Position', [FigWeigth_border FigHeigth_border FigWeigth FigHeigth]);
set(handles.pan_Tracking, 'Position', get(handles.pan_Orbital_Elements, 'Position'));

% Set widget's function
for i = 1:5
    for j = 1:5
        set(eval(['handles.axes_OE_' num2str(i) num2str(j)]), 'DrawMode', 'fast');
        plot_axes_OE(handles, 0, [num2str(i) num2str(j)]);
        set(eval(['handles.axes_OE_' num2str(i) num2str(j)]), 'Tag', ['axes_OE_' num2str(i) num2str(j)]);
        set(eval(['handles.axes_OE_' num2str(i) num2str(j)]), 'ButtonDownFcn', str2func('@(hObject,eventdata)fig_main(''axes_OE_ButtonDownFcn'',hObject,eventdata,guidata(hObject))'));
    end
end
globals;

% Set constants
start_handle = handles;


% UIWAIT makes fig_main wait for user response (see UIRESUME)
% uiwait(handles.fig_main);


% --- Outputs from this function are returned to the command line.
function varargout = fig_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_Run.
function pb_Run_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_axes_Earth(handles, 0);

   

%> ======================================================================
%> @brief 3D-Earth, SV's orbits
%> @param handles Main handles struct
%> @param hF Flag (handle) for separate window
% ======================================================================
function plot_axes_Earth(handles, hF)
globals;

R_earth_pol = 6356777; % Polar radius of Earth
R_earth_equa = 6378160; % Equatorial radius of Earth
Earth_axe_angl = deg2rad(23); 

[x,y,z] = sphere(50);
x = R_earth_equa * x; y = R_earth_equa * y; z = R_earth_pol * z;
load('topo.mat','topo','topomap1');
% cla reset
% axis square off
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;
alpha_earth = pi/1;
x_new = x*cos(alpha_earth) + y*sin(alpha_earth);
y_new = -x*sin(alpha_earth) + y*cos(alpha_earth);

if hF 
    figure(hF); 
    hA = gca;
else
    set(handles.fig_main,'CurrentAxes', handles.axes_3D)
    hA = handles.axes_3D;
end

cla(hA);
set(hA, 'FontSize', Font_Size);
surface(x_new, y_new, z, props);
light('position',[-1 0 1*tan(Earth_axe_angl)]*R_earth_equa*2);
view(3)        
hold(hA, 'on');
for i = 1:length(SV_GLO_List)
    plot3(hA, eval([SV_GLO_List{i} '.SVOrbit.X0']), ...
              eval([SV_GLO_List{i} '.SVOrbit.Y0']), ...
              eval([SV_GLO_List{i} '.SVOrbit.Z0']));   
end
hold(hA, 'off');
xlim([-3e7 3e7]); ylim([-3e7 3e7]); zlim([-3e7 3e7]);
title(hA, 'Space View');
    

% --- Executes on button press in pb_addSV.
function pb_addSV_Callback(hObject, eventdata, handles)
% hObject    handle to pb_addSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
globals;

p_orb = 19100e3 + 6500e3;
e_orb = 0.01;
i_orb = pi/3;
omega_p_orb = 0;

SV_GLO = SV(p_orb, e_orb, 0, 2*pi/3, omega_p_orb, i_orb, 'SV GLONASS 1', 'SV_GLO(1)');
SV_GLO(2) = SV(p_orb, e_orb, 0, pi/3, omega_p_orb, i_orb, 'SV GLONASS 2', 'SV_GLO(2)');



% --- Executes on mouse press over axes background.
function axes_3D_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
globals;
hF_cont = hF_cont + 1;
plot_axes_Earth(handles, hF_cont);


% --- Executes on mouse press over axes background.
function axes_OE_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_OE_ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dig = get(hObject, 'Tag'); dig = dig(end-1:end);
plot_axes_OE(handles, next_hF(), dig)


%> ======================================================================
%> @brief plot Orbit elements graphs
%> @param handles Main handles struct
%> @param hF Flag (handle) for separate window
%> @param hA Indexies of axes (2-el string)
% ======================================================================
function plot_axes_OE(handles, hF, ij)
globals;

if hF 
    figure(hF); 
    hA = gca;
    XLab = 't, s';
else
    hA = eval(['handles.axes_OE_' ij]);
    set(handles.fig_main,'CurrentAxes', hA)
    XLab = '';
    YLab = '';
end

X1 = tmod; X2 = tmod; X3 = tmod; X4 = tmod;
Y1 = nan(1, Nmod); Y2 = nan(1, Nmod); Y3 = nan(1, Nmod); Y4 = nan(1, Nmod);
switch ij
    case '11'
        Y1 = Xist.sqrtA;
        YLab = 'a^{1/2}, m^{1/2}';        
    case '12'
        Y1 = Xist.e;
        YLab = 'e';
    case '13'
        Y1 = Xist.Crc;
        Y2 = Xist.Crs;
        YLab = 'C_{rc}, C_{rs}, m';        
    case '14'
        Y1 = Xist.r;
        Y4 = Xist.A;
        YLab = 'r, a, m'; 
    case '15'
    case '21'
        Y1 = Xist.i0;
        YLab = 'i_0, rad';
    case '22'
        Y1 = Xist.i_dot;
        YLab = 'i_{dot}, rad/s';        
    case '23'
        Y1 = Xist.Cic;
        Y2 = Xist.Cis;
        YLab = 'C_{ic}, C_{is}, rad';        
    case '24'
        Y1 = Xist.i;
        YLab = 'i, rad';
    case '25'
    case '31'
        Y1 = Xist.M0;
        Y2 = Xist.E;
        Y3 = Xist.theta;
        YLab = 'M_0, E, \theta, rad';        
    case '32'
        Y1 = Xist.omega;
        YLab = '\omega, rad';
    case '33'
        Y1 = Xist.Cuc;
        Y2 = Xist.Cus;
        YLab = 'C_{uc}, C_{us}, rad';         
    case '34'
        Y1 = Xist.u;
        YLab = 'u, rad'; 
    case '35'
    case '41'
        Y1 = Xist.Omega;
        YLab = '\Omega, rad';
    case '42'
        Y1 = Xist.Omega_dot;
        YLab = '\Omega_{dot}, rad/s';        
    case '43'
    case '44'
        Y1 = Xist.x0;
        Y2 = Xist.y0;
        Y3 = Xist.z0;
        YLab = 'x_0, y_0, z_0, m'; 
    case '45'
        
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
ly = ylabel(YLab);
lx = xlabel(XLab);
grid(hA, 'on');

if ~hF
    set(hA, 'XTick', []);
    set(hA, 'YTick', []);
    set(ly, 'FontSize', Font_Size);
end


% --- Executes on button press in pb_replot.
function pb_replot_Callback(hObject, eventdata, handles)
% hObject    handle to pb_replot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i = 1:5
    for j = 1:5
        plot_axes_OE(handles, 0, [num2str(i) num2str(j)]);
    end
end


% --- Executes on button press in pb_Test1.
function pb_Test1_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Test1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
globals;
figure
plot3(Xist.x0, Xist.y0, Xist.z0)

% --- Executes on button press in pb_Test2.
function pb_Test2_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Test2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
globals;
figure
plot(tmod, Xist.r - Xist.r(1), tmod, Xist.A - Xist.A(1))


% --- Executes on button press in pb_Track.
function pb_Track_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

globals;


% --- Executes on button press in pb_True.
function pb_True_Callback(hObject, eventdata, handles)
% hObject    handle to pb_True (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pan_Tracking, 'Visible', 'off');
set(handles.pan_Orbital_Elements, 'Visible', 'on');

% --- Executes on button press in pb_Tracking.
function pb_Tracking_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pan_Orbital_Elements, 'Visible', 'off');
set(handles.pan_Tracking, 'Visible', 'on');
