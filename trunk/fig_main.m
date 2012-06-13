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

% Last Modified by GUIDE v2.5 13-Jun-2012 12:18:01

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
set(handles.pan_Sch, 'Position', get(handles.pan_Orbital_Elements, 'Position'));

% Draw functional scheme
fig_main_pictureData = imread('Sch1.png');
image(fig_main_pictureData, 'Parent', handles.axes_Sch1);
set(handles.axes_Sch1, 'XTick', []);
set(handles.axes_Sch1, 'YTick', []);

% Set widget's function
for i = 1:5
    for j = 1:5
        set(eval(['handles.axes_OE_' num2str(i) num2str(j)]), 'DrawMode', 'fast');
        plot_axes_OE(handles, 0, [num2str(i) num2str(j)]);
        set(eval(['handles.axes_OE_' num2str(i) num2str(j)]), 'Tag', ['axes_OE_' num2str(i) num2str(j)]);
        set(eval(['handles.axes_OE_' num2str(i) num2str(j)]), 'ButtonDownFcn', str2func('@(hObject,eventdata)fig_main(''axes_OE_ButtonDownFcn'',hObject,eventdata,guidata(hObject))'));
    end
end
for i = 1:5
    for j = 1:5
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'DrawMode', 'fast');
        plot_axes_Track(handles, 0, [num2str(i) num2str(j)]);
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'Tag', ['axes_Track_' num2str(i) num2str(j)]);
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'ButtonDownFcn', str2func('@(hObject,eventdata)fig_main(''axes_Track_ButtonDownFcn'',hObject,eventdata,guidata(hObject))'));
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


% --- Executes on mouse press over axes background.
function axes_Track_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_OE_ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dig = get(hObject, 'Tag'); dig = dig(end-1:end);
plot_axes_Track(handles, next_hF(), dig)


% --- Executes on button press in pb_Track.
function pb_Track_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
std_x = str2double(get(handles.ed_stdX, 'String'));
if ~((std_x >= 0)&&(std_x < 1000))
    msgbox('Coordinate''s RMS is incorrect', 'Error', 'error');
    set(handles.ed_stdX, 'String', '10');
    return
end
std_V = str2double(get(handles.ed_stdV, 'String'));
if ~((std_V >= 0)&&(std_V < 5))
    msgbox('Velocity''s RMS is incorrect', 'Error', 'error');
    set(handles.ed_stdV, 'String', '0.01');
    return
end
track_with_noise(handles, std_x, std_V);

% --- Executes on button press in pb_True.
function pb_True_Callback(hObject, eventdata, handles)
% hObject    handle to pb_True (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pan_Orbital_Elements, 'Visible', 'on');
set(handles.pan_Tracking, 'Visible', 'off');
set(handles.pan_Sch, 'Visible', 'off');


% --- Executes on button press in pb_Tracking.
function pb_Tracking_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pan_Orbital_Elements, 'Visible', 'off');
set(handles.pan_Tracking, 'Visible', 'on');
set(handles.pan_Sch, 'Visible', 'off');


% --- Executes on button press in pb_Solve.
function pb_Solve_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Solve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
solve_wo_noise(handles);



function ed_stdX_Callback(hObject, eventdata, handles)
% hObject    handle to ed_stdX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_stdX as text
%        str2double(get(hObject,'String')) returns contents of ed_stdX as a double


% --- Executes during object creation, after setting all properties.
function ed_stdX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_stdX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_stdV_Callback(hObject, eventdata, handles)
% hObject    handle to ed_stdV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_stdV as text
%        str2double(get(hObject,'String')) returns contents of ed_stdV as a double


% --- Executes during object creation, after setting all properties.
function ed_stdV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_stdV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_Sch.
function pb_Sch_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Sch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.pan_Orbital_Elements, 'Visible', 'off');
set(handles.pan_Tracking, 'Visible', 'off');
set(handles.pan_Sch, 'Visible', 'on');
