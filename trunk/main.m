try 
    close(handle_fig_main);
end

close all
clear 
clc
globals;

addpath([pwd '/func/tracking']);
addpath([pwd '/func/interface']);

Tmod = 3.2*60*60;  %[s], duration of the simulation

hF_cont = 0; % Figure's handles
Font_Size = 8; % Font size for output interface
mu_earth = 3.9860044e14; % Gravity constant
options_solve = optimset('Display','off');  % Turn off display for fsolve

load TrueTrajectory.mat
Nmod = fix(Tmod/T);
if Nmod > Nmod_max
    Nmod = Nmod_max;
    Tmod = (Nmod-1)*T;
end
resize_arrays; % resize arrays for new Tmod

handle_fig_main = fig_main(); % open GUI form
