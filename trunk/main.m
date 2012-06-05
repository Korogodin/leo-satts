try 
    close(handle_fig_main);
end

close all
clear
clc

globals;

h = 600e3;
r_e = 6371e3;
N_r = 2.6561e+07 / (r_e + h);

Tmod = 24*60*60;  %[s]
T = 1;
dTmod = T * N_r^(3/2); % [s]
tmod = 0:dTmod:Tmod;
Nmod = length(tmod);

Init_arrays; % Memory allocation
Font_Size = 8; % Font size for output interface
mu_earth = 3.9860044e14; % Gravity constant
options_solve = optimset('Display','off');  % Turn off display for fsolve

true_traectory;

handle_fig_main = fig_main();
