try 
    close(handle_fig_main);
end

close all
clear
clc

globals;

Tmod = 6*60*60;  %[s]
dTmod = 15; % [s]
tmod = 0:dTmod:Tmod;
Nmod = length(tmod);

Init_arrays; % Memory allocation
Font_Size = 8; % Font size for output interface
mu_earth = 3.9860044e14; % Gravity constant
options_solve = optimset('Display','off');  % Turn off display for fsolve

true_traectory;

handle_fig_main = fig_main();
