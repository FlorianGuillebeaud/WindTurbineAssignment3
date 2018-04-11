% Assignment 3 
% Question 1
close all
clear all
clc

%% Read Blade and airfoil Data %%
blade_data = xlsread('Blade_data');

modes = importdata('modeshapes.txt');
W3_100 = importdata('cylinder_ds.txt'); %100% CILINDER
W3_60  = importdata('FFA-W3-600_ds.txt'); %600
W3_48  = importdata('FFA-W3-480_ds.txt'); %480
W3_36  = importdata('FFA-W3-360_ds.txt'); %360
W3_30  = importdata('FFA-W3-301_ds.txt'); %301
W3_24  = importdata('FFA-W3-241_ds.txt'); %241

N_element = length(blade_data) ;

% Data about WT
H = 119 ; % hub height (m)
Ls = 7.1 ; % m
R = 89.17 ; % [m] Rotor radius
B = 3 ; % Number of blades
Vcut_in = 4 ; % [m/s] Cut in speed
Vcut_out = 25 ; % [m/s] Cut out speed

Theta_pitch = degtorad(-3.34) ; % [rad]
Theta_cone = 0 ; % [rad]
Theta_tilt = 0 ; % [rad]
Theta_yaw = 0 ; % [rad]
rho = 1.225 ; % [kg/m3] air mass density
% time data
delta_t = 0.02 ; % [s]
N = 100 ; % [s]
omega_1f=3.93; %rad/s
omega_1e=6.1;
omega_2f=11.28;
dr(1)=blade_data(1,1);

for i=2:N_element
    dr(i)=blade_data(i,1)-blade_data(i-1,1);
end

uy_1f=modes(:,2);
uz_1f=modes(:,3);
uy_1e=modes(:,4);
uz_1e=modes(:,5);
uy_2f=modes(:,6);
uz_2f=modes(:,7);
V_0=8;


%%


[py_1f,pz_1f,time_1f]=TURB_BEM_turb(H, Ls, R, B, omega_1f, V_0, rho, delta_t, N, N_element, Theta_pitch, Theta_cone, Theta_tilt, Theta_yaw);
[py_1e,pz_1e,time_1e]=TURB_BEM_turb(H, Ls, R, B, omega_1e, V_0, rho, delta_t, N, N_element, Theta_pitch, Theta_cone, Theta_tilt, Theta_yaw)
[py_2f,pz_2f,time_2f]=TURB_BEM_turb(H, Ls, R, B, omega_2f, V_0, rho, delta_t, N, N_element, Theta_pitch, Theta_cone, Theta_tilt, Theta_yaw)

GF1=trapz(py_1f.*uy_1f,dr)+trapz(pz_1f.*uz_1f,dr);
GF2=trapz(py_1e.*uy_1e,dr)+trapz(pz_1e.*uz_1e,dr);
GF3=trapz(py_2f.*uy_2f,dr)+trapz(pz_2f.*uz_2f,dr);
GF=[GF1;GF2;GF3];