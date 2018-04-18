% Assignment 3 
% Question 1
close all
clear all
clc

%% Read Blade and airfoil Data %%
blade_data = xlsread('Blade_data');
global blade_data W3_100 W3_60 W3_48 W3_36 W3_30 W3_24
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

delta = 0.03 ; %damping factor
% time data
delta_t = 0.02 ; % [s]
N = 100 ; % [points]
omega_1f = 3.93; %rad/s
omega_1e=6.1;
omega_2f=11.28;
lambda=8;
omega0 = lambda*V_0/R ; 
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

m=blade_data(:,5)'.*dr;
%% GF
[py,pz,time]=TURB_BEM_turb(H, Ls, R, B, omega0, V_0, rho, delta_t, N, N_element, Theta_pitch, Theta_cone, Theta_tilt, Theta_yaw);
% [py_1e,pz_1e,time_1e]=TURB_BEM_turb(H, Ls, R, B, omega_1e, V_0, rho, delta_t, N, N_element, Theta_pitch, Theta_cone, Theta_tilt, Theta_yaw)
% [py_2f,pz_2f,time_2f]=TURB_BEM_turb(H, Ls, R, B, omega_2f, V_0, rho, delta_t, N, N_element, Theta_pitch, Theta_cone, Theta_tilt, Theta_yaw)

for i=1:length(time)
    GF1(i)=trapz(py(i,:)'.*uy_1f,dr)+trapz(pz(i,:)'.*uz_1f,dr);
    GF2(i)=trapz(py(i,:)'.*uy_1e,dr)+trapz(pz(i,:)'.*uz_1e,dr);
    GF3(i)=trapz(py(i,:)'.*uy_2f,dr)+trapz(pz(i,:)'.*uz_2f,dr);
    GF(:,i)=[GF1(i);GF2(i);GF3(i)];
end

%% M

M1=trapz(uy_1f'.*m.*uy_1f',dr)+trapz(uz_1f'.*m.*uz_1f',dr);
M2=trapz(uy_1e'.*m.*uy_1e',dr)+trapz(uz_1e'.*m.*uz_1e',dr);
M3=trapz(uy_2f'.*m.*uy_2f',dr)+trapz(uz_2f'.*m.*uz_2f',dr);
M=eye(3).*[M1 M2 M3];
%K
K=eye(3).*[omega_1f^2*M1 omega_1e^2*M2 omega_2f^2*M3];
D=eye(3)*delta.*[omega_1f*M1 omega_1e*M2 omega_2f*M3]./pi;

%GF(varying in time)=M*x_dot_dot+D*x_dot+K*x;
%% 
global Uy_dot Uz_dot
Uy_dot=x_dot(1)*uy_1f+x_dot(2)*uy_1e+x_dot(3)*uy_2f;
Uz_dot=x_dot(1)*uz_1f+x_dot(2)*uz_1e+x_dot(3)*uz_2f;

[Vrel_y, Vrel_z] = TURB_Vrel(H, Ls, R, B, omega0, V_0, rho, delta_t, N, N_element, Theta_pitch, Theta_cone, Theta_tilt, Theta_yaw, Uy_dot, Uz_dot);