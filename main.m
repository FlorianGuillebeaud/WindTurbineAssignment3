% Assignment 3 
% Question 1
close all
clear all
clc

%% Read Blade and airfoil Data %%
blade_data = xlsread('Blade_data');
global blade_data W3_100 W3_60 W3_48 W3_36 W3_30 W3_24 N_element 
modes = importdata('modeshapes.txt');
W3_100 = importdata('cylinder_ds.txt'); %100% CILINDER
W3_60  = importdata('FFA-W3-600_ds.txt'); %600
W3_48  = importdata('FFA-W3-480_ds.txt'); %480
W3_36  = importdata('FFA-W3-360_ds.txt'); %360
W3_30  = importdata('FFA-W3-301_ds.txt'); %301
W3_24  = importdata('FFA-W3-241_ds.txt'); %241
N_element = length(blade_data) ;

% Data about WT
global H Ls R blades
H = 119 ; % hub height (m)
Ls = 7.1 ; % m
R = 89.17 ; % [m] Rotor radius
blades = 3 ; % Number of blades

global Theta_pitch Theta_cone Theta_tilt Theta_yaw V_0 rho N delta_t omega_1f omega_1e omega_2f omega0
Theta_pitch = degtorad(-3.34) ; % [rad]
Theta_cone = 0 ; % [rad]
Theta_tilt = 0 ; % [rad]
Theta_yaw = 0 ; % [rad] 
V_0=8;
rho = 1.225 ; % [kg/m3] air mass density
N = 8000 ; % [points]
delta = 0.03 ; %damping factor
% time data
delta_t = 0.02 ; % [s]
omega_1f = 3.93; %rad/s
omega_1e=6.1;
omega_2f=11.28;
lambda=7.5;
omega0 = lambda*V_0/R ; 
dr(1)=blade_data(1,1);

for i=2:N_element
    dr(i)=blade_data(i,1)-blade_data(i-1,1);
end

global uy_1f uz_1f uy_1e uz_1e uy_2f uz_2f
uy_1f=modes(:,2);
uz_1f=modes(:,3);
uy_1e=modes(:,4);
uz_1e=modes(:,5);
uy_2f=modes(:,6);
uz_2f=modes(:,7);


m = blade_data(:,5)';

%% M
global M K D

M1=trapz(blade_data(:,1), uy_1f'.*m.*uy_1f')+trapz(blade_data(:,1), uz_1f'.*m.*uz_1f');
M2=trapz(blade_data(:,1), uy_1e'.*m.*uy_1e')+trapz(blade_data(:,1), uz_1e'.*m.*uz_1e');
M3=trapz(blade_data(:,1), uy_2f'.*m.*uy_2f')+trapz(blade_data(:,1), uz_2f'.*m.*uz_2f');
M=eye(3).*[M1 M2 M3];
%K
K=eye(3).*[omega_1f^2*M1 omega_1e^2*M2 omega_2f^2*M3];
D=eye(3)*delta.*[omega_1f*M1 omega_1e*M2 omega_2f*M3]./pi;
%D = 0 ;
%% what output do we want have ? so far we have py and pz but might not be relevant
N_blade=1;
[Vrel_y,Vrel_z, x, x_dotdot, M_edge, M_flap, time,py,pz]=BEM_turb(N_blade);


%% Q2
Uyf_tip=x(:,1)*uy_1f';
Uzf_tip=x(:,1)*uz_1f';
Uye_tip=x(:,2)*uy_1e';
Uze_tip=x(:,2)*uz_1e';
Uy2f_tip=x(:,3)*uy_2f';
Uz2f_tip=x(:,3)*uz_2f';

figure()
plot(time, Uyf_tip(2:end,18)) %maximum one: last element
hold on
plot(time, Uye_tip(2:end,18))
hold on
plot(time, Uy2f_tip(2:end,18))
xlabel('Time (s)')
ylabel('Deformation (y axis)')
legend('Flapwise Def.', 'Edgewise Def.','2nd Flapwise Def.')
hold off

figure()
plot(time, Uzf_tip(2:end,18))
hold on
line([0, time(length(time))], [mean(Uzf_tip(2:end,18)), mean(Uzf_tip(2:end,18))])
hold on
plot(time, Uze_tip(2:end,18))
hold on
plot(time, Uy2f_tip(2:end,18))
xlabel('Time (s)')
ylabel('Deformation (z axis)')
legend('Flapwise Def.', 'Edgewise Def.', '2nd Flapwise Def.')
hold off
%% 
Fs=1/delta_t;
f = Fs*(0:(N/2))/N;
Y = fft(Uze_tip(1:end-1,18));
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure()
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of Uze_tip')
xlabel('f (Hz)')
ylabel('|P1(f)|')
%Mean deflections
% EIs=load ('bladestruc.txt');
% global r EI rglobal EI1global EI2global betaRadglobal
% r=blade_data(:,1);
% % EI_flap=EIs(:,4);
% EI_edge=EIs(:,5);
% rglobal = r;
% EI1global = EI_edge*ones(1,length(r)); 
% EI2global = EI_edge*ones(1,length(r));
% betaRadglobal = zeros(1,length(r));

% [py,pz]=TURB_BEM_def(N_blade, H, Ls, R, blades, omega0, V_0, rho, delta_t, N, N_element, Theta_pitch, Theta_cone, Theta_tilt, Theta_yaw);
% [uy, uz]=beam_deflection(py,pz,Theta_pitch);

%bending moment
% figure()
% plot(time, M_flap)
% ylabel('Flapwise bending moment (Nm)')
% xlabel('Time (s)')
% 
% figure()
% plot(time, M_edge)
% ylabel('Edgewise bending moment (Nm)')
% xlabel('Time (s)')

%% Q3 
global  M10dof D10dof K10dof
N_blade=3;

M1=trapz(blade_data(:,1), uy_1f'.*m.*uy_1f')+trapz(blade_data(:,1), uz_1f'.*m.*uz_1f');
M2=trapz(blade_data(:,1), uy_1e'.*m.*uy_1e')+trapz(blade_data(:,1), uz_1e'.*m.*uz_1e');
M3=trapz(blade_data(:,1), uy_2f'.*m.*uy_2f')+trapz(blade_data(:,1), uz_2f'.*m.*uz_2f');
Mnacelle = 446000 ; %kg

M10dof = [Mnacelle+3*sum(m) trapz(blade_data(:,1), uz_1f'.*m) trapz(blade_data(:,1), uz_1e'.*m) trapz(blade_data(:,1),uz_2f'.*m) trapz(blade_data(:,1),uz_1f'.*m) trapz(blade_data(:,1), uz_1e'.*m) trapz(blade_data(:,1), uz_2f'.*m) trapz(blade_data(:,1), uz_1f'.*m) trapz(blade_data(:,1), uz_1e'.*m) trapz(blade_data(:,1), uz_2f'.*m);
        trapz(blade_data(:,1), m.*uz_1f') M1 0 0 0 0 0 0 0 0;
        trapz(blade_data(:,1), m.*uz_1e') 0 M2 0 0 0 0 0 0 0;
        trapz(blade_data(:,1), m.*uz_2f') 0 0 M3 0 0 0 0 0 0;
        trapz(blade_data(:,1), m.*uz_1f') 0 0 0 M1 0 0 0 0 0;
        trapz(blade_data(:,1), m.*uz_1e') 0 0 0 0 M2 0 0 0 0;
        trapz(blade_data(:,1), m.*uz_2f') 0 0 0 0 0 M3 0 0 0;
        trapz(blade_data(:,1), m.*uz_1f') 0 0 0 0 0 0 M1 0 0;
        trapz(blade_data(:,1), m.*uz_1e') 0 0 0 0 0 0 0 M2 0;
        trapz(blade_data(:,1), m.*uz_2f') 0 0 0 0 0 0 0 0 M3];

%K
K1=omega_1f^2*M1;
K2=omega_1e^2*M2;
K3=omega_2f^2*M3;
k=1.7*10^6; %N/m stifness of the tower

K10dof = [1.7*10^6 0 0 0 0 0 0 0 0 0 ;
           0 K1 0 0 0 0 0 0 0 0;
           0 0 K2 0 0 0 0 0 0 0;
           0 0 0 K3 0 0 0 0 0 0;
           0 0 0 0 K1 0 0 0 0 0;
           0 0 0 0 0 K2 0 0 0 0;
           0 0 0 0 0 0 K3 0 0 0;
           0 0 0 0 0 0 0 K1 0 0;
           0 0 0 0 0 0 0 0 K2 0;
           0 0 0 0 0 0 0 0 0 K3];


D1=delta/pi*omega_1f*M1;
D2=delta/pi*omega_1e*M2;
D3=delta/pi*omega_2f*M3;

D10dof = [0 0 0 0 0 0 0 0 0 0 ;
    0 D1 0 0 0 0 0 0 0 0 ;
    0 0 D2 0 0 0 0 0 0 0 ; 
    0 0 0 D3 0 0 0 0 0 0 ;
    0 0 0 0 D1 0 0 0 0 0 ;
    0 0 0 0 0 D2 0 0 0 0 ;
    0 0 0 0 0 0 D3 0 0 0 ;
    0 0 0 0 0 0 0 D1 0 0 ;
    0 0 0 0 0 0 0 0 D2 0 ;
    0 0 0 0 0 0 0 0 0 D3 ];
%%
N_blade = 3 ;
[Vrel_y,Vrel_z, x, M_edge, M_flap, time,py,pz]=BEM_turb10dof(N_blade);


%% 
% Blade 1
Uyf_tip=x(:,2)*uy_1f';
Uzf_tip=x(:,2)*uz_1f';
Uye_tip=x(:,3)*uy_1e';
Uze_tip=x(:,3)*uz_1e';
Uy2f_tip=x(:,4)*uy_2f';
Uz2f_tip=x(:,4)*uz_2f';

% Blade 2
Uyf_tip=x(:,5)*uy_1f';
Uzf_tip=x(:,5)*uz_1f';
Uye_tip=x(:,6)*uy_1e';
Uze_tip=x(:,6)*uz_1e';
Uy2f_tip=x(:,7)*uy_2f';
Uz2f_tip=x(:,7)*uz_2f';

% Blade 3
Uyf_tip=x(:,8)*uy_1f';
Uzf_tip=x(:,8)*uz_1f';
Uye_tip=x(:,9)*uy_1e';
Uze_tip=x(:,9)*uz_1e';
Uy2f_tip=x(:,10)*uy_2f';
Uz2f_tip=x(:,10)*uz_2f';
% 
figure()
plot(time, Uyf_tip(1:end-1,18)) %maximum one: last element
hold on
plot(time, Uye_tip(1:end-1,18))
xlabel('Time (s)')
ylabel('Deformation (y axis)')
legend('Flapwise Def.', 'Edgewise Def.')
hold off
% 
figure()
plot(time, Uzf_tip(1:end-1,18))
hold on
plot(time, Uze_tip(1:end-1,18))
xlabel('Time (s)')
ylabel('Deformation (z axis)')
legend('Flapwise Def.', 'Edgewise Def.')

% hold off
% 
% 
% %bending moment
% figure()
% plot(time3, M_flap3)
% ylabel('Flapwise bending moment (Nm)')
% xlabel('Time (s)')
% 
% figure()
% plot(time3, M_edge3)
% ylabel('Edgewise bending moment (Nm)')
% xlabel('Time (s)')
