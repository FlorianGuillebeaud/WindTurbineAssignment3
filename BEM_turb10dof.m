function [Vrel_y, Vrel_z, x, x_dotdot, M_edge, M_flap, time, py, pz] = BEM_turb(N_blade)

global V_0 Theta_pitch0 blades N_element N omega0 Theta_tilt Theta_yaw Theta_cone delta_t u

%% Read the binary files
% WE NEED TO HAVE DIFFERENT FILES FOR EACH VELOCITY
if V_0==15
    fid1=fopen('sim1.bin');% (fluctuating u component)
    n1=16384;
elseif V_0==7
    fid1=fopen('sim2.bin');
    n1=16384;
elseif V_0==8
    fid1=fopen('simV08.bin');
    n1=8192;
    fprintf('Read file simV08 \n')
else
    fprintf('No file for this wind velocity \n')
end
n2=32;
n3=32;
uraw=fread(fid1,'single');
itael=0;
for i=1:n1
 for j=1:n2
 for k=1:n3
 itael=itael+1;
 u(i,j,k)=uraw(itael);
 end
 end
end

% grid
global Mx My
n_point=32;
size_grid=200;
dxy=size_grid/n_point;

for i=1:n_point
    Mx(i)=-size_grid/2+i*dxy;
    My(i)=-size_grid/2+i*dxy;
end

%Initialization 
global M10dof K10dof D10dof blade_data uy_1f uz_1f uy_1e uz_1e uy_2f uz_2f V0y V0z
Theta_pitch = Theta_pitch0; % [rad]
V0y = 0 ;
V0z = V_0 ;
Wy = zeros(blades,N_element,N) ;
Wz =  zeros(blades,N_element,N) ;
a_rem = zeros(blades,N_element,N) ; 

global a_12 a_21 a_34 
a_12 = [cos(Theta_tilt) sin(Theta_yaw)*sin(Theta_tilt) -cos(Theta_yaw)*sin(Theta_tilt) ;
    0 cos(Theta_yaw) sin(Theta_yaw) ;
    sin(Theta_tilt) -sin(Theta_yaw)*cos(Theta_tilt) cos(Theta_yaw)*cos(Theta_tilt)] ;

a_21 = a_12' ;

a_34 = [cos(Theta_cone) 0 -sin(Theta_cone) ;
    0 1 0 ;
    sin(Theta_cone) 0 cos(Theta_cone)];

a_43 = a_34' ;


time(1) = 0 ;

% Wind position initialization
Theta_wing1(1) = 0 ; % blade 1
Theta_wing2(1) = Theta_wing1(1) + 2*pi/3 ; % blade 2
Theta_wing3(1) = Theta_wing1(1) + 4*pi/3 ; % blade 3

% blade element
dr(1)=blade_data(1,1);
for i=2:N_element
    dr(i)=blade_data(i,1)-blade_data(i-1,1);
end

if N_blade==3
    x_dotdot = zeros(N, 10);
    x_dot = zeros(N, 10) ; 
    x = zeros(N, 10) ;
    x(2,2) = 1;
    x(2,5) = 1 ;
    x(2,8) = 1 ;
    Uy_dot = zeros(N, N_element, N_blade) ; 
    Uz_dot = zeros(N, N_element, N_blade) ;
end

py=zeros(N,N_element,N_blade);
pz=zeros(N,N_element,N_blade);
M_edge=zeros(N,N_element, N_blade);
M_flap=zeros(N,N_element, N_blade);
    %% Loop
for i=2:N
    i 
    time(i) = time(i-1) + delta_t ;
    Theta_wing1(i) = Theta_wing1(i-1) + omega0*delta_t ; % blade 1
    Theta_wing2(i) = Theta_wing1(i) + 2*pi/3 ; % blade 2
    Theta_wing3(i) = Theta_wing1(i) + 4*pi/3 ; % blade 3
    
    % Step 1 KUNGA
    [GF, Vrel_y, Vrel_z, M_edge(i,:,:), M_flap(i,:,:), py(i,:,:), pz(i,:,:)] = GF_compute(i, Uy_dot(i-1,:,:), Uz_dot(i-1,:,:), N_blade, Theta_wing1, Theta_wing2, Theta_wing3, Wy, Wz) ;
    GF_loc = GF(:,i);
    x_dotdot(i,:) = (inv(M10dof)*(GF_loc-D10dof*x_dot(i,:)'-K10dof*x(i,:)'))' ; 
    A = 0.5*delta_t*x_dotdot(i,:) ; 
    b = 0.5*delta_t*(x_dot(i,:)+0.5*A);
    x_dotnew = x_dot(i,:)+A;
    x_new = x(i,:)+b;
    
    Uy_dot(i, :,1)=x_dotnew(2)'.*uy_1f+x_dotnew(3)'.*uy_1e+x_dotnew(4)'.*uy_2f;
    Uz_dot(i, :,1)=x_dotnew(2)'.*uz_1f+x_dotnew(3)'.*uz_1e+x_dotnew(4)'.*uz_2f;
    Uy_dot(i,:,2)=x_dotnew(5)'.*uy_1f+x_dotnew(6)'.*uy_1e+x_dotnew(7)'.*uy_2f;
    Uz_dot(i,:,2)=x_dotnew(5)'.*uz_1f+x_dotnew(6)'.*uz_1e+x_dotnew(7)'.*uz_2f;
    Uy_dot(i,:,3)=x_dotnew(8)'.*uy_1f+x_dotnew(9)'.*uy_1e+x_dotnew(10)'.*uy_2f;
    Uz_dot(i,:,3)=x_dotnew(8)'.*uz_1f+x_dotnew(9)'.*uz_1e+x_dotnew(10)'.*uz_2f;
       
    
    [GF, Vrel_y, Vrel_z, M_edge(i,:,:), M_flap(i,:,:) , py(i,:,:), pz(i,:,:)] = GF_compute(i, Uy_dot(i,:,:), Uz_dot(i,:,:), N_blade, Theta_wing1, Theta_wing2, Theta_wing3, Wy, Wz) ;
    GF_loc = GF(:,i);
    
    % Step 2 Kunga
    x_dotdot_new = (inv(M10dof)*(GF_loc-D10dof*x_dotnew'-K10dof*x_new'))';
    B = 0.5*delta_t*x_dotdot_new ;
    x_dotnew = x_dot(i,:)+B;
    
    Uy_dot(i, :,1)=x_dotnew(2)'.*uy_1f+x_dotnew(3)'.*uy_1e+x_dotnew(4)'.*uy_2f;
    Uz_dot(i, :,1)=x_dotnew(2)'.*uz_1f+x_dotnew(3)'.*uz_1e+x_dotnew(4)'.*uz_2f;
    Uy_dot(i,:,2)=x_dotnew(5)'.*uy_1f+x_dotnew(6)'.*uy_1e+x_dotnew(7)'.*uy_2f;
    Uz_dot(i,:,2)=x_dotnew(5)'.*uz_1f+x_dotnew(6)'.*uz_1e+x_dotnew(7)'.*uz_2f;
    Uy_dot(i,:,3)=x_dotnew(8)'.*uy_1f+x_dotnew(9)'.*uy_1e+x_dotnew(10)'.*uy_2f;
    Uz_dot(i,:,3)=x_dotnew(8)'.*uz_1f+x_dotnew(9)'.*uz_1e+x_dotnew(10)'.*uz_2f;
       
    [GF, Vrel_y, Vrel_z, M_edge(i,:,:), M_flap(i,:,:) , py(i,:,:), pz(i,:,:)] = GF_compute(i, Uy_dot(i,:,:), Uz_dot(i,:,:), N_blade, Theta_wing1, Theta_wing2, Theta_wing3, Wy, Wz) ;
    GF_loc = GF(:,i);
    
    
    % Step 3 Kunga
    x_dotdot_new = (inv(M10dof)*(GF_loc-D10dof*x_dotnew'-K10dof*x_new'))';
    C = 0.5*delta_t*x_dotdot_new ;
    d = delta_t*(x_dot(i,:)+C);
    x_dotnew = x_dot(i,:)+2*C;
    x_new = x(i,:)+d;
    
    Uy_dot(i, :,1)=x_dotnew(2)'.*uy_1f+x_dotnew(3)'.*uy_1e+x_dotnew(4)'.*uy_2f;
    Uz_dot(i, :,1)=x_dotnew(2)'.*uz_1f+x_dotnew(3)'.*uz_1e+x_dotnew(4)'.*uz_2f;
    Uy_dot(i,:,2)=x_dotnew(5)'.*uy_1f+x_dotnew(6)'.*uy_1e+x_dotnew(7)'.*uy_2f;
    Uz_dot(i,:,2)=x_dotnew(5)'.*uz_1f+x_dotnew(6)'.*uz_1e+x_dotnew(7)'.*uz_2f;
    Uy_dot(i,:,3)=x_dotnew(8)'.*uy_1f+x_dotnew(9)'.*uy_1e+x_dotnew(10)'.*uy_2f;
    Uz_dot(i,:,3)=x_dotnew(8)'.*uz_1f+x_dotnew(9)'.*uz_1e+x_dotnew(10)'.*uz_2f;
       
    
    [GF, Vrel_y, Vrel_z, M_edge(i,:,:), M_flap(i,:,:) , py(i,:,:), pz(i,:,:)] = GF_compute(i, Uy_dot(i,:,:), Uz_dot(i,:,:), N_blade, Theta_wing1, Theta_wing2, Theta_wing3, Wy, Wz) ;
    GF_loc = GF(:,i);
    
    
    x_dotdot_new = (inv(M10dof)*(GF_loc-D10dof*x_dotnew'-K10dof*x_new'))';
    DD = 0.5*delta_t*x_dotdot_new ; 
    x(i+1,:) = x(i,:) + delta_t*(x_dot(i,:)+(1/3)*(A+B+C));
    x_dot(i+1,:) = x_dot(i,:) + (1/3)*(A+2*B+2*C+DD);
    
%     Uy_dot(i, :,1)=x_dotnew(2)'.*uy_1f+x_dotnew(3)'.*uy_1e+x_dotnew(4)'.*uy_2f;
%     Uz_dot(i, :,1)=x_dotnew(2)'.*uz_1f+x_dotnew(3)'.*uz_1e+x_dotnew(4)'.*uz_2f;
%     Uy_dot(i,:,2)=x_dotnew(5)'.*uy_1f+x_dotnew(6)'.*uy_1e+x_d otnew(7)'.*uy_2f;
%     Uz_dot(i,:,2)=x_dotnew(5)'.*uz_1f+x_dotnew(6)'.*uz_1e+x_dotnew(7)'.*uz_2f;
%     Uy_dot(i,:,3)=x_dotnew(8)'.*uy_1f+x_dotnew(9)'.*uy_1e+x_dotnew(10)'.*uy_2f;
%     Uz_dot(i,:,3)=x_dotnew(8)'.*uz_1f+x_dotnew(9)'.*uz_1e+x_dotnew(10)'.*uz_2f;
    Uy_dot(i, :,1)=x_dot(i+1,2)'.*uy_1f+x_dot(i+1,3)'.*uy_1e+x_dot(i+1,4)'.*uy_2f;
    Uz_dot(i, :,1)=x_dot(i+1,2)'.*uz_1f+x_dot(i+1,3)'.*uz_1e+x_dot(i+1,4)'.*uz_2f;
    Uy_dot(i,:,2)=x_dot(i+1,5)'.*uy_1f+x_dot(i+1,6)'.*uy_1e+x_dot(i+1,7)'.*uy_2f;
    Uz_dot(i,:,2)=x_dot(i+1,5)'.*uz_1f+x_dot(i+1,6)'.*uz_1e+x_dot(i+1,7)'.*uz_2f;
    Uy_dot(i,:,3)=x_dot(i+1,8)'.*uy_1f+x_dot(i+1,9)'.*uy_1e+x_dot(i+1,10)'.*uy_2f;
    Uz_dot(i,:,3)=x_dot(i+1,8)'.*uz_1f+x_dot(i+1,9)'.*uz_1e+x_dot(i+1,10)'.*uz_2f;
    
 

end
% end
