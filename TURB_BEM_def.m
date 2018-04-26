function [py, pz]=TURB_BEM_def(N_blade, H, Ls, R, B, omega0, V_0, rho, delta_t, N, N_element, Theta_pitch0, Theta_cone, Theta_tilt, Theta_yaw)
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
else
    fprintf('No file for this wind velocity \n')
end

% fid2=fopen('sim2.bin'); % (fluctuating v component)
% fid3=fopen('sim3.bin'); % (fluctuating w component)
global u

n2=32;
n3=32;
uraw=fread(fid1,'single');
% vraw=fread(fid2,'single');
% wraw=fread(fid3,'single');
itael=0;
for i=1:n1
 for j=1:n2
 for k=1:n3
 itael=itael+1;
 u(i,j,k)=uraw(itael);
%  v(i,j,k)=vraw(itael);
%  w(i,j,k)=wraw(itael);
 end
 end
end

%EXCEEDS MATRIX DIMENSIONS

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
global W3_100 W3_60 W3_48 W3_36 W3_30 W3_24 blade_data %M_G omega_list 
global M K D
global uy_1f uz_1f uy_1e uz_1e uy_2f uz_2f


Uy_dot = zeros(N_element, N) ; 
Uz_dot = zeros(N_element, N) ;

% Theta_pitch_i(1)=Theta_pitch0;
Theta_pitch = Theta_pitch0; % [rad]

% omega_max=deg2rad(10); %rad/s
% omega_ref=1.08; %rad/s
omega = omega0 ;
% Theta_max=deg2rad(35);
% Theta_min=deg2rad(0); %He said -2 as a possible value in class

V0y = 0 ;
V0z = V_0 ;
Wy = zeros(B,N_element,N) ;
Wz =  zeros(B,N_element,N) ;
a_rem = zeros(B,N_element,N) ; 



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

x_dotdot = zeros(N, 3);
x_dot = zeros(N, 3) ; 
x = zeros(N, 3) ;

    %% Loop
for i=2:N

    time(i) = time(i-1) + delta_t ;
    Theta_wing1(i) = Theta_wing1(i-1) + omega*delta_t ; % blade 1
    Theta_wing2(i) = Theta_wing1(i) + 2*pi/3 ; % blade 2
    Theta_wing3(i) = Theta_wing1(i) + 4*pi/3 ; % blade 3
    
    % loop over each blade B
    for b=1:N_blade
        % b
        % loop over each element N_element
        for k=1:N_element

           
                Theta_wing=eval(['Theta_wing',num2str(b)]);
                u_turb = velocity_turbulence(blade_data(k),Theta_wing(i),i);

            
             if k==9 && b==1
                  u_turb9(i)=u_turb;
             end
            
            [Vrel_y, Vrel_z] = velocity_compute_turb_2(u_turb,b, blade_data(k), H, Ls, Wy(b,k,i-1), Wz(b,k,i-1), Theta_wing1(i), Theta_wing2(i), Theta_wing3(i),omega,V_0, Theta_cone, Uy_dot(k,i), Uz_dot(k,i));
            
            phi = atan(real(-Vrel_z)/real(Vrel_y)) ;
            alpha = radtodeg(phi - (-degtorad(blade_data(k,3)) + Theta_pitch)) ;
            % alpha

            % first method (doesn't take into account the dynamic stall 
            thick = [100, 60, 48, 36, 30.1, 24.1] ; 
            
            % Cl interpolation 
            cl1 = interp1(W3_100(:,1), W3_100(:,2), alpha) ;
            cl2 = interp1(W3_60(:,1), W3_60(:,2), alpha) ;
            cl3 = interp1(W3_48(:,1), W3_48(:,2), alpha);
            cl4 = interp1(W3_36(:,1), W3_36(:,2), alpha);
            cl5 = interp1(W3_30(:,1), W3_30(:,2), alpha);
            cl6 = interp1(W3_24(:,1), W3_24(:,2), alpha);
            cl_union = [cl1 cl2 cl3 cl4 cl5 cl6] ; 
            Cl = interp1(thick, cl_union, blade_data(k,4)) ;
           
            
            % Cd interpolation 
            cd1 = interp1(W3_100(:,1), W3_100(:,3), alpha) ;
            cd2 = interp1(W3_60(:,1), W3_60(:,3), alpha) ;
            cd3 = interp1(W3_48(:,1), W3_48(:,3), alpha);
            cd4 = interp1(W3_36(:,1), W3_36(:,3), alpha);
            cd5 = interp1(W3_30(:,1), W3_30(:,3), alpha);
            cd6 = interp1(W3_24(:,1), W3_24(:,3), alpha);
            cd_union = [cd1 cd2 cd3 cd4 cd5 cd6] ; 
            Cd = interp1(thick, cd_union, blade_data(k,4)) ;
                  
            Vrel_abs = sqrt(Vrel_y^2+Vrel_z^2) ;
            Lift = 0.5*rho*Vrel_abs^2*Cl*blade_data(k,2) ;
            Drag = 0.5*rho*Vrel_abs^2*Cd*blade_data(k,2) ;
            
            pz(i,k) = Lift*cos(phi) + Drag*sin(phi) ; % normal
            py(i,k) = Lift*sin(phi) - Drag*cos(phi) ; % tangential

            % without Yaw, a can be calculate as follow : 
            a = abs(Wz(b,k,i-1))/V_0 ;
            a_rem(b,k,i) = a ;
            % with yaw : need to be implemented 
            if a<=1/3
                fg = 1 ;
            else
                fg = 1/4*(5-3*a) ;
            end
            
            % Prand
            f = (B/2)*(R-blade_data(k))/(blade_data(k)*abs(sin(phi)));
            F= 2*acos(exp(-f))/pi;
             
           
            % We add this if statement otherwise the last element if NaN 
            % (F = 0 !! )
            if k==N_element
                Wz(b,k,i) = 0 ; 
                Wy(b,k,i) = 0 ; 
            else   
                Wz(b,k,i) = - B*Lift*cos(phi)/(4*pi*rho*blade_data(k)*F*(sqrt(V0y^2+(V0z+u_turb+fg*Wz(b,k,i-1))^2))) ;
                Wy(b,k,i) = - B*Lift*sin(phi)/(4*pi*rho*blade_data(k)*F*(sqrt(V0y^2+(V0z+u_turb+fg*Wz(b,k,i-1))^2))) ;
            end
          dm_edge(k,i) = blade_data(k)*py(i,k) ;
          dm_flap(k,i)= blade_data(k)*pz(i,k);
            dP(k) = omega*dm_edge(k) ;
            
        end
        pz(i,N_element) = 0 ;
        py(i,N_element) = 0 ; 
        dm_edge(N_element,N) = 0 ;
        dm_flap(N_element,N)=0;
        dP(N_element) = 0 ;
        M_edge(:,i)=trapz(blade_data(:,1), real(dm_edge(:,i))) ;
        M_flap(:,i)=trapz(blade_data(:,1), real(dm_flap(:,i)));
    end
    
    GF1(i)=trapz(py(i,:)'.*uy_1f,dr)+trapz(pz(i,:)'.*uz_1f,dr);
    GF2(i)=trapz(py(i,:)'.*uy_1e,dr)+trapz(pz(i,:)'.*uz_1e,dr);
    GF3(i)=trapz(py(i,:)'.*uy_2f,dr)+trapz(pz(i,:)'.*uz_2f,dr);
    GF(:,i)=[GF1(i);GF2(i);GF3(i)];

    GF_loc = GF(:,i);
    x_dotdot(i,:) = (inv(M)*(GF_loc-D*x_dot(i,:)'-K*x(i,:)'))' ; 
    A = 0.5*delta_t*x_dotdot(i,:) ; 
    b = 0.5*delta_t*(x_dot(i)+0.5*A);
    
    x_dotnew = x_dot(i,:)+A;
    x_new = x(i,:)+b;
    x_dotdot_new = (inv(M)*(GF_loc-D*x_dotnew'-K*x_new'))';
    B = 0.5*delta_t*x_dotdot_new ;
    
    x_dotnew = x_dot(i,:)+B;
    x_dotdot_new = (inv(M)*(GF_loc-D*x_dotnew'-K*x_new'))';
    C = 0.5*delta_t*x_dotdot_new ;
    
    d = delta_t*(x_dot(i,:)+C);
    x_dotnew = x_dot(i,:)+2*C;
    x_new = x(i,:)+d;
    x_dotdot_new = (inv(M)*(GF_loc-D*x_dotnew'-K*x_new'))';
    
    DD = 0.5*delta_t*x_dotdot_new ; 
    
    x(i+1,:) = x(i,:) + delta_t*(x_dot(i,:)+(1/3)*(A+B+C));
    x_dot(i+1,:) = x_dot(i,:) + (1/3)*(A+2*B+2*C+DD);
    
%% 
Uy_dot(:,i)=x_dot(i+1,1)'.*uy_1f+x_dot(i+1,2)'.*uy_1e+x_dot(i+1,3)'.*uy_2f;
Uz_dot(:,i)=x_dot(i+1,1)'.*uz_1f+x_dot(i+1,2)'.*uz_1e+x_dot(i+1,3)'.*uz_2f;

end
end