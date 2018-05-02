function [GF, Vrel_y, Vrel_z, M_edge, M_flap ,py,pz] = GF_compute(i, Uy_dot, Uz_dot, N_blade, Theta_wing1, Theta_wing2, Theta_wing3, Wy, Wz)

global V_0 blades N_element N Theta_pitch blade_data rho R omega0 V0y V0z
global W3_100 W3_60 W3_48 W3_36 W3_30 W3_24
global uy_1f uz_1f uy_1e uz_1e uy_2f uz_2f

for b=1:N_blade
    % loop over each element N_element
    for k=1:N_element
        
        Theta_wing=eval(['Theta_wing',num2str(b)]);
        u_turb = velocity_turbulence(k, Theta_wing(i),i);
       
        [Vrel_y, Vrel_z] = velocity_compute_turb(u_turb, b, k, Wy(b,k,i-1), Wz(b,k,i-1), Theta_wing1(i), Theta_wing2(i), Theta_wing3(i), Uy_dot(k), Uz_dot(k));
        
        phi = atan((-Vrel_z)/(Vrel_y)) ;
        alpha = radtodeg(phi - (-degtorad(blade_data(k,3)) + Theta_pitch)) ;
        
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
        
        pz(i,k,b) = Lift*cos(phi) + Drag*sin(phi) ; % normal
        py(i,k,b) = Lift*sin(phi) - Drag*cos(phi) ; % tangential
        
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
        f = (blades/2)*(R-blade_data(k))/(blade_data(k)*sin(abs(phi)));
        F= 2*acos(exp(-f))/pi;
        
        
        % We add this if statement otherwise the last element if NaN
        % (F = 0 !! )
        if k==N_element
            Wz(b,k,i) = 0 ;
            Wy(b,k,i) = 0 ;
        else
            Wz(b,k,i) = - blades*Lift*cos(phi)/(4*pi*rho*blade_data(k)*F*(sqrt(V0y^2+(V0z+u_turb+fg*Wz(b,k,i-1))^2))) ;
            Wy(b,k,i) = - blades*Lift*sin(phi)/(4*pi*rho*blade_data(k)*F*(sqrt(V0y^2+(V0z+u_turb+fg*Wz(b,k,i-1))^2))) ;
        end
        
        dm_edge(k,i) = blade_data(k)*py(i,k,b) ;
        dm_flap(k,i)= blade_data(k)*pz(i,k,b);
        dP(k) = omega0*dm_edge(k) ;
        
    end
    pz(i,N_element,b) = 0 ;
    py(i,N_element,b) = 0 ;
    dm_edge(N_element,N) = 0 ;
    dm_flap(N_element,N)=0;
    dP(N_element) = 0 ;
    M_edge(:,i)=trapz(blade_data(:,1), real(dm_edge(:,i))) ;
    M_flap(:,i)=trapz(blade_data(:,1), real(dm_flap(:,i)));
    thrust_blade(b) = trapz(blade_data(:,1), real(pz(i,:,b)));
end

if N_blade==3
    %blade1
    GF11(i)=trapz(blade_data(:,1), py(i,:,1)'.*uy_1f)+trapz(blade_data(:,1),pz(i,:,1)'.*uz_1f);
    GF12(i)=trapz(blade_data(:,1), py(i,:,1)'.*uy_1e)+trapz(blade_data(:,1),pz(i,:,1)'.*uz_1e);
    GF13(i)=trapz(blade_data(:,1), py(i,:,1)'.*uy_2f)+trapz(blade_data(:,1),pz(i,:,1)'.*uz_2f);
    %blade2
    GF21(i)=trapz(blade_data(:,1),py(i,:,2)'.*uy_1f)+trapz(blade_data(:,1),pz(i,:,2)'.*uz_1f);
    GF22(i)=trapz(blade_data(:,1),py(i,:,2)'.*uy_1e)+trapz(blade_data(:,1),pz(i,:,2)'.*uz_1e);
    GF23(i)=trapz(blade_data(:,1),py(i,:,2)'.*uy_2f)+trapz(blade_data(:,1),pz(i,:,2)'.*uz_2f);
    %blade3
    GF31(i)=trapz(blade_data(:,1),py(i,:,3)'.*uy_1f)+trapz(blade_data(:,1),pz(i,:,3)'.*uz_1f);
    GF32(i)=trapz(blade_data(:,1),py(i,:,3)'.*uy_1e)+trapz(blade_data(:,1),pz(i,:,3)'.*uz_1e);
    GF33(i)=trapz(blade_data(:,1),py(i,:,3)'.*uy_2f)+trapz(blade_data(:,1),pz(i,:,3)'.*uz_2f);
    % tower
    GFtower(i)= sum(thrust_blade); %sum of the thrust and the load of the wind on the tower ?
    GF(:,i)=[GFtower(i);GF11(i);GF12(i);GF13(i);GF21(i);GF22(i);GF23(i);GF31(i);GF32(i);GF33(i)];
end

if N_blade==1
    %blade1
    GF11(i)=trapz(blade_data(:,1),py(i,:,1)'.*uy_1f)+trapz(blade_data(:,1),pz(i,:,1)'.*uz_1f);
    GF12(i)=trapz(blade_data(:,1),py(i,:,1)'.*uy_1e)+trapz(blade_data(:,1),pz(i,:,1)'.*uz_1e);
    GF13(i)=trapz(blade_data(:,1),py(i,:,1)'.*uy_2f)+trapz(blade_data(:,1),pz(i,:,1)'.*uz_2f);
    
    GF(:,i)=[GF11(i);GF12(i);GF13(i)];
end
end
