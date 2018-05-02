%% Compute the relative wind velocity in Coordinate system 4 in y and z direction %%

function [Vrel_y, Vrel_z] = velocity_compute_turb(u_turb, b, k, Wy, Wz, Theta_wing1, Theta_wing2, Theta_wing3, Uy_dot, Uz_dot)
         
global a_12 a_21 a_34 blade_data H Ls omega0 Theta_cone V_0
        
        rt = [H 0 0] ; 
        rs = a_21*[0 0 -Ls]' ;
        
        if b==1 % blade 1  
            
            a_23_1 = [cos(Theta_wing1) sin(Theta_wing1) 0 ; 
                       -sin(Theta_wing1) cos(Theta_wing1) 0 ;
                       0 0 1] ;
            a_14_1 = a_34*a_23_1*a_12 ;
            a_41_1 = a_14_1' ;
            
            V0_4 = a_14_1*[0 0 V_0+u_turb]' ;
            Vrel_y = V0_4(2) + Wy - omega0*blade_data(k)*cos(Theta_cone) - Uy_dot ;
            Vrel_z = V0_4(3) + Wz - Uz_dot;
            
        elseif b==2 % blade 2
            a_23_2 = [cos(Theta_wing2) sin(Theta_wing2) 0 ; 
                      -sin(Theta_wing2) cos(Theta_wing2) 0 ;
                      0 0 1] ;
            a_14_2 = a_34*a_23_2*a_12 ;
            a_41_2 = a_14_2' ;            
         
            V0_4 = a_14_2*[0 0 V_0+u_turb]' ;
            Vrel_y = V0_4(2) + Wy - omega0*blade_data(k)*cos(Theta_cone) - Uy_dot ;
            Vrel_z = V0_4(3) + Wz - Uz_dot;
            
        else % blade 3 
            a_23_3 = [cos(Theta_wing3) sin(Theta_wing3) 0 ; 
                      -sin(Theta_wing3) cos(Theta_wing3) 0 ;
                      0 0 1] ;
            a_14_3 = a_34*a_23_3*a_12 ;
            a_41_3 = a_14_3' ; 
            
            V0_4 = a_14_3*[0 0 V_0+u_turb]' ;
            Vrel_y = V0_4(2) + Wy - omega0*blade_data(k)*cos(Theta_cone)- Uy_dot ;
            Vrel_z = V0_4(3) + Wz - Uz_dot ;
             
        end
    
end