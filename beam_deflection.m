function [uy, uz] = beam_deflection(py,pz,pitch)

global EI1global EI2global betaRadglobal rglobal 

    Ty = zeros(1,length(rglobal));
    Tz = zeros(1,length(rglobal));
    My = zeros(1,length(rglobal));
    Mz = zeros(1,length(rglobal));

    for i = 1:(length(rglobal)-1)
        j = length(rglobal)-i+1;
        Ty(j-1) = Ty(j) + 0.5*(py(j-1)+py(j))*(rglobal(j)-rglobal(j-1));
        Tz(j-1) = Tz(j) + 0.5*(pz(j-1)+pz(j))*(rglobal(j)-rglobal(j-1));
        My(j-1) = My(j) - Tz(j)*(rglobal(j)-rglobal(j-1)) - (pz(j-1)/6+pz(j)/3)*(rglobal(j)-rglobal(j-1))^2;
        Mz(j-1) = Mz(j) + Ty(j)*(rglobal(j)-rglobal(j-1)) + (py(j-1)/6+py(j)/3)*(rglobal(j)-rglobal(j-1))^2;
    end
    
    M1 = My.*cos(betaRadglobal + pitch)-Mz.*sin(betaRadglobal + pitch);
    M2 = My.*sin(betaRadglobal + pitch)+Mz.*cos(betaRadglobal + pitch);
    kappa1 = M1./EI1global;
    kappa2 = M2./EI2global;
    kappaz = -kappa1.*sin(betaRadglobal + pitch)+kappa2.*cos(betaRadglobal + pitch);
    kappay = kappa1.*cos(betaRadglobal + pitch)+kappa2.*sin(betaRadglobal + pitch);
    
    thetay = zeros(1,length(rglobal));
    thetaz = zeros(1,length(rglobal));
    uy = zeros(1,length(rglobal));
    uz = zeros(1,length(rglobal));

    for i = 1:(length(rglobal)-1)
        thetay(i+1) = thetay(i) + 0.5*(kappay(i+1)+kappay(i))*(rglobal(i+1)-rglobal(i));
        thetaz(i+1) = thetaz(i) + 0.5*(kappaz(i+1)+kappaz(i))*(rglobal(i+1)-rglobal(i));
        uy(i+1) = uy(i) + thetaz(i)*(rglobal(i+1)-rglobal(i)) + (kappaz(i+1)/6+kappaz(i)/3)*(rglobal(i+1)-rglobal(i))^2;
        uz(i+1) = uz(i) - thetay(i)*(rglobal(i+1)-rglobal(i)) - (kappay(i+1)/6+kappay(i)/3)*(rglobal(i+1)-rglobal(i))^2;
    end
end