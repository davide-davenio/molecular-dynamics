function [vx,vy,vz] = fixVelocities(vx_i, vy_i, vz_i, T,nat, mass, kb)
        
    [xavg, yavg, zavg] = deal(0);
    vsquared = 0;
    for i=1:nat
        xavg = xavg + vx_i(i);
        yavg = yavg + vy_i(i);
        zavg = zavg + vz_i(i);

       
    end
    xavg = xavg/nat;
    yavg = yavg/nat;
    zavg = zavg/nat;

    vx = (vx_i-xavg);
    vy = (vy_i-yavg);
    vz = (vz_i-zavg);

    for i=1:nat
         vsquared = vsquared + (vx(i)^2 + vy(i)^2 + vz(i)^2);
    end
    Tprime = (mass*vsquared)/(3*kb*nat);

    vx = vx*sqrt(T/Tprime);
    vy = vy*sqrt(T/Tprime);
    vz = vz*sqrt(T/Tprime);
end