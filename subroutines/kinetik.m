function [E, t] = kinetik(vx,vy,vz,kb,mass,nat)
    %   computes kinetik energy for a set of velocities
   

    E = 0;
    for i=1:nat
    E = E + ((vx(i)^2)+(vy(i)^2)+(vz(i)^2));
    end
    
    E = E*mass/2;
    t = E /(3/2 * nat *kb);
end