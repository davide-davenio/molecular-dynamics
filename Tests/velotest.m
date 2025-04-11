%test out assigning velocities
% compute temperature from assigned velocities
seed = 170901;
rng(seed)

mass = 108*(1.66*10.)^-37 * (1/16);  
nat = 256;
kb = 1/11603;
tkelv = 300;


[vx,vy,vz] = assignv(nat, tkelv);

kinetik = 0;
for i=1:nat
    kinetik = kinetik + ((vx(i).^2)+(vy(i).^2)+(vz(i).^2));
end

kinetik = kinetik*mass/2;
tresult = kinetik /(3/2 * nat *kb);

