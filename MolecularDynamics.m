% Simulate using velocity verlet algorithm a crystal given a set of
% parameters and an initial configuration
% param file is a *.m file with parameters assigned as: mass = 49593;
% initial configuration file is a three column file (separated by spaces)
% that contains only xyz coordinates of all atoms of the simulation cell


addpath 'C:\Users\david\Desktop\Montalens\subroutines'

%input
initconf = "fcc111a336+1_minimized.txt";
outputname = 'Ex1_C_3_111.txt';
params1C3_111; %this line runs the params file and loads all variables into the environment


tic;



%load data
[x,y,z] = loadData(initconf);
nat = numel(x);


%assign seed 
rng(seed)




   
[vx, vy, vz] = assignv(nat, tkelv,mass, kb);
[vx, vy, vz] = fixVelocities(vx, vy, vz, tkelv,nat, mass, kb);

%compute initial forces
[labels, nn]= neigh(x,y,z,rc,Lx,Ly,Lz);
[fx,fy,fz] = forcesPoly(x,y,z,nn,labels,nat,rp, rc, eps,sig,Lx,Ly,Lz);


foutput = fopen(outputname,'w');
%velocity verlet algo
for j=1:MDsteps+floor(thermalization/dt)
    %position update
    for i=1:nat
        x(i)=x(i)+vx(i)*dt+0.5*(dt^2)*(fx(i)/mass);
        y(i)=y(i)+vy(i)*dt+0.5*(dt^2)*(fy(i)/mass);
        z(i)=z(i)+vz(i)*dt+0.5*(dt^2)*(fz(i)/mass);
    end
    
    [labels, nn] = neigh(x,y,z,rc,Lx,Ly,Lz);
    [fx_next,fy_next,fz_next] = forcesPoly(x,y,z,nn,labels,nat,rp,rc,eps,sig,Lx,Ly,Lz);
    
    % velocity update
    for i=1:nat
        vx(i) = vx(i) + ((fx(i)+fx_next(i))/(2*mass))*dt;
        vy(i) = vy(i) + ((fy(i)+fy_next(i))/(2*mass))*dt;
        vz(i) = vz(i) + ((fz(i)+fz_next(i))/(2*mass))*dt;
    end
             

    
    if j > floor(thermalization/dt)
    
    %to save configurations step by step
         fname = sprintf('Snap%03i.xyz',j);
         writeXYZ(fname,x,y,z, nat);
    

    % to save potential energy and temperature step by step
        [Ekin, tresult] = kinetik(vx,vy,vz,kb,mass,nat);
        [Epot]= energyPoly(x, y, z, labels, nn, sig, eps, rp, rc,Lx,Ly,Lz);
        Etot = Ekin + Epot;
        fprintf(foutput,"%i %f %f\n",j, Etot, tresult);
    %

    end



    %update forces
    fx = fx_next;
    fy = fy_next;
    fz = fz_next;
end


fclose(foutput);
elapsedtime = toc;
fprintf("elapsed time %g seconds", elapsedtime);