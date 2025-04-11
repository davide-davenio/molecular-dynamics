% EXERCIZE 1 - PART B.1: for 1200K, save positions step by step 

addpath 'C:\Users\david\Desktop\Montalens\subroutines'
%parameters
datafile = "fcc100a"+"256"+".txt";
mass = 108*1.66*10^(-27) * (1/16);  
kb = 1/11603;
tkelv = 1200;
rc = 4.5;
rp= 4.2;
MDsteps = 12500;
dt_step = 2*10^(-15);
thermalization = 3*10^(-12);  %amount of seconds of simulation to discard
eps = 0.345;
sig = 2.644;
dt = dt_step;

tic;

%assign seed 
seed = 42;
rng(seed)



%load data
[x,y,z] = loadData(datafile);
nat = numel(x);
   
[vx,vy,vz] = assignv(nat, tkelv,mass, kb);
[vx, vy, vz] = fixVelocities(vx, vy, vz, tkelv,nat, mass, kb);

%compute initial forces
[labels, nn]= neigh(x,y,z,rc);
[fx,fy,fz] = forcesPoly(x,y,z,nn,labels,nat,rp, rc, eps,sig);

outputname = sprintf('Ex1_A_3movie.txt');
foutput = fopen(outputname,'w');


%velocity verlet algo
for j=1:MDsteps+floor(thermalization/dt)
    %position update
    for i=1:nat
        x(i)=x(i)+vx(i)*dt+0.5*(dt^2)*(fx(i)/mass);
        y(i)=y(i)+vy(i)*dt+0.5*(dt^2)*(fy(i)/mass);
        z(i)=z(i)+vz(i)*dt+0.5*(dt^2)*(fz(i)/mass);
    end
    
    [labels, nn,] = neigh(x,y,z,rc);
    [fx_next,fy_next,fz_next] = forcesPoly(x,y,z,nn,labels,nat,rp, rc, eps,sig);
    
    % velocity update
    for i=1:nat
        vx(i) = vx(i) + ((fx(i)+fx_next(i))/(2*mass))*dt;
        vy(i) = vy(i) + ((fy(i)+fy_next(i))/(2*mass))*dt;
        vz(i) = vz(i) + ((fz(i)+fz_next(i))/(2*mass))*dt;
    end

    if j > floor(thermalization/dt)
         
         fname = sprintf('Snap%03i.xyz',j);
         writeXYZ(fname,x,y,z, nat);
    end

    %update forces
    fx = fx_next;
    fy = fy_next;
    fz = fz_next;
end


fclose(foutput);
elapsedtime = toc;
fprintf("elapsed time %g seconds", elapsedtime);