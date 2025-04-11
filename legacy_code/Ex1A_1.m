% EXERCIZE 1 - PART A.1: using sharp cutoff approach find largest time step
% possible while keeping energy fluctuations under 10^-5 for T_i = 20K

addpath 'C:\Users\david\Desktop\Montalens\subroutines'
%parameters
datafile = "fcc100a"+"256"+".txt";
mass = 108*1.66*10^(-27) * (1/16);  
kb = 1/11603;
tkelv = 20;
rc = 4.5;
MDsteps = 5000;
dt_step = 5*10^(-15);
thermalization = 3*10^(-12);  %amount of seconds of simulation to discard
eps = 0.345;
sig = 2.644;


tic;
%load data
[x,y,z] = loadData(datafile);
nat = numel(x);

%assign seed 
seed = 170901;
rng(seed)


%looping over different timesteps to find the one that best satisfies solution
timesteps = 6;
dt = dt_step;
for k=1:timesteps
    
    %save timestep in a separate file
    if k==1
        foutput = fopen("timestepslist1A1.txt", "w");
        fprintf(foutput, "%g\n", dt);
        fclose(foutput);
    else
        foutput = fopen("timestepslist1A1.txt", "a");
        fprintf(foutput, "%g\n", dt);
        fclose(foutput);
    end

    %load data
    [x,y,z] = loadData(datafile);
    nat = numel(x);
    
    [vx,vy,vz] = assignv(nat, tkelv,mass, kb);

    %compute initial forces
    [nlabels, nn] = neigh(x,y,z,rc);
    [fx,fy,fz] = forces(x,y,z,nn,nlabels,nat,eps,sig);

    outputname = sprintf('Ex1_A_1%05i.txt', k);
    foutput = fopen(outputname,'w');


    %velocity verlet algo
    for j=1:MDsteps+floor(thermalization/dt)
        %position update
        for i=1:nat
            x(i)=x(i)+vx(i)*dt+0.5*(dt^2)*(fx(i)/mass);
            y(i)=y(i)+vy(i)*dt+0.5*(dt^2)*(fy(i)/mass);
            z(i)=z(i)+vz(i)*dt+0.5*(dt^2)*(fz(i)/mass);
        end
    
        [nlabels, nn] = neigh(x,y,z,rc);
        [fx_next,fy_next,fz_next] = forces(x,y,z,nn,nlabels,nat,eps,sig);
    
        %velocity update
        for i=1:nat
            vx(i) = vx(i) + ((fx(i)+fx_next(i))/(2*mass))*dt;
            vy(i) = vy(i) + ((fy(i)+fy_next(i))/(2*mass))*dt;
            vz(i) = vz(i) + ((fz(i)+fz_next(i))/(2*mass))*dt;
        end

        if j > floor(thermalization/dt)
            [Ekin, tresult] = kinetik(vx,vy,vz,kb,mass,nat);
            [Epot]= energy(x,y,z,nlabels,nn, sig, eps);
            Etot = Ekin + Epot;
            fprintf(foutput,"%i %f %f\n",j, Etot, tresult);
        end
        %update forces
        fx = fx_next;
        fy = fy_next;
        fz = fz_next;
    end
    fclose(foutput);
    dt = dt+dt_step;
    fprintf("iteration: %i done\n", k);
end
elapsedtime = toc;
fprintf("elapsed time %g seconds", elapsedtime);