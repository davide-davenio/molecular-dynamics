% produces xyz configuration files with
% potential enrgy minimized through steepest descent

addpath 'C:\Users\david\Desktop\Montalens\subroutines'
%parameters

datafile = "fcc100a256.txt";
outputname = "fcc100a256_minimized.txt";
params1B2_100;

%stats
output_stats = "fcc100a256_minstats.txt";
foutput_stats = fopen(output_stats,'w');



%produce minimized 
[x,y,z] = loadData(datafile);
nat = numel(x);

flag = 0;
maxStep = 10000;
for it=1:maxStep

    [x,y,z,flag, maxf]=steepestDescentStep(x,y,z,rc,rp, eps, sig, nat, Csteep, ftoll,Lx,Ly,Lz);
    
    %save positions
    %fname = sprintf('Snap%03i.xyz',it);
    %writeXYZ(fname,x,y,z, nat);
    
    %save stats 
    [labels, nn] = neigh(x,y,z,rc,Lx,Ly,Lz);
    [Epot]= energyPoly(x, y, z, labels, nn, sig, eps, rp, rc,Lx,Ly,Lz);
    fprintf(foutput_stats,"%i %f %f\n",it, Epot, maxf);

    if flag == 1
        break;
    end
end
writeXYZ(outputname,x,y,z, nat);

