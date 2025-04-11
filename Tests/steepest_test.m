%testing steepest descent

addpath 'C:\Users\david\Desktop\Montalens\subroutines'
%parameters
datafile = "fcc100a"+"256"+".txt";
rc = 4.5;
rp= 4.2;
Csteep = 0.005;
ftoll = 1.e-3;
eps = 0.345;
sig = 2.644; 
%load
[x,y,z] = loadData(datafile);
nat = numel(x);

%test steepest descent

maxit = 10000;
flag = 0;

maxf = zeros(maxit, 1);
Epot = zeros(maxit, 1);

for it= 1:maxit
    [x,y,z,flag, maxf(it)]=steepestDescentStep(x,y,z,rc,rp, eps, sig, nat, Csteep, ftoll);

    [labels, nn]= neigh(x,y,z,rc);
    Epot(it) = energyPoly(x, y, z, labels, nn, sig, eps, rp, rc);
    

        
    if flag == 1
        break
    end
end


%plot and save plots
figure;
plot(Epot(1:it));
title('potential energy vs iterations');
xlabel('iterations');
ylabel('E_{pot} [eV]');
grid on;
exportgraphics(gcf, 'Epot_vs_iterations_1B2.pdf', 'ContentType', 'vector');  

figure;
plot(maxf(1:it));
title('max force vs iterations');
xlabel('iterations');
ylabel('F_{max} [eV/Å]');
grid on;
exportgraphics(gcf, 'MaxF_vs_iterations_1B2.pdf', 'ContentType', 'vector');  

