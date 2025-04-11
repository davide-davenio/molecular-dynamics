%compute lattice potential with lennard jones with a nearest neighbour
%approximation

function [x, y, z] = loadData(dataFile) 
    data = load(dataFile);
    x = data(:, 1);
    y = data(:, 2);
    z = data(:, 3);
end 

function [d] = dist(Ax,Ay,Az,Bx,By,Bz)
        d = sqrt((Ax-Bx)*(Ax-Bx) + (Ay-By)*(Ay-By) + (Az-Bz)*(Az-Bz));
end

function [d] = nearNeigh(point, vx, vy, vz)
    d = realmax;
    for i=1:numel(vx)
        if i ~= point
            dtemp = dist(vx(point),vy(point),vz(point), vx(i), vy(i), vz(i));
            if dtemp < d
                d = dtemp;
            end
        end
    end
end 

function[nnlabels,nnbr] = neigh(vx, vy, vz, cutoff, maxnnn)
    derror = 0.1; %this is an arbitrary value that has to be changed based on the scale of the problem
    nat = numel(vx);
    nnbr = zeros(nat,1);
    nnlabels = zeros (nat, maxnnn);
    for i=1:nat
        for j=1:nat
        
            dtemp = dist(vx(i), vy(i), vz(i), vx(j),vy(j),vz(j));
            if dtemp < cutoff+derror && dtemp ~= 0
                nnbr(i) = nnbr(i) +1;
                nnlabels(i, nnbr(i)) = j;
            end
        end
    end
end

function [pot] = LJpotential(r)
        sig = 2.644;
        eps = 0.345;
        pot = 4 * eps * ((sig/r).^12 - (sig/r).^6);
end

function [Etot] = energy(x, y, z, nlabels)
    Etot = 0;
    nat = numel(x);
    maxnnn = 64;
    for i=1:nat
        for j=1:maxnnn
            label = nlabels(i,j);
            if label ~= 0 
                Etot = Etot + LJpotential(dist(x(i), y(i), z(i), x(label), y(label), z(label)));
            end
        end
    end
    Etot = 0.5 * Etot;
end


%load data
[x,y,z] = loadData("fcc100a"+"108"+".txt");

nat = numel(x);
maxnn = 64;
rc = 3;

%compute nndist and nnnbr

[nlabels, nn] = neigh(x,y,z,rc,maxnn);

Etot = energy(x,y,z,nlabels);