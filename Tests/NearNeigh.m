% find the nearest neighbour distance and second neighbour distance 
% for each atom compute number of neighbours and save in a vector



function [d] = computeDistance(Ax,Ay,Az,Bx,By,Bz)
        d = sqrt((Ax-Bx)*(Ax-Bx) + (Ay-By)*(Ay-By) + (Az-Bz)*(Az-Bz));
end

function [d] = nearNeigh(point, vx, vy, vz)
    d = realmax;
    for i=1:numel(vx)
        if i ~= point
            dtemp = computeDistance(vx(point),vy(point),vz(point), vx(i), vy(i), vz(i));
            if dtemp < d
                d = dtemp;
            end
        end
    end
end 

function [d] = secondNeigh(point, vx,vy,vz, nnDistance)
    derror = 0.1; %this is and arbitrary parameter based on the scale of the problem
    d = realmax;
    for i=1:numel(vx)
        if i ~= point
           dtemp = computeDistance(vx(point),vy(point),vz(point), vx(i), vy(i), vz(i));
            if dtemp < d && dtemp > nnDistance + derror
                d = dtemp;
            end 
        end
    end
end

function[number] = neighNumber(point, vx, vy, vz, nDistance)
    derror = 0.1; %this is an arbitrary value that has to be changed based on the scale of the problem
    number = 0;
    for i=1:numel(vx)
        dtemp = computeDistance(vx(point), vy(point), vz(point), vx(i),vy(i),vz(i));
        if dtemp > nDistance-derror && dtemp < nDistance+derror
            number = number +1;
        end
    end
end

%load data
load fcc100a108.txt
x = fcc100a108(:, 1);
y = fcc100a108(:, 2);
z = fcc100a108(:, 3);

%testing funtions
nnd = nearNeigh(1, x, y, z);

snd = secondNeigh(1, x, y, z, nnd);

nnn = zeros(numel(x),1);
for i=1:numel(x)
    nnn(i) = neighNumber(i,x,y,z,nnd); 
end

