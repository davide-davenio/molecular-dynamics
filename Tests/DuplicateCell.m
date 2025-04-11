% duplicate the lattice in any direction

function [x, y, z] = loadData(dataFile) 
    data = load(dataFile);
    x = data(:, 1);
    y = data(:, 2);
    z = data(:, 3);
end 

function [newv] = duplicateCrystal(v)
 
    vdist = abs(diff(v));
    vdist(vdist == 0) = [];
    side = max(vdist);
    a = min(vdist);
    newv = v+side+a;
 
end 
[x,y,z] = loadData("fcc100a108.txt");

xnew = duplicateCrystal(x);

x = [x; xnew];
y = [y;y];
z = [z;z];




