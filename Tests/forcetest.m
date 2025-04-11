%create a routine that computes forces on each atom using a cutoff rc

%load data
[x,y,z] = loadData("fcc100a"+"256"+".txt");

nat = numel(x);
rc = 3;

%compute nndist and nnnbr

[nlabels, nn] = neigh(x,y,z,rc);

[fx,fy,fz] = forces(x,y,z,nn,nlabels);


