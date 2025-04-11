%repeat cutoff computation using newly made library

%load data
[x,y,z] = loadData("fcc100a"+"108"+".txt");

nat = numel(x);
rc = 3;

%compute nndist and nnnbr

[nlabels, nn] = neigh(x,y,z,rc);

Etot = energy(x,y,z,nlabels);

