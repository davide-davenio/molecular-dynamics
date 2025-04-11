%create a routine that computes forces on each atom using a cutoff rc and
%poly junction

%load data
[x,y,z] = loadData("fcc100a"+"256"+".txt");

nat = numel(x);
rc = 4.5;
rpoly=4.2;
sig = 2.644;
eps = 0.345;

%compute forces

[labels,nn,labels2, nn2] = neighPoly(x, y, z,rpoly,rc);

[fx_base,fy_base,fz_base] = forces(x,y,z,nn,labels,nat, eps, sig);

[fx,fy,fz] = forcesPoly(x,y,z,nn,labels,nn2,labels2,nat,rpoly, rc, eps,sig);


