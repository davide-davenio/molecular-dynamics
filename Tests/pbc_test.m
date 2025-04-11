

Lx = 16.6416;
Ly = Lx;
Lz = Lx;

rc = 4.5;
rp = 4.2;
eps = 0.345;
sig = 2.644;

initconf = "fcc100a"+"256"+".txt";
[x,y,z] = loadData(initconf);
nat = numel(x);


[labels, nn]= neigh(x,y,z,rc,Lx,Ly,Lz);
[fx,fy,fz] = forces(x,y,z,nn,labels,nat,eps,sig,Lx,Ly,Lz);