% testing the computation of energy with the polynomial junction

[x,y,z] = loadData("fcc100a2048.txt");

rc = 4.5;
rpoly = 4.2;
sig = 2.644;
eps = 0.345;

[nlabels,nn,nlabels2,nn2] = neighPoly(x, y, z, rpoly, rc);
Epot = energyPoly(x, y, z, nlabels, nn, nlabels2, nn2, sig, eps, rpoly, rc);
Epot_base = energy(x, y, z, nlabels, nn, sig, eps);

