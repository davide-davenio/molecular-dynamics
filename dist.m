%computes distance in 3D
function [d] = dist(Ax,Ay,Az,Bx,By,Bz,Lx,Ly,Lz)
        d = sqrt((dist1(Ax,Bx,Lx)*dist1(Ax,Bx,Lx)) + (dist1(Ay,By,Ly)*dist1(Ay,By,Ly)) + (dist1(Az,Bz,Lz)*dist1(Az,Bz,Lz)));
end