function [cmx, cmy, cmz] = CoM(x,y,z)
%compute center of mass of lattice
    cmx=sum(x);
    cmy=sum(y);
    cmz=sum(z);
end



