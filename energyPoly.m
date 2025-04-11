%computes energy of the whole crystal given the neighbour matrix
function [Epot] = energyPoly(x, y, z, nlabels, nn, sig, eps, rpoly, rc,Lx,Ly,Lz)

    Epot = 0;
    nat = numel(x);
    for i=1:nat

        for j=1:nn(i)
            label = nlabels(i,j);
            dtemp = dist(x(i), y(i), z(i), x(label), y(label), z(label),Lx,Ly,Lz);
            if dtemp < rpoly
                Epot = Epot + LJpotential(dtemp, sig, eps); 
            else
                Epot = Epot + polyPotential(dtemp, rpoly, rc,eps, sig);
            end
        end

    end
    Epot = 0.5 * Epot;
end