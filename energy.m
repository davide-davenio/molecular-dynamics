%computes energy of the whole crystal given the neighbour matrix
function [Epot] = energy(x, y, z, nlabels, nn, sig, eps,Lx,Ly,Lz)

    Epot = 0;
    nat = numel(x);
    for i=1:nat
        for j=1:nn(i)
            label = nlabels(i,j);
            Epot = Epot + LJpotential(dist(x(i), y(i), z(i), x(label), y(label), z(label),Lx,Ly,Lz), sig, eps);
          
        end
    end
    Epot = 0.5 * Epot;
end