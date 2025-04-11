% returns a vector with the number of neighbours (under a certain cutoff)
% for each point and matrix that encodes the indices (of the original 
% vector) of each neighbour of each point

function[nnlabels,nnbr] = neigh(vx, vy, vz, cutoff,Lx,Ly,Lz)
  
    maxnn = 64;  % Default value for maxnn
    derror = 0.01; %this is an arbitrary, scale dependant, value 
    
    
    nat = numel(vx);
    nnbr = zeros(nat,1);
    nnlabels = zeros (nat, maxnn);
    for i=1:nat
        for j=1:nat
        
            dtemp = dist(vx(i), vy(i), vz(i), vx(j),vy(j),vz(j),Lx,Ly,Lz);
            if dtemp < cutoff+derror && dtemp ~= 0
                nnbr(i) = nnbr(i) +1;
                nnlabels(i, nnbr(i)) = j;
            end
        end
    end
end