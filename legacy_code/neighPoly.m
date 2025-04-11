function[nnlabels,nnbr,nnlabels2, nnbr2] = neighPoly(vx, vy, vz,rpoly,rc, maxnn)
    % returns a vector with the number of neighbours (under a certain 
    % cutoff) and "second neighbours" between the first and second cutoff 
    % for each point and matrix that encodes the indices (of the original 
    % vector) of each neighbour of each point


    %define default value for maxnn
    if nargin < 6
        maxnn = 64;  % Default value for maxnn
    end
    
    derror = 0.01; %this is an arbitrary, scale dependant, value 
    nat = numel(vx);
    nnbr = zeros(nat,1);
    nnlabels = zeros (nat, maxnn);
    nnbr2 = zeros(nat,1);
    nnlabels2 = zeros (nat, maxnn);
    for i=1:nat
        for j=1:nat
            
            dtemp = dist(vx(i), vy(i), vz(i), vx(j),vy(j),vz(j));
            if dtemp < rpoly+derror && dtemp > derror
                nnbr(i) = nnbr(i) +1;
                nnlabels(i, nnbr(i)) = j;
            end
            if dtemp < rc+derror && dtemp > rpoly+derror
                nnbr2(i) = nnbr2(i) +1;
                nnlabels2(i, nnbr2(i)) = j;
            end
        end
    end
end