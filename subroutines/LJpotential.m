% computes LJpotential 
function [pot] = LJpotential(r, sig, eps)
    % computes LJpotential 
    pot = 4 * eps * ((sig/r)^12 - (sig/r)^6);
end


