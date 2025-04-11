function [x,y,z]=steepestDescent(x,y,z,rc,rp, eps, sig, nat, Csteep, ftoll,Lx,Ly,Lz) 
   
    maxStep = 10000;
    for it=1:maxStep % we don't want to run the loop forever if not converging
        
        [labels, nn]= neigh(x,y,z,rc,Lx,Ly,Lz);
        [fx,fy,fz] = forcesPoly(x,y,z,nn,labels,nat,rp, rc, eps,sig,Lx,Ly,Lz);
        
        maxf = max(sqrt(sum(fx.^2)+sum(fy.^2)+sum(fz.^2)));
    
        if(maxf < ftoll)
            break % converged! We are done and exit the loop
    
        else
            for i=1:nat
                x(i) =  x(i) + Csteep * fx(i);
                y(i) =  y(i) + Csteep * fy(i);
                z(i) =  z(i) + Csteep * fz(i);
            end
        end
    end

end
