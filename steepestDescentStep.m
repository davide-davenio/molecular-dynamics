function [x,y,z,flag, maxf]=steepestDescentStep(x,y,z,rc,rp, eps, sig, nat, Csteep, ftoll,Lx,Ly,Lz) 
   
    flag = 0;
    [labels, nn]= neigh(x,y,z,rc,Lx,Ly,Lz);
    [fx,fy,fz] = forcesPoly(x,y,z,nn,labels,nat,rp, rc, eps,sig,Lx,Ly,Lz);
        
    maxf = max(sqrt(sum(fx.^2)+sum(fy.^2)+sum(fz.^2)));
        
    if(maxf < ftoll)
        flag = 1;
    
    else
        for i=1:nat
            x(i) =  x(i) + Csteep * fx(i);
            y(i) =  y(i) + Csteep * fy(i);
            z(i) =  z(i) + Csteep * fz(i);
        end
    end
end


