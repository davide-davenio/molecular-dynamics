function [fx,fy,fz] = forcesPoly(x,y,z,nn,labels,nat,rpoly, rc, eps,sig,Lx,Ly,Lz)
%computes forces on all atoms with polynomial junction

    fx = zeros(nat,1);
    fy = zeros(nat,1);
    fz = zeros(nat,1);

    %for each atom compute its coefficients 
    for i=1:nat

        for j=1:nn(i)

            label = labels(i,j);

            if label ~= 0

                d = dist(x(i), y(i),z(i), x(label),y(label),z(label),Lx,Ly,Lz);
                if d < rpoly
                    coef1 = 24*(sig^6)*eps*(1/(d^8));
                    coef2 = (2*(sig^6)/(d^6)) - 1;
                    
                    dx = dist1(x(label),x(i), Lx);
                    dy = dist1(y(label),y(i), Ly);
                    dz = dist1(z(label),z(i), Lz);
                
                    fx(i) = fx(i) - dx*coef1*coef2;
                    fy(i) = fy(i) - dy*coef1*coef2;
                    fz(i) = fz(i) - dz*coef1*coef2;
                else
                    cpoly = polyDerivative(d, rpoly, rc,eps, sig);
                
                    dx = dist1(x(label),x(i), Lx);
                    dy = dist1(y(label),y(i), Ly);
                    dz = dist1(z(label),z(i), Lz);
                    
                    fx(i) = fx(i) + dx*cpoly*(1/d);
                    fy(i) = fy(i) + dy*cpoly*(1/d);
                    fz(i) = fz(i) + dz*cpoly*(1/d);
            
                end
            end
        end
    end
end