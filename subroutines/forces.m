function [fx,fy,fz] = forces(x,y,z,nn,labels,nat,eps,sig,Lx,Ly,Lz)
%computes forces on all atoms

    
    %make force matrix
   
    fx = zeros(nat,1);
    fy = zeros(nat,1);
    fz = zeros(nat,1);

    %for each atom compute its coefficients 
    for i=1:nat

        for j=1:nn(i)
          
            if labels(i,j) ~= 0

                d = dist(x(i), y(i),z(i), x(labels(i,j)),y(labels(i,j)),z(labels(i,j)),Lx,Ly,Lz);

                coef1 = 24*(sig^6)*eps*(1/(d^8));
                coef2 = (2*(sig^6)/(d^6)) - 1;
                
                dx = x(labels(i,j))-x(i);
                dy = y(labels(i,j))-y(i);
                dz = z(labels(i,j))-z(i);
                
                fx(i) = fx(i) - dx*coef1*coef2;
                fy(i) = fy(i) - dy*coef1*coef2;
                fz(i) = fz(i) - dz*coef1*coef2;

            end
            
        end
    
    end

end