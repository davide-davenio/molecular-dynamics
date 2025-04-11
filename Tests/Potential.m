% cmpute total pot energy with lenard jones of the whole cell
% E = 1/2 * sumij(phi(rij)
% phi = 4 eps * [(sig/r)^12-(sig/r)^6]
% sig = 2.644 A
% eps = 0.345 eV
%plot E/N vs N


function [d] = dist(Ax,Ay,Az,Bx,By,Bz)
        d = sqrt((Ax-Bx)*(Ax-Bx) + (Ay-By)*(Ay-By) + (Az-Bz)*(Az-Bz));
end

function [pot] = LJpotential(r)
        sig = 2.644;
        eps = 0.345;
        pot = 4 * eps * ((sig/r).^12 - (sig/r).^6);
end

function [Etot] = energy(x, y, z)
    Etot = 0;
    for i=1:numel(x)
        for j=1:numel(x)
            if j ~= i
            Etot = Etot + LJpotential(dist(x(i), y(i), z(i), x(j), y(j), z(j)));
            end
        end
    end
    Etot = 0.5 * Etot;
end


function [x, y, z] = loadData(dataFile) 
    data = load(dataFile);
    x = data(:, 1);
    y = data(:, 2);
    z = data(:, 3);
end 


%compute all energies
Num = [108, 256, 500, 864, 1372, 2048];
Etot = zeros(size(Num));

for i=1:numel(Num)
    
    [x,y,z] = loadData("fcc100a"+Num(i)+".txt");
    Etot(i) = energy(x,y,z)/Num(i);

end
