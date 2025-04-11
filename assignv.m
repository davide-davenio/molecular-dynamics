function[vx,vy,vz] = assignv(nat, tkelvin, mass, kb)
    
    const = sqrt(3*kb*tkelvin/mass); 

    vx = zeros(nat, 1);
    vy = zeros(nat, 1);
    vz = zeros(nat,1);
    for i=1:nat
        vx(i) = (2 * rand() - 1)*const;
        vy(i) = (2 * rand() - 1)*const;
        vz(i) = (2 * rand() - 1)*const;
    end
end