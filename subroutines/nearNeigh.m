% finds nearest neighbour distance for a point
function [d] = nearNeigh(point, vx, vy, vz,Lx,Ly,Lz)
    d = realmax;
    for i=1:numel(vx)
        if i ~= point
            dtemp = dist(vx(point),vy(point),vz(point), vx(i), vy(i), vz(i),Lx,Ly,Lz);
            if dtemp < d
                d = dtemp;
            end
        end
    end
end 