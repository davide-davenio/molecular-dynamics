function dx = dist1(x1,x2, Lx)
    dx = x1-x2;
    dx = dx-Lx*round(dx/Lx);
end