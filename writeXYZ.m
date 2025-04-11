function writeXYZ(fname,x,y,z, nat)
    %saves positions for each timestep
  
    filesnap=fopen(fname,'w');
    for i=1:nat
        fprintf(filesnap, '%g %g %g\n', x(i),y(i),z(i) );
    end
    fclose(filesnap);
end