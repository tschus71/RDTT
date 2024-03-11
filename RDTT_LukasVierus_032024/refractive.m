function ngradn = refractive (x)

speed = Soundspeed(x,0);

c     = speed(1);
dc_x1 = speed(2);
dc_x2 = speed(3);

%n, gradn
ngradn = [1./c; -dc_x1./c.^2; -dc_x2./c.^2];
end