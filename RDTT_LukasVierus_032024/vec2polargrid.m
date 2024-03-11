function [field1_grid, field2_grid] = vec2polargrid(field, anz_r, anz_phi)
% Computes f on grid

field1_grid = zeros(anz_r, anz_phi);
field2_grid = zeros(anz_r, anz_phi);

rho = (1/anz_r:1/anz_r:1);
phi = (2*pi/anz_phi:2*pi/anz_phi:2*pi);

for i=1:anz_r
    for j=1:anz_phi
        x = rho(i)*[cos(phi(j));sin(phi(j))];
        fx = field(x);
        field1_grid(i,j)= fx(1);
        field2_grid(i,j)= fx(2);
    end
end

end

