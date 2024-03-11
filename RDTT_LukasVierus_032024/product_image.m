function value = product_image(boundary_values)
% This function computes the L^2-norm on \partial\Omega M

anz_phi = size(boundary_values,1);
anz_theta = size(boundary_values,2);

value = 4*pi^2/(anz_phi*anz_theta) * sum(sum(boundary_values.^2));
value = sqrt(value);
end

