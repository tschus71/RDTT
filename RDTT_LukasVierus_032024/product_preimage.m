function value = product_preimage(arrayA1, arrayA2, arrayB1, arrayB2) % L^2(Omega)-norm for vector fields
integrand = arrayA1 .* arrayB1 + arrayA2 .* arrayB2; % =(f*g)(r,phi)

anz_rad = size(integrand,1);
anz_ang = size(integrand,2);

h_rad = 1/anz_rad;
h_ang = 2*pi/anz_ang;

r   = (h_rad:h_rad:1);

rg = r'.*integrand;
% trapezoidal rule in 2D; value at (0,0) not necessary since functional determinant r=0
r_integral = h_rad*(sum(rg)-0.5*integrand(anz_rad,:)); 

value = abs(sqrt(h_ang * sum(r_integral))); % norm
