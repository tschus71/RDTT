% Computes tau such that ||x+tau*xi||=1, function works in 2D and 3D

function time = ComputeLastStep(x,xi)
a=-dot(x,xi)/norm(xi)^2;
b= (norm(x)^2-1)/norm(xi)^2;
time= a+ sqrt(a^2-b);
