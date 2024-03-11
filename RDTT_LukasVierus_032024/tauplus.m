function tau = tauplus(x,xi)
tau = -dot(x,xi)+sqrt(dot(x,xi)^2+1-dot(x,x));
end

