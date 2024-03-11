function v = RightSideODE (u, speed)
% u is the vector [gamma1(t), gamma2(t), gamma1'(t), gamma2'(t)]
% Soundspeed contains the speed of sound and its gradient

% If we are close enough to the boundary of the domain, we neglect refraction
% effects and approximate the geodesics by eucledian curves

if norm(u(1:2)) < 0.95
   help=speed;
   c=help(1);
   dc_dx=help(2);
   dc_dy=help(3);
else
    c=1;
    dc_dx=0;
    dc_dy=0;
end

v=zeros(4,1);
v(1)=u(3);
v(2)=u(4);

norm_xi          = power(u(3),2) + power(u(4),2);
product_xi_gradc = u(3)*dc_dx + u(4)*dc_dy;

v(3)= -(dc_dx*norm_xi - 2*u(3)*product_xi_gradc)/c;
v(4)= -(dc_dy*norm_xi - 2*u(4)*product_xi_gradc)/c;
end
