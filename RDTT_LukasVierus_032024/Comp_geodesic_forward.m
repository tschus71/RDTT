% This programm computes one geodesic from tau = 0 until tau = tau_+
function geodesic = Comp_geodesic_forward(x,xi,step,type)
    geodesic=[];

    u=[x;xi]; %u_in
    u(abs(u)<10^(-10))=0;
    geodesic=[geodesic,u];
           
    u = RungeKuttaODE(u,step,Soundspeed(u(1:2),type));
    u(abs(u)<10^(-10))=0;
    geodesic=[geodesic,u];

    while (norm(u(1:2))<1)
          u = RungeKuttaODE(u,step,Soundspeed(u(1:2),type));
          u(abs(u)<10^(-10))=0;
          geodesic=[geodesic,u];
    end
    n=size(geodesic,2);

    %Last step to arrive on the boundary of M
    x   = geodesic(1:2,n-1);
    xi  = geodesic(3:4,n-1);
    tau = ComputeLastStep(x,xi);
    tau = real(tau);

    xi_normed     = 1/norm(xi) * xi; 
    nplus1        = [x+tau*xi; xi_normed];
    geodesic      = [geodesic,nplus1];
    geodesic(:,n) = [];
    geodesic      = [geodesic,[tau;0;0;0]];
end
