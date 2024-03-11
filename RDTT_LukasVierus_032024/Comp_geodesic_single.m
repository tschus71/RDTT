function geodesic = Comp_geodesic_single(x,xi,step,type)
    geodesic=[];

    u_out=[x,xi];    
        
    u_out(abs(u_out)<10^(-10))=0;
    if dot(u_out(1:2),u_out(3:4))>0
       u=[x;-xi]; % -xi points inside

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

       xi_normed = 1/norm(xi) * xi; 
       nplus1    = [x+tau*xi; xi_normed];
       geodesic  = [geodesic,nplus1];

       % Compute last distance
       last_step = step*norm(geodesic(:,n-1)-geodesic(:,n+1))/norm(geodesic(:,n-1)-geodesic(:,n));
       % insert distance as last entry
       geodesic = [geodesic, [last_step;0;0;0]];  
       % delete point that is not needed
       geodesic(:,n)=[];    
    end
end
