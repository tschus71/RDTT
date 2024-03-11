% Combined adjoined operators% 
function BProj = RTT_adjoint(boundary_values,anz_r,alpha,approx, operator,epsi,step,type)
mode = strcat(approx, '_', operator);

switch mode
     %% Analytical & integration
    case 'ana_int'
    [anz_phi, anz_theta] = size(boundary_values);

    wxi1 = zeros(anz_r,anz_phi,anz_theta); % Integrand
    wxi2 = zeros(anz_r,anz_phi,anz_theta);
    
    h_r     = 1   /anz_r;
    h_phi   = 2*pi/anz_phi;
    h_theta = 2*pi/anz_theta;
    
    rho =  (h_r    :h_r    :1   );
    phi =  (h_phi  :h_phi  :2*pi);
    theta= (h_theta:h_theta:2*pi);
    
    parfor l=1:anz_r
        for p=1:anz_phi
            x  = rho(l)*[cos(phi(p))  ;sin(phi(p))  ];
            for q=1:anz_theta   
    
                xi =    [cos(theta(q));sin(theta(q))];
                
                tauP = tauplus(x,xi);
                
                endpoint = real((x+tauP*xi)/norm(x+tauP*xi));
                enddirection = xi;
    
                [ang, ~] = cart2pol(endpoint(1),endpoint(2));
                if ang<0
                    ang = ang +2*pi;
                end
                if ang>=2*pi
                   ang = ang -2*pi;
                end
                ang_klein = floor(ang/h_phi); % Smallest index of adjacent grid point
    
                if ang_klein == 0 || ang_klein == anz_phi % Euclidean: xi const, no interpolation over q 
                    weight        = ang/h_phi;
                    h_x_plus_tauP = (1-weight)*boundary_values(anz_phi,q) + weight*boundary_values(1,q);
                else
                    weight        = (ang-phi(ang_klein))/h_phi;
                   
                    h_x_plus_tauP = (1-weight)*(boundary_values(ang_klein,q))+weight*boundary_values(ang_klein+1,q);
                end
                if dot(endpoint,enddirection)<1e-4%^2+1-dot(endpoint,endpoint)<1e-4% Radikand klein -> x*xi ~0
                    w=0;
                else
                    
                   w = h_x_plus_tauP *exp(-alpha*tauP)/sqrt(dot(x,xi)^2+1-dot(x,x)); %Exponential ist nicht richtig für beliebige Metriken
                end
                
                wxi1(l,p,q) = w*xi(1);
                wxi2(l,p,q) = w*xi(2);
            end
        end
    end
    
    BProj = zeros(2,anz_r,anz_phi);
    
    parfor l=1:anz_r
        for p = 1:anz_phi
            BProj(:,l,p) = 2*pi/anz_theta * [sum(wxi1(l,p,:)); sum(wxi2(l,p,:))];
        end
    end

    %% Numerical & integration -> geodesic    
    case 'num_int'
        [anz_phi, anz_theta] = size(boundary_values);

        wxi1 = zeros(anz_r,anz_phi,anz_theta); % später Integrand
        wxi2 = zeros(anz_r,anz_phi,anz_theta);
        
        h_r     = 1   /anz_r;
        h_phi   = 2*pi/anz_phi;
        h_theta = 2*pi/anz_theta;
        
        rho =  (h_r    :h_r:1);
        phi =  (h_phi  :h_phi:2*pi);
        theta= (h_theta:h_theta:2*pi);
        
        parfor l=1:anz_r
            for p=1:anz_phi
                x  = rho(l)*[cos(phi(p));sin(phi(p))];
                for q=1:anz_theta
                    
                    xi = [cos(theta(q));sin(theta(q))];
                    
                    geodesic = Comp_geodesic_forward(x,xi,step,type); % Computes geodesic including gamma(tau_+)
                    tauP = (size(geodesic,2)-1)*step+geodesic(1,end);
                    geodesic(:,end)=[];
                    
                    endpoint     = geodesic(1:2,end); % Last point of geodesic
                    enddirection = geodesic(3:4,end); % Last direction of geodesic
                    
                    %% Error avoiding for x
                    if imag(endpoint(1))>0.001 || imag(endpoint(2))>0.001
                        fprintf('Something wrong with endpoint in Iaf_adjoint\n');
                    else 
                        endpoint(1)=real(endpoint(1)); 
                        endpoint(2)=real(endpoint(2));
                    end
                    [ang_phi, ~] = cart2pol(endpoint(1),endpoint(2));
                    if ang_phi<=0
                       ang_phi = ang_phi +2*pi;
                    end
                    if ang_phi>=2*pi
                       ang_phi = ang_phi-2*pi;
                    end
                    %% Error avoiding for xi
                    if imag(enddirection(1))>0.001 || imag(enddirection(2))>0.001
                        fprintf('Something wrong with endpoint in Iaf_adjoint\n');
                    else 
                        enddirection(1)=real(enddirection(1)); 
                        enddirection(2)=real(enddirection(2));
                    end
                    [ang_theta, ~] = cart2pol(enddirection(1),enddirection(2));
                    if ang_theta<=0
                        ang_theta = ang_theta +2*pi;
                    end
                    if ang_theta>=2*pi
                       ang_theta = ang_theta-2*pi;
                    end
                    %% Interpolation
                    ang_phi_small   = floor(ang_phi/h_phi); % Smalles index of adjecent grid point
                    ang_theta_small = floor(ang_theta/h_phi);

                    if (ang_phi_small == 0 || ang_phi_small == anz_phi) && (ang_theta_small == 0 || ang_theta_small == anz_theta)
                        
                        weight_phi   = ang_phi/h_phi;
                        weight_theta = ang_theta/h_theta;

                        weight00 = (1-weight_phi) * (1-weight_theta) * boundary_values(anz_phi,anz_theta);
                        weight01 = (1-weight_phi) * weight_theta     * boundary_values(anz_phi,1        ); 
                        weight10 = weight_phi     * (1-weight_theta) * boundary_values(1      ,anz_theta); 
                        weight11 = weight_phi     * weight_theta     * boundary_values(1      ,1        ); 
                        h_x_plus_tauP = weight00+weight01+weight10+weight11;

                    elseif (ang_phi_small == 0 || ang_phi_small == anz_phi) % theta is checked

                        weight_phi   = ang_phi/h_phi;
                        weight_theta = (ang_theta-theta(ang_theta_small))/h_theta;

                        weight00 = (1-weight_phi) * (1-weight_theta) * boundary_values(anz_phi,ang_theta_small  );
                        weight01 = (1-weight_phi) * weight_theta     * boundary_values(anz_phi,ang_theta_small+1); 
                        weight10 = weight_phi     * (1-weight_theta) * boundary_values(1      ,ang_theta_small  ); 
                        weight11 = weight_phi     * weight_theta     * boundary_values(1      ,ang_theta_small+1); 
                        h_x_plus_tauP = weight00+weight01+weight10+weight11;

                    elseif (ang_theta_small == 0 || ang_theta_small == anz_phi) % phi is checked
                        weight_phi   = (ang_phi-phi(ang_phi_small))/h_phi;
                        weight_theta = ang_theta/h_theta;

                        weight00 = (1-weight_phi) * (1-weight_theta) * boundary_values(ang_phi_small  ,anz_theta);
                        weight01 = (1-weight_phi) * weight_theta     * boundary_values(ang_phi_small  ,1        ); 
                        weight10 = weight_phi     * (1-weight_theta) * boundary_values(ang_phi_small+1,anz_theta); 
                        weight11 = weight_phi     * weight_theta     * boundary_values(ang_phi_small+1,1        ); 
                        h_x_plus_tauP = weight00+weight01+weight10+weight11;

                    else    %phi and theta checked                  
                        weight_phi   = (ang_phi   -     phi(ang_phi_small))/h_phi;
                        weight_theta = (ang_theta - theta(ang_theta_small))/h_theta;

                        weight00 = (1-weight_phi) * (1-weight_theta) * boundary_values(ang_phi_small  ,ang_theta_small  );
                        weight01 = (1-weight_phi) * weight_theta     * boundary_values(ang_phi_small  ,ang_theta_small+1); 
                        weight10 = weight_phi     * (1-weight_theta) * boundary_values(ang_phi_small+1,ang_theta_small  ); 
                        weight11 = weight_phi     * weight_theta     * boundary_values(ang_phi_small+1,ang_theta_small+1); 
                        h_x_plus_tauP = weight00+weight01+weight10+weight11;

                    end
                    if dot(endpoint,enddirection)<1e-5      
                        w=0;
                    else      
                        %disp(norm(geodesic(1:2,1)-geodesic(1:2,2)));
                        n_endpoint = refractive(endpoint); %n(x+), nabla(n(x+))
                        
                        step_star = norm(geodesic(1:2,end)-geodesic(1:2,end-1)); % last step size

                        XI = zeros(size(geodesic,2)); % Initialize Xi_n

                        for j = 1:size(geodesic,2)
                            ngradn = refractive(geodesic(1:2,j));
                            n      = ngradn(1);
                            gradn  = ngradn(2:3);
                            XI(j)  = 0.5*power(n,-1)*dot(gradn,geodesic(3:4,j));
                        end                 
                        expXI = exp(-0.5*step*(sum(XI(1:end-2)) + sum(XI(2:end-1))) - 0.5*step_star*(XI(end-1)+XI(end)));
                        w     = h_x_plus_tauP *exp(-alpha*tauP)*expXI/(power(n_endpoint(1),-2)*dot(endpoint,enddirection));
                    end
                    wxi1(l,p,q) = w*xi(1);
                    wxi2(l,p,q) = w*xi(2);                    
                end
            end
        end
        
        BProj = zeros(2,anz_r,anz_phi);
        
        for l=1:anz_r
            for p = 1:anz_phi
                BProj(:,l,p) = 2*pi/anz_theta * [sum(wxi1(l,p,:)); sum(wxi2(l,p,:))];
            end
        end

    %% Analytical & integration
    case 'ana_pde'
        [anz_phi, anz_theta] = size(boundary_values);

        ww = zeros(anz_phi,anz_theta);
        
        h_r     = 1   /anz_r;
        h_phi   = 2*pi/anz_phi;
        h_theta = 2*pi/anz_theta;
        
        rho =  (h_r    :h_r    :1   );
        phi =  (h_phi  :h_phi  :2*pi);
        theta= (h_theta:h_theta:2*pi);
        
        %% Computation of boundary values
            for p=1:anz_phi
                x  = [cos(phi(p))  ;sin(phi(p))  ];
                for q=1:anz_theta   
        
                    xi =    [cos(theta(q));sin(theta(q))];
                    
                    tauP = tauplus(x,xi); 
                    
                    endpoint = real((x+tauP*xi)/norm(x+tauP*xi));
                    enddirection = xi;
        
                    [ang, ~] = cart2pol(endpoint(1),endpoint(2));
                    if ang<0
                        ang = ang +2*pi;
                    end
                    if ang>=2*pi
                       ang = ang -2*pi;
                    end
                    ang_klein = floor(ang/h_phi); % Smallest index of adjacent grid point
        
                    if ang_klein == 0 || ang_klein == anz_phi % Euclidean: xi const, no interpolation over q 
                        weight        = ang/h_phi;
                        h_x_plus_tauP = (1-weight)*boundary_values(anz_phi,q) + weight*boundary_values(1,q);
                    else
                        weight        = (ang-phi(ang_klein))/h_phi;
                       
                        h_x_plus_tauP = (1-weight)*(boundary_values(ang_klein,q))+weight*boundary_values(ang_klein+1,q);
                    end
                    if dot(endpoint,enddirection)<1e-4  % Radicand small -> x*xi ~0
                        w=0;
                    else
                        
                       w = h_x_plus_tauP *exp(-alpha*tauP)/sqrt(dot(x,xi)^2+1-dot(x,x)); %Exponential ist nicht richtig für beliebige Metriken
                    end
                    
                    ww(p,q)=w;
                end
            end
        
        
        
        % Beschreiben und der Matrix und Lösung des LGS
          w_collect1 = zeros(anz_r,anz_phi,anz_theta);
          w_collect2 = zeros(anz_r,anz_phi,anz_theta);
        
        parfor q=1:anz_theta
            A = zeros(anz_r*anz_phi);
            b = zeros(anz_r*anz_phi,1);
          
            for r = 1:anz_r
                for p=1:anz_phi
                    if r==anz_r
                        A((r-1)*anz_phi+p, (r-1)*anz_phi+p) = 1;            % Allocation of boundary values
                        b((r-1)*anz_phi+p)                  = ww(p,q);
                    elseif r == 1 
        
                        cosine = cos(theta(q)-phi(p))/h_r;
                        sine   = sin(theta(q)-phi(p))/(2*rho(r)*h_phi);        
                       
                        A((r-1)*anz_phi+p,  r*anz_phi+p)        =  - cosine - epsi/(rho(r)*h_r) - 2*epsi/(3*h_r^2);                             
                        A((r-1)*anz_phi+p,  (r-1)*anz_phi+p)    =  + cosine + alpha +epsi/(rho(r)*h_r) + 1*epsi/(h_r^2) + 2*epsi/(rho(r)^2*h_r^2);       %(r,p)

                        index = mod(p+int8(anz_phi/2),anz_phi);
                        index(index ==0)=anz_phi;
                        A((r-1)*anz_phi+p,  index)              =  -epsi/(3*h_r^2); 
        
                        if p==anz_phi
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +1)       =   -sine - epsi/(rho(r)^2*h_phi^2);            %(r,p+1)
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +p-1)     =   +sine + epsi/(rho(r)^2*h_phi^2);            %(r,p-1)
                        elseif p==1
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +p+1)     =   -sine - epsi/(rho(r)^2*h_phi^2);            %(r,p+1)
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +anz_phi) =   +sine + epsi/(rho(r)^2*h_phi^2);            %(r,p-1)
                        else
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +p+1)     =   -sine - epsi/(rho(r)^2*h_phi^2);            %(r,p+1)
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +p-1)     =   +sine + epsi/(rho(r)^2*h_phi^2);            %(r,p-1)
                        end
        
                    else % r ~= 1
                        cosine = cos(theta(q)-phi(p))/h_r;
                        sine   = sin(theta(q)-phi(p))/(2*rho(r)*h_phi);
        
                        if r == anz_r -1 
                            b((r-1)*anz_phi+p)                      =  +(cosine + epsi/(rho(r)*h_r) + epsi/(h_r^2)) * ww(p,q);
                        else
                            A((r-1)*anz_phi+p,  r*anz_phi+p)        =  - cosine - epsi/(rho(r)*h_r) - epsi/(h_r^2);        
                        end
                        A((r-1)*anz_phi+p,  (r-1)*anz_phi+p)        =  +cosine + alpha +epsi/(rho(r)*h_r) + 2*epsi/(h_r^2) + 2*epsi/(rho(r)^2*h_r^2);       %(r,p)
                        A((r-1)*anz_phi+p,  (r-2)*anz_phi+p)        =  -epsi/(h_r^2);
                        
                        if p==anz_phi
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +1)       =   -sine - epsi/(rho(r)^2*h_phi^2);            %(r,p+1)
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +p-1)     =   +sine + epsi/(rho(r)^2*h_phi^2);            %(r,p-1)
                        elseif p==1
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +p+1)     =   -sine - epsi/(rho(r)^2*h_phi^2);            %(r,p+1)
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +anz_phi) =   +sine + epsi/(rho(r)^2*h_phi^2);            %(r,p-1)
                        else
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +p+1)     =   -sine - epsi/(rho(r)^2*h_phi^2);            %(r,p+1)
                            A((r-1)*anz_phi+p,  (r-1)*anz_phi +p-1)     =   +sine + epsi/(rho(r)^2*h_phi^2);            %(r,p-1)
                        end
        
                    end
        
                end
            end
            f = lsqminnorm(A,b);
            for r=1:anz_r
                for p=1:anz_phi
                    w_collect  = f((r-1)*anz_phi+p); % Computes w 
                    w_collect1(r,p,q) = cos(theta(q)) * w_collect;
                    w_collect2(r,p,q) = sin(theta(q)) * w_collect;
                end
            end
        end
        
        % Extraction of iterated
        adj1 = 2*pi/anz_theta  *  sum(w_collect1,3);
        adj2 = 2*pi/anz_theta  *  sum(w_collect2,3);
        
        BProj = zeros(2,anz_r,anz_phi);
        BProj(1,:,:) = adj1;
        BProj(2,:,:) = adj2;

    %% Numerical & PDE -> geodesic    
    case 'num_pde' %geodesic
         disp('sinnlos, da riesige Laufzeiten')
    otherwise
        disp('Error in adjoint!');
        
end