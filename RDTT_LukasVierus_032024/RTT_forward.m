% Combine all forward operators

function boundary_values = RTT_forward(set_array1, set_array2,anz_r,anz_phi,anz_theta,alpha,approx,step,type)

    switch approx
        %% Analytical & integration
        case 'ana' 
            h_r     = 1/anz_r;
            h_phi   = 2*pi/anz_phi;
            h_theta = 2*pi/anz_theta;
            
            phi   = (h_phi:h_phi:2*pi);
            theta = (h_theta:h_theta:2*pi);
            
            boundary_values = zeros(anz_phi, anz_theta);
            
            for iter_xi = 1:anz_theta
                xi = [cos(theta(iter_xi));sin(theta(iter_xi))];
            
                    for j = 1:anz_phi
                         x = [cos(phi(j));sin(phi(j))];
                          
                         tauM      = tauminus(x,xi);
                         delta_tau = linspace(tauM,0);
                         geodesic  = zeros(4,numel(delta_tau));
                         geodesic(1:2,:) = x+delta_tau.*xi;               
                         
                         if size(geodesic,2)~=0
                            argument_f  =  geodesic(1:2,:);
                            no_points   = size(geodesic,2);%
                            step = abs(tauM)/(no_points-1);              
                            integrand1 = zeros(1,no_points);
                            integrand2 = zeros(1,no_points);
                            for k = 1:no_points
                                vecc = argument_f(:,k);
                                                 
                                expo = exp(alpha*delta_tau(k));
                                 
                                integrand1(k) = findpixel_weighted(vecc,set_array1,h_r,h_phi)*xi(1)*expo;
                                integrand2(k) = findpixel_weighted(vecc,set_array2,h_r,h_phi)*xi(2)*expo;     
                            end
                            
                            Q1 = step * sum((integrand1(1:end-1) + integrand1(2:end)))/2 ; % Value of integral in 1st component
                            Q2 = step * sum((integrand2(1:end-1) + integrand2(2:end)))/2 ; % Value of integral in 2nd component
                            Q = Q1+ Q2;
                            boundary_values(j,iter_xi)= Q;
                         end
                    end
            end
            boundary_values(abs(boundary_values)<0.0001)=0;%0.001 

        %% Numerical  & integration
        case 'num' 
            h_r     = 1/anz_r;
            h_phi   = 2*pi/anz_phi;
            h_theta = 2*pi/anz_theta;
            
            phi   = (h_phi:h_phi:2*pi);
            theta = (h_theta:h_theta:2*pi);
            
            boundary_values = zeros(anz_phi, anz_theta);
            
            parfor iter_xi = 1:anz_theta
                xi = [cos(theta(iter_xi));sin(theta(iter_xi))];
            
                    for j = 1:anz_phi
                         x = [cos(phi(j));sin(phi(j))];
                          
                         geodesic = Comp_geodesic_single(x,xi,step,type); % gamma(tau_minus,...,0)
                         geodesic = flip(geodesic')';  % flips entries of matrix columnwise
                         if size(geodesic,2)~=0
                            last_step = geodesic(1,1); %Extracts step_last, it is now at first place
                            geodesic(:,1)=[]; % Delete last column
                            argument_f  =  geodesic(1:2,:);
                            gamma_dot   = -geodesic(3:4,:); % Minus inserted because integration goes in other direction
                            no_points   = size(geodesic,2);
                            alpha_value = alpha*ones(1,no_points); 

                            integrand1 = zeros(1,no_points);
                            integrand2 = zeros(1,no_points);
                            for k = 1:no_points
                                vecc = argument_f(:,k);
                                used_alphas = alpha_value(k:no_points); % for expo use all tau from tau_i to outflow (1)
                                adjust_expo=0;
                                if k == 1
                                   adjust_expo = (step-last_step)*(used_alphas(1)+used_alphas(2))/2;
                                end
                
                                expo = exp(-step * sum(used_alphas(1:end-1) + used_alphas(2:end))/2 + adjust_expo);
                                 
                                integrand1(k) = findpixel_weighted(vecc,set_array1,h_r,h_phi)*gamma_dot(1,k)*expo;
                                integrand2(k) = findpixel_weighted(vecc,set_array2,h_r,h_phi)*gamma_dot(2,k)*expo;
                            
                            end
                            adjust1 = (step-last_step)*(integrand1(1)+integrand1(2))/2; % Compensation for last step size
                            adjust2 = (step-last_step)*(integrand2(1)+integrand2(2))/2;
                            Q1 = step * sum((integrand1(1:end-1) + integrand1(2:end)))/2 - adjust1; % Value of integral in 1st component
                            Q2 = step * sum((integrand2(1:end-1) + integrand2(2:end)))/2 - adjust2; % Value of integral in 2nd component
                            Q = Q1+ Q2;
                            boundary_values(j,iter_xi)= Q;
                         end
                    end
            end
            boundary_values(abs(boundary_values)<0.0001)=0;
        otherwise
            disp('Error!');
    end
end


