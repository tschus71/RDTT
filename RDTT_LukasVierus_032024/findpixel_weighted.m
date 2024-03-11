function value = findpixel_weighted(x,field,h_rad,h_phi)
if nargin == 4
    x=real(x);
    [ang, rad] = cart2pol(x(1),x(2));

    if ang<0
        ang = ang +2*pi;
    end
    while rad >=0.9999 % points on the boundary are shifted to the outer ring
        rad = 0.9999*rad;
    end
    if ang>= 2*pi
        ang = ang-2*pi;
    end
    
    ang_klein = floor(ang/h_phi); % Smallest angular index of adjacent grid point
    r_klein   = floor(rad/h_rad); % Smallest radial index of adjacent grid point
    
    
    r = (h_rad:h_rad:1); 
    phi = (h_phi: h_phi:2*pi); 
    phi_end = size(phi,2);
    
    if ang_klein > 0 
        if r_klein > 0
           
            k1 = (rad-r(r_klein))/h_rad;
            k2 = (ang - phi(ang_klein))/h_phi;
    
            Q11 = [r_klein,ang_klein];          % point(:,1);
    
            if ang_klein == size(field,2)       % Problem bei Übersprung 2pi -> 0
                Q12 = [r_klein,1];              % point(:,2); %D
                Q22 = [r_klein+1,1];            % point(:,3); %C
            else
                Q12 = [r_klein,ang_klein+1];    % point(:,2); %D
                Q22 = [r_klein+1,ang_klein+1];  % point(:,3); %C
            end
            Q21 = [r_klein+1,ang_klein];        % point(:,4); %B
           
            order0 = field(Q11(1),Q11(2));
            order1 = (field(Q21(1),Q21(2))-field(Q11(1),Q11(2)))*k1 + (field(Q12(1),Q12(2))-field(Q11(1),Q11(2)))*k2;
            order2 = (field(Q11(1),Q11(2))+field(Q22(1),Q22(2))-field(Q12(1),Q12(2))-field(Q21(1),Q21(2)))*k1*k2;
            value = order0 +order1+ order2;
    
        elseif r_klein == 0 % only 3 grid points
            value_in_0 = mean(field(1,:));
            k1 = rad/h_rad;
            k2 = (ang - phi(ang_klein))/h_phi;
        
            if ang_klein == size(field,2)
                value = k1*value_in_0  + (1-k1)*(k2*field(1,ang_klein) + (1-k2)*field(1,1));
            else
                value = k1*value_in_0  + (1-k1)*(k2*field(1,ang_klein) + (1-k2)*field(1,ang_klein+1));
            end
        else
            fprintf('Something wrong in findpixel\n');
        end
    elseif ang_klein == 0
      
        if r_klein > 0
    
            k1 = (rad-r(r_klein))/h_rad;
            k2 = ang/h_phi;
    
            Q11 = [r_klein,phi_end];    %point(:,1); %A
            Q12 = [r_klein,1 ];         %point(:,2); %D
            Q22 = [r_klein+1,1];        %point(:,3); %C
            Q21 = [r_klein+1,phi_end];  %point(:,4); %B
            
            order0 = field(Q11(1),Q11(2)); 
            order1 = (field(Q21(1),Q21(2))-field(Q11(1),Q11(2)))*k1 + (field(Q12(1),Q12(2))-field(Q11(1),Q11(2)))*k2;
            order2 = (field(Q11(1),Q11(2))+field(Q22(1),Q22(2))-field(Q12(1),Q12(2))-field(Q21(1),Q21(2)))*k1*k2; 
            value  = order0 +order1 +order2;
    
        elseif r_klein == 0
            value = field(1,ang_klein+1);
        else
            fprintf('Something wrong in findpixel\n');
        end    
    end
end