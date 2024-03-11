% Combined Landweber
alpha = 0;    % attenuation
delta = 0;    % Noise level
epsi  = 0;    % visco parameter
relax = 0.1;  % relaxation in Landweber
step = 0.005; % step size Runke-Kutta
beta = 0;     % Tikhonov shift
type = 0;     % type of metric
%% 

approx   = 'ana'; % 'ana' for analyitcal or 'num' for numerical
operator = 'pde'; % 'int' or 'pde' for adjoint

field1 =@(x) [1;0]; field2 =@(x) [x(1);-x(2)]; field3 =@(x) [x(1)+x(2);x(1)-x(2)]; field5 =@(x) [x(1)^2-2*x(2)^2; -2*x(1)*x(2)]; 
field = field5;

[anz_r,anz_phi,anz_theta] = deal(34,106,106); % Set grid 
[fieldPOL1, fieldPOL2] = vec2polargrid(field, anz_r, anz_phi); % Compute field on grid

normF = product_preimage(fieldPOL1,fieldPOL2,fieldPOL1,fieldPOL2); % ||f||

y_delta = RTT_forward(fieldPOL1, fieldPOL2, anz_r,anz_phi, anz_theta, alpha,approx,step,type); % generating data
y_delta = y_delta .*( 2*delta * rand(size(y_delta))+(1-delta).*ones(size(y_delta))); % noising

F_iter1 = zeros(size(fieldPOL1)); F_iter2 = zeros(size(fieldPOL2)); % initial guess f=0
    
e_new =2; e_old =3;iter = 0;
    
F1_new = F_iter1; % Set dummy variables for Netserov accelaration
F1_old = F_iter1;
F2_new = F_iter2;
F2_old = F_iter2;

while e_new <= e_old  && abs(e_new-e_old)/normF>0.00001
      iter = iter +1;

      lambda = 0;%(iter-1)/(iter+2);                       % Nesterov step
      Z_iter1 = F1_new + lambda * (F1_new - F1_old);
      Z_iter2 = F2_new + lambda * (F2_new - F2_old);

      Sf = RTT_forward(F_iter1, F_iter2, anz_r,anz_phi, anz_theta, alpha,approx,step,type);
    
      e_old = product_preimage(F_iter1-fieldPOL1,F_iter2-fieldPOL2,F_iter1-fieldPOL1,F_iter2-fieldPOL2);

       fprintf("Relative error ||f-f_exact||/||f_exact|| = %f\n", e_old/normF);
        
      Residual = Sf - y_delta;
        
      norm_residual = product_image(Residual); % ||Sf_n - y^delta||
        
      R = RTT_adjoint(Residual,anz_r,alpha,approx, operator,epsi,step,type);
        
      R1 = reshape(R(1,:,:), anz_r,anz_phi);
      R2 = reshape(R(2,:,:), anz_r,anz_phi);

      F_iter1 = (1-relax*beta)*Z_iter1 - relax*R1; % Iteration step for Landweber/Tikhonov 
      F_iter2 = (1-relax*beta)*Z_iter2 - relax*R2;
      
      F1_old = F1_new ; F2_old = F2_new;  % Set old iterated
      F1_new = F_iter1; F2_new = F_iter2; % Set new iterated

      e_new = product_preimage(F_iter1-fieldPOL1,F_iter2-fieldPOL2,F_iter1-fieldPOL1,F_iter2-fieldPOL2);    
end   
fprintf("Relative error ||f-f_exakt||/||f_exakt|| = %f\n", e_old/normF);