function speed = Soundspeed(x, type)

switch(type)
    case 0
        c=1;
        dc_x1=0;
        dc_x2=0;

    case 1
        x_0 = [0;0];
        nu = 0.02; % nu(1+|x_0|) !< alpha

        c     =  1/(nu*power(norm(x-x_0),2)+1);
        dc_x1 = -2*nu/(nu*power(norm(x-x_0),2)+1)*(x(1)-x_0(1));
        dc_x2 = -2*nu/(nu*power(norm(x-x_0),2)+1)*(x(2)-x_0(2));

    case 2
        x_0  = [0;0];
        lamb = 0.01;
        mu   = 0.01;

        c     = 1/(1+lamb*exp(-mu*power(norm(x-x_0),2)));
        dc_x1 = 2*lamb*mu*exp(-mu*power(norm(x-x_0),2))*(x(1)-x_0(1))/(1+lamb*exp(-mu*power(norm(x-x_0),2)))^2;
        dc_x2 = 2*lamb*mu*exp(-mu*power(norm(x-x_0),2))*(x(2)-x_0(2))/(1+lamb*exp(-mu*power(norm(x-x_0),2)))^2;

    case 3 % n = 2-x1*x2
        c     = 1/(2-x(1)*x(2));
        dc_x1 = x2/(2-x(1)*x(2))^2;
        dc_x2 = x1/(2-x(1)*x(2))^2;
 end

speed=[c, dc_x1, dc_x2];
end