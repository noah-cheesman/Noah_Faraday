function [sol]=integrator(epsilon,I__d,h,alpha,k,T,c,c_theta,omega,delta,g,isotropy,Wic)

% if no input parameters use these defaults
if nargin == 0
    epsilon=0.1;I__d=0.001;h=10;alpha=1;k=0;T=5000;c=0.01;c_theta=0.0000;omega=1.4;delta=0.4;g=0;
    Wic=[0.5;0.01;0.01;0.1;0;0];
end

% set ODE options structure: use mass matrix and events functions found
% below
opts1 = odeset('Mass',@(t,y)mass(t,y,epsilon,I__d),'Events',@(t,y)event_section(t,y));
opts2 = odeset('Mass',@(t,y)mass(t,y,epsilon,I__d),'Events',@(t,y)event_section(t,y),'AbsTol',1e-8,'RelTol',1e-8);

% solve using stiff solver time used is 5000 (I should either 
% programatically decide a good time span or use an event function to give 
% some idea of convergence to a periodic orbit
sol_t=ode15s(@(t,X)rhs(t,X,epsilon,h,alpha,k,T,c,c_theta,omega,delta,g,isotropy),[0,5000],Wic,opts1);
IC=sol_t.y(:,end);
sol=ode15s(@(t,X)rhs(t,X,epsilon,h,alpha,k,T,c,c_theta,omega,delta,g,isotropy),[0,2000],IC,opts2);


% mass matrix 
    function Mass=mass(~,y,epsilon,I__d)
        Mass=[...
            1,0,0,0,0,0;...
            0,1,0,0,0,0;...
            0,0,1,0,0,0;...
            0,0,0,1,0,-epsilon*sin(y(3));...
            0,0,0,0,1,+epsilon*cos(y(3));...
            0,0,0,y(2)-epsilon*sin(y(3)),-y(1)+epsilon*cos(y(3)),I__d+epsilon^2-epsilon*y(1)*cos(y(3))-epsilon*y(2)*sin(y(3))];
    end

    function [dX]=rhs(~,X,epsilon,h,alpha,k,T,c,c_theta,omega,delta,g,isot)
        heav=(1+sign(sqrt(X(1)^2+X(2)^2)-delta))/2;
        %f=alpha*heav*(1-delta/sqrt(X(1)^2+X(2)^2))*h^2/(h^2-X(1)^2-X(2)^2);
        f=alpha*heav*(1-delta/sqrt(X(1)^2+X(2)^2));
        dX=[0,0,0,1,0,0;...
            0,0,0,0,1,0;...
            0,0,0,0,0,1;...
            -k*isot-f*isot,0,0,-c,0,0;...
            0,-k-f,0,0,-c,0;...
            0,0,0,0,0,-c_theta-T]*X+epsilon*X(6)^2*[0;0;0;cos(X(3));sin(X(3));-X(1)*sin(X(3))+X(2)*cos(X(3))]+[0;0;0;-g;0;T*omega];
    end



    function [value,isterminal,direction] = event_section(~,y)
        value = sin(y(3)); % Poincaré section at sin(y(3))
        isterminal = 0;    % Stop the integration
        direction = 1;     % Only when increasing (i.e. once per period)
    end

end