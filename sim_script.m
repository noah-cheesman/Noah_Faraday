pars.epsilon=0.1;
pars.beta=0.01;
pars.h=0;
pars.isotropy=1;
pars.T=5;
pars.k=0.1/1.4;
pars.k=[1/1.7,1/1.4,1,1.4,1.7,1.4^2,1.7^2]';
%pars.alpha=linspace(1/1.7,1/1.65,10)';
%pars.alpha=linspace(0.59814,0.61028,5)';
%pars.alpha=0.9/1.4;
pars.alpha=0;
pars.c=0.05;
pars.c_theta=0;
pars.isotropy=1;
pars.omega=1;
pars.delta=0.16235;
pars.delta=0;
%pars.delta=0.12;
%pars.delta=fliplr(linspace(0.08,0.09,6))';
pars.g=0;
pars.epsilon=0.1;

%vary one parameter
pars.epsilon=0.1;

disp(['delta_crit=',num2str((pars.epsilon.*sqrt(1-pars.c.^2)./(pars.k+pars.alpha))')]);
Wic=[];

NUM=max([length(pars.c),length(pars.k),length(pars.epsilon),length(pars.omega),length(pars.delta),length(pars.alpha)]);
for num=1:NUM
    epsilon=pars.epsilon(min([num,length(pars.epsilon)]));
    omega=pars.omega(min([num,length(pars.omega)]));
    c=pars.c(min([num,length(pars.c)]));
    k=pars.k(min([num,length(pars.k)]))+pars.alpha(min([num,length(pars.alpha)]));
x0 = epsilon*(-omega^2+k)*omega.^2/(omega^4+(c^2-2*k)*omega^2+k^2);
y0 = -omega^3/(omega^4+(c^2-2*k)*omega^2+k^2)*c*epsilon;
u0 = epsilon*omega^4*c/(omega^4+(c^2-2*k)*omega^2+k^2);
v0 = epsilon*(-omega^2+k)*omega^3/(omega^4+(c^2-2*k)*omega^2+k^2);
theta0=0;
thetadot0=omega;

%Wic=[Wic;x0,y0,theta0,u0,v0,thetadot0; x0+0.01,y0,theta0,u0,v0+0.01,thetadot0; x0+0.2,y0,theta0,u0,v0+0.02,thetadot0;x0,y0,theta0,u0,v0,thetadot0+0.02; x0+0.05,y0,theta0,u0,v0+0.05,thetadot0;];
% Wic=[Wic;...
%     x0,y0,theta0,u0,v0*3.20,thetadot0;...
%     x0,y0,theta0,u0,v0*1.00,thetadot0;...
%     x0,y0,theta0,u0,v0*1.20,thetadot0;...
%     x0,y0,theta0,u0,v0*1.40,thetadot0;...
%     x0,y0,theta0,u0,v0*1.60,thetadot0;...
%     x0,y0,theta0,u0,v0*1.80,thetadot0;...
%     x0,y0,theta0,u0,v0*2.00,thetadot0...
%     ];
% Wic=[Wic;x0,y0,theta0,10*u0,10*v0,thetadot0;...
%     -0.1658,-0.2747,0,0.0961,-0.2159,1.0000;...
%     2*x0,2*y0,2*theta0,2*u0,v0,2*thetadot0;...
%     ];
%Wic=[x0,y0,theta0,10*u0,10*v0,thetadot0];
Wic=[      0.1190,   -0.1701,         0,    0.0633,    0.1814,    1.0000];
%Wic=[      0.3,   0,         0,    0,    0.22,    1.0000;...
%    0.1,0,0,0,0.1,1;...
%    -0.1,0,0,0,0.05,1;...
%    -0.176,0,0,0,0,1;...
%    -0.176,0,0,0,-0.1,1];

end
close all

% Run integrator
% inputs are the parameter values for the integration (only one can be a
% vector
% outputs are:
%   sols - a cell array of ode solution structures
%   Pars - a structure of parameter values for a parameter sweep
[sols,Pars]=sim_wrap(pars,Wic);

% Run analysis on the data and save to folder name
% be careful not to overwrite
foldername='linear_is_boring3';
mkdir(foldername)
save_parameters(foldername,pars)
writematrix(Wic,[foldername,'/ICS.txt']);
%make_data_figures(sols,Pars,foldername);
if 1==1
save_data_figures(sols,Pars,foldername)
end
disp('finished')