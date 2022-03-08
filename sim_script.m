pars.epsilon=0.1;
pars.I__d=0.01;
pars.h=20;
pars.alpha=1;
pars.k=0;
pars.T=500;
pars.c=0.01;
pars.c_theta=0;
pars.omega=1.4;
pars.delta=0.0;
pars.g=0;
%vary one parameter
pars.delta=linspace(0,0.4,11)';

% pars.epsilon=0.1;pars.I__d=0.01;pars.h=20;pars.alpha=1;pars.k=0;
% pars.T=500;pars.c=0.01;pars.c_theta=0;pars.omega=1.4;pars.delta=0.2;pars.g=0;
% vary one parameter
% pars.omega=[1.58973459023457,1.389684894,1.3236765479,1.202356,1.1496545678]';

close all
% Run integrator
% inputs are the parameter values for the integration (only one can be a
% vector
% outputs are:
%   sols - a cell array of ode solution structures
%   Pars - a structure of parameter values for a parameter sweep
[sols,Pars]=sim_wrap(pars);
% Run analysis on the data
make_data_figures(sols,Pars)