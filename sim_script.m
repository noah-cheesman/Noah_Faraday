pars.epsilon=0.1;
pars.beta=0.01;
pars.h=20;
pars.alpha=0.9;
pars.k=0.1;
pars.T=500;
pars.c=0.01;
pars.c_theta=0;
pars.omega=1.4;
pars.delta=0.12;
pars.g=0;
%vary one parameter
%pars.delta=linspace(0,0.4,11)';
pars.delta=linspace(0,0.1,5)';
%pars.k=linspace(0.5,1.3,5)';
Wic=[   0.01,0,0,0,0.01,pars.omega;...
        0.1,0,0,0,0.1,pars.omega;...
        0.5,0,0,0,0.5,pars.omega;];
Wic=[   0.01,0,0,0,0.01,pars.omega;...
        0.1,0,0,0,0.1,pars.omega];
    Wic=[0.09,0,0,0,0.095,pars.omega];
    
    
%pars.epsilon=0.1;pars.I__d=0.01;pars.h=20;pars.alpha=1;pars.k=0;
%pars.T=500;pars.c=0.01;pars.c_theta=0;pars.omega=1.4;pars.delta=0.2;pars.g=0;
% vary one parameter
%pars.omega=[1.58973459023457,1.389684894,1.3236765479,1.202356,1.1496545678]';

close all
% Run integrator
% inputs are the parameter values for the integration (only one can be a
% vector
% outputs are:
%   sols - a cell array of ode solution structures
%   Pars - a structure of parameter values for a parameter sweep
[sols,Pars]=sim_wrap(pars,Wic);
% Run analysis on the data
currDate = strrep(datestr(datetime), ':', '_');
mkdir(currDate)
save_parameters(currDate,pars)

for i =1:length(sols)
    sol=sols{i};
    [omega,omega_s]=omegas(pars,i);
    [ires,freq_vec,Spec,t]=perfect(sol,omega);
    freq_vec=freq_vec/omega;
    freq=omega_s;
    freq2=freq_vec;
    Spec=Spec(freq2<5);
    freq2=freq2(freq2<5);
    writetable(array2table([t(ceil(0.9*end:end))',ires(ceil(0.9*end:end))'],'VariableNames',{'t','ires'}),[currDate,'/IRES_',num2str(i),'.txt']);
    writetable(array2table([freq2',Spec'],'VariableNames',{'freq','Spec'}),[currDate,'/IRESfft_',num2str(i),'.txt']);
end

for i =1:length(sols)
    sol=sols{i};
    [omega,omega_s]=omegas(pars,i);
    [freq_vec,Spec,Z,cfreq,x,t]=perfect_x(sol,omega);
    freq_vec=freq_vec/omega;
    freq=omega_s;
    freq2=freq_vec;
    Spec=Spec(freq2<5);
    freq2=freq2(freq2<5);
    writetable(array2table([t(ceil(0.9*end:end))',x(ceil(0.9*end:end))'],'VariableNames',{'t','x'}),[currDate,'/X_',num2str(i),'.txt']);
    writetable(array2table([freq2',Spec'],'VariableNames',{'freq','Spec'}),[currDate,'/Xfft_',num2str(i),'.txt']);
end

%make_data_figures(sols,Pars,currDate); 
%save_data_figures(sols,Pars,currDate)
%