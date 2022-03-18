function [sols,Pars]=sim_wrap(pars,wics)
%close all
diffwix=0;
if nargin==0
%constant values of each parameter
pars.epsilon=0.1;
pars.beta=0.001;
pars.h=2;
pars.alpha=1;
pars.k=0;
pars.T=100;
pars.c=0.1;
pars.c_theta=0;
pars.omega=1.4;
pars.delta=0.2;
pars.g=0;
%vary one parameter (MUST BE A COLUMN VECTOR)
pars.omega=[1.6,1.3,1.2,1.15]';
end

%my grid makes all
Pars=mygrid(pars);
dims=size(Pars.delta);N2=max(dims);N3=size(wics,1);
N=max([N2,N3]);

% Make empty cell array for solution structures
sols={};
if max(dims(dims<N))>1
    disp('only vary one parameter');
else
    for itemp=1:N
        i=min([itemp,N2]);
        % if you'd like to hard code the IC
        if itemp==1 || diffwix==0
        %Wic=[0.01;0;0;0;0.01;Pars.omega(i)];
        Wic=wics(min([itemp,N3]),:)';
        end

        sol=integrator(Pars.epsilon(i),Pars.beta(i),Pars.h(i),...
            Pars.alpha(i),Pars.k(i),Pars.T(i),Pars.c(i),Pars.c_theta(i),...
            Pars.omega(i),Pars.delta(i),Pars.g(i),Wic);
        W=sol.y;
        sols{itemp}=sol;
        Wic=W(:,end);
        
        
        % to check on integration progress
        disp([num2str(100*itemp/N),'% complete'])
        
    end
end


disp('done');

end

function parsout=mygrid(parsin)
l=0;
list_of_fields=fieldnames(parsin);
for i=1:length(list_of_fields)
    l=max([l,length(parsin.(list_of_fields{i}))]);
end
for i=1:length(list_of_fields)
    parsout.(list_of_fields{i})=parsin.(list_of_fields{i}).*ones(l,1);
end
end
