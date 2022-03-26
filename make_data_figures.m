function make_data_figures(sols,pars,currDate)
[parname,parvalues]=find_varying_parameter(pars,sols);
xtix=cellstr(num2str((parvalues)))';

colours=[[0, 0.4470, 0.7410];...
[0.8500, 0.3250, 0.0980];...
[0.9290, 0.6940, 0.1250];...
[0.4940, 0.1840, 0.5560];...
[0.4660, 0.6740, 0.1880];...
[0.3010, 0.7450, 0.9330];...
[0.6350, 0.0780, 0.1840]];



figures(1)=figure('Units','Centimeters','Position',[7,7,7,7]);hold on
for i =1:length(sols)
    W=sols{i}.y;
    % how much would you like to plot? Last 0.1 of time?
    last=0.3;
    W=W(:,ceil(end*(1-last)):end);
    % plot orbits in yz with parameter iterate in x axis
    plot3(i*ones(size(W(1,:))),W(1,:),W(2,:),'-')
end
xlabel(parname)
xticks(1:length(parvalues));
xticklabels(xtix)
view(-45,30)
ylabel('$x$','Interpreter','Latex')
zlabel('$y$','Interpreter','Latex')
title('Orbits in stationary frame','Interpreter','Latex')

figures(2)=figure('Units','Centimeters','Position',[7,7,7,7]);hold on
for i =1:length(sols)
    W=sols{i}.y;
    % how much would you like to plot? Last 0.1 of time?
    last=0.1;
    W=W(:,ceil(end*(1-last)):end);
    % plot orbits in yz with parameter iterate in x axis
    plot3(i*ones(size(W(1,:)))+0.25*(W(4,:))/max(W(4,:)),W(1,:),W(2,:),'-')
end
xlabel(parname)
xticks(1:length(parvalues));
xticklabels(xtix)
view(-45,30)
ylabel('$x$','Interpreter','Latex')
zlabel('$y$','Interpreter','Latex')
title('Orbits in stationary frame','Interpreter','Latex')

figures(3)=figure('Units','Centimeters','Position',[7,7,7,7]);hold on
for i =1:length(sols)
    W=sols{i}.y;
    last=0.1;
    W=W(:,ceil(end*(1-last)):end);
    % transform to rotating frame
    xs=W(1,:).*cos(W(3,:))+W(2,:).*sin(W(3,:));
    ys=-W(1,:).*sin(W(3,:))+W(2,:).*cos(W(3,:));
    % as with previous plot but now in rotating frame
    plot3(i*ones(size(xs)),xs,ys,'-','color',colours(i,:))
    plot3(i,xs(end),ys(end),'.','color',colours(i,:))
end
xlabel(parname)
xticks(1:length(parvalues));
xticklabels(xtix)
view(-45,30)
ylabel('$\tilde{x}$','Interpreter','Latex')
zlabel('$\tilde{y}$','Interpreter','Latex')
title('Orbits in rotating frame','Interpreter','Latex')

% plots hi low bifurcation diagram
figures(4)=figure('Units','Centimeters','Position',[7,7,7,7]);hold on
for i =1:length(sols)
    W=sols{i}.y;
    last=0.1;
    W=W(:,ceil(end*(1-last)):end);
    % transform to rotating frame
    xs=W(1,:).*cos(W(3,:))+W(2,:).*sin(W(3,:));
    %ys=-W(1,:).*sin(W(3,:))+W(2,:).*cos(W(3,:));
    % finds maxima and minima of lateral x in rotating frame
    Whi=+findpeaks(+xs);
    Wlo=-findpeaks(-xs);
    plot(i*ones(size(Whi)),Whi,'x')
    plot(i*ones(size(Wlo)),Wlo,'o')
end
xlabel(parname)
xticks(1:length(parvalues));
xticklabels(xtix)
ylabel('$\tilde{x}$','Interpreter','Latex')
title('min/max $\tilde{x}$ bif diagram','Interpreter','Latex')


figures(5)=figure('Units','Centimeters','Position',[7,7,7,7]);hold on
for i =1:length(sols)
    W=sols{i}.y;
    t=sols{i}.x;
    last=0.1;
    W=W(:,ceil(end*(1-last)):end);
    t=t(ceil(end*(1-last)):end);
    plot(t,((W(6,:))),'-');
    if length(pars.T)>1
    omega_s=pars.omega(i)*pars.T(i)/(pars.T(i)+pars.c_theta(i));
    else
        omega_s=pars.omega*pars.T/(pars.T+pars.c_theta);
    end
    plot([min(t),max(t)],omega_s*[1,1],'--')
end
ylabel('$\dot{\theta}$','Interpreter','Latex');
xlabel('$t$','Interpreter','Latex');
title('angular velocity of rotor compared to driven speed','Interpreter','Latex')


figures(6)=figure('Units','Centimeters','Position',[7,7,7,7]);hold on
for i =1:length(sols)
    We=sols{i}.ye;
    last=0.1;
    We=We(:,ceil(end*(1-last)):end);
    plot(i*ones(size(We(1,:))),We(1,:),'o')
end
ylabel('$x$','Interpreter','Latex');
xlabel(parname)
xticks(1:length(parvalues));
xticklabels(xtix)
title('Poincare section bif diagram','Interpreter','Latex')


figures(7)=figure('Units','Centimeters','Position',[7,7,7,7]);hold on
maxi=[];
mini=[];
xmax=1;
for i =1:length(sols)
    sol=sols{i};
    if length(pars.omega)>1
        omega=pars.omega(i);
    else
        omega=pars.omega;
    end
    [~,freq_vec,Spec]=perfect(sol,omega);
    freq=omega^2;
    plot3(i*ones(size(freq_vec)),freq_vec/freq,log(Spec),'-')
    maxi=max([maxi,max(log(Spec))]);
    mini=min([mini,min(log(Spec))]);
    ylim([ 0 xmax])
end
plot_farey_3d(mini,maxi,ceil(xmax))
ylabel('$\omega/\Omega$','Interpreter','Latex');
zlabel('log$|F(\omega/\Omega|$','Interpreter','Latex');
xlabel(parname)
xticks(1:length(parvalues));
xticklabels(xtix)
view(90,00);
title('fft of $I_\mathrm{res}$','Interpreter','Latex');

figures(8)=figure('Units','Centimeters','Position',[7,7,7,7]);hold on
maxi=[];
mini=[];
for i =1:length(sols)
    sol=sols{i};
    if length(pars.omega)>1
        omega=pars.omega(i);
    else
        omega=pars.omega;
    end
    [freq_vec,Spec,Z,cfreq]=perfect_x(sol,omega);
    freq=omega^2;
    plot3(i*ones(size(freq_vec)),freq_vec/freq,log(Spec),'-')
    maxi=max([maxi,max(log(Spec))]);
    mini=min([mini,min(log(Spec))]);
    ylim([ 0 xmax])
end
%plot_farey_3d(mini,maxi,ceil(xmax))
ylabel('$\omega/\Omega$','Interpreter','Latex');
title('fft of $x$','Interpreter','Latex');
view(90,00)

figures(9)=figure('Units','Centimeters','Position',[7,7,7,7]);hold on;
maxi=[];
mini=[];
for i =1:length(sols)
    sol=sols{i};
    [omega,omega_s]=omegas(pars,i);
    [freq_vec,Spec,Z,cfreq]=perfect_x(sol,omega);
    cfreq=cfreq/omega;
    freq=omega_s;
    freq2=cfreq/freq;
    Z=Z(abs(freq2)<3);
    freq2=freq2(abs(freq2)<3);
    plot3(i*ones(size(freq2)),freq2,log(abs(Z)),'-')
    maxi=max([maxi,max(log(abs(Z)))]);
    mini=min([mini,min(log(abs(Z)))]);
     ylim([ -3 3])
end
% plot_farey_3d(mini,maxi,ceil(xmax))
ylabel('$\omega/\Omega$','Interpreter','Latex');
title('fft of $x+ j y $','Interpreter','Latex');
view(90,00)


for figcount=1:length(figures)
    savefig(figures(figcount),[currDate,'/fig',num2str(figcount),'.fig'])
    saveas(figures(figcount),[currDate,'/fig',num2str(figcount),'.pdf'])
end

end

function [I_res,freq_vec,Spec]=perfect(ODESol_struct,Omega)
ts=ODESol_struct.x;
t_1=ts(end);
t_0=0.75*ts(end);

Fs = 200;            % Sampling frequency                    
T = 2*pi*Omega/Fs;    % Sampling period       
L = 10^5;             % Length of signal
t_sample = (t_0:T:t_1);
L=length(t_sample);
if mod(L,2)==1
t_sample=t_sample(1:end-1);
L=length(t_sample);
end
%dot(theta) is sixth component of W
dottheta = deval(ODESol_struct,t_sample,6);
% figure;hold on
% plot(t_sample,dottheta);
% plot(ts(ts>t_0),ODESol_struct.y(6,ts>t_0));

% Calculate Residual current
I_res = (Omega- dottheta).*sin(Omega*t_sample);
spec_a = fft(I_res);
spec_b = abs(spec_a/L); % This is the two-sided spectrum - absolute taken as this FS is in complex form.

% Construct the one-sided by taking the first half of the vector (this is
% the Nyquist frequency) and then doubling amplitudes.
Spec = spec_b(1:L/2+1);
Spec(2:end-1) = 2*Spec(2:end-1);

freq_vec = Fs*(0:(L/2))/L;
end


function [freq_vec,Spec,Z,cfreq]=perfect_x(ODESol_struct,Omega)
ts=ODESol_struct.x;
%t_0 = ts(end)*0.75; %Where to start calculating Fourier transform from
t_1=ts(end);
t_0=0.75*ts(end);

Fs = 200;            % Sampling frequency                    
T = 2*pi*Omega/Fs;    % Sampling period       
L = 10^5;             % Length of signal
t_sample = (t_0:T:t_1);
L=length(t_sample);
if mod(L,2)==1
t_sample=t_sample(1:end-1);
L=length(t_sample);
end  


%x is the first component of W
x = deval(ODESol_struct,t_sample,1);

% Calculate Residual current
spec_a = fft(x);
spec_b = abs(spec_a/L); % This is the two-sided spectrum - absolute taken as this FS is in complex form.

% Construct the one-sided by taking the first half of the vector (this is
% the Nyquist frequency) and then doubling amplitudes.
Spec = spec_b(1:L/2+1);
Spec(2:end-1) = 2*Spec(2:end-1);

freq_vec = Fs*(0:(L/2))/L;

% construct complex FFT 
y = deval(ODESol_struct,t_sample,2);
z=x+1j*y;
Z = fft(z)/L;
cfreq = Fs*[  (0:(L/2-1)) , (-L/2):-1 ]/L;
%sort into ascending frequency order
[cfreq,ix]=sort(cfreq);
Z=Z(ix);

end

function plot_farey(mini,maxi,xmax)
seq=farey_sequence(8);
seqs=seq;
for j=1:xmax
seqs=[seqs,seq+j];
end
% i tried to do farey sequence grid but it didn't seem to work
% ax=gca;
% grid minor
% ax.YAxis.MinorTick='on';
% ax.YGrid = 'on';
% ax.XGrid = 'off';
% ax.ZGrid = 'off';
% ax.YAxis.MinorTickValues=(unique(seqs,'sorted'));
plot([seqs;seqs],[maxi*ones(size(seqs));mini*ones(size(seqs))],'-','Color',[0.9,0.9,0.9]);
end

function plot_farey_3d(mini,maxi,xmax)
seq=farey_sequence(7);
seqs=seq;
for j=1:xmax
seqs=[seqs,seq+j];
end
% i tried to do farey sequence grid but it didn't seem to work
% ax=gca;
% grid minor
% ax.YAxis.MinorTick='on';
% ax.YGrid = 'on';
% ax.XGrid = 'off';
% ax.ZGrid = 'off';
% ax.YAxis.MinorTickValues=(unique(seqs,'sorted'));
plot3(zeros(size([seqs;seqs])),[seqs;seqs],[maxi*ones(size(seqs));mini*ones(size(seqs))],'-','Color',[0.9,0.9,0.9]);
end

function [parname,parvalues]=find_varying_parameter(pars,sols)
list_of_fields=fieldnames(pars);
found=0;
for i=1:length(list_of_fields)
    if ~ range(pars.(list_of_fields{i}))==0
        parname=list_of_fields{i};
        parvalues=pars.(list_of_fields{i});
        found=1;
    end
end
if found==0
        parname='ICs';
        parvalues=1:length(sols);
end
end

