function save_data_figures(sols,pars,currDate)
[parname,parvalues]=find_varying_parameter(pars,sols);
xtix=cellstr(num2str((parvalues)))';



for i =1:length(sols)
    W=sols{i}.y;
    % how much would you like to plot? Last 0.1 of time?
    last=0.05;
    W=W(:,ceil(end*(1-last)):end);
    writetable(array2table([W',i*ones(max(size(W)),1)],'VariableNames',{'w1','w2','w3','w4','w5','w6','i'}),[currDate,'/dat1_',num2str(i),'.txt']);
end


for i =1:length(sols)
    W=sols{i}.y;
    % how much would you like to plot? Last 0.1 of time?
    last=0.01;
    W=W(:,ceil(end*(1-last)):end);
    % plot orbits in yz with parameter iterate in x axis
    writetable(array2table([W',i*ones(max(size(W)),1)+0.25*(W(4,:)')/max(W(4,:))],'VariableNames',{'w1','w2','w3','w4','w5','w6','it'}),[currDate,'/dat2_',num2str(i),'.txt']);
end


for i =1:length(sols)
    W=sols{i}.y;
    last=0.01;
    W=W(:,ceil(end*(1-last)):end);
    % transform to rotating frame
    xs=W(1,:).*cos(W(3,:))+W(2,:).*sin(W(3,:));
    ys=-W(1,:).*sin(W(3,:))+W(2,:).*cos(W(3,:));
    % as with previous plot but now in rotating frame
    writetable(array2table([xs',ys',i*ones(size(xs'))],'VariableNames',{'xt','yt','i'}),[currDate,'/dat3_',num2str(i),'.txt']);
end

% plots hi low bifurcation diagram

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
    writetable(array2table([Whi',ones(size(Whi'))],'VariableNames',{'hi','it'}),[currDate,'/dat4_hi',num2str(i),'.txt']);
    writetable(array2table([Wlo',ones(size(Wlo'))],'VariableNames',{'hi','it'}),[currDate,'/dat4_lo',num2str(i),'.txt']);

end



for i =1:length(sols)
    W=sols{i}.y;
    t=sols{i}.x;
    last=0.01;
    W=W(:,ceil(end*(1-last)):end);
    t=t(ceil(end*(1-last)):end);
    writetable(array2table([t',W(6,:)'],'VariableNames',{'t','thetadot'}),[currDate,'/dat5_',num2str(i),'.txt']);
    
    %plot([min(t),max(t)],omega_s*[1,1],'--')
end


for i =1:length(sols)
    We=sols{i}.ye;
    last=0.1;
    We=We(:,ceil(end*(1-last)):end);
    writetable(array2table([We',i*ones(size(We(1,:)'))],'VariableNames',{'we1','we2','we3','we4','we5','we6','it'}),[currDate,'/dat6_',num2str(i),'.txt']);
end


for i =1:length(sols)
    sol=sols{i};
    [omega,omega_s]=omegas(pars,i);
    [~,freq_vec,Spec]=perfect(sol,omega);
    freq_vec=freq_vec/omega;
    freq=omega_s;
    freq2=freq_vec/freq;
    Spec=Spec(freq2<5);
    freq2=freq2(freq2<5);
    writetable(array2table([freq2',Spec'],'VariableNames',{'freqs','specs'}),[currDate,'/dat7_',num2str(i),'.txt']);
end
for i =1:length(sols)
    sol=sols{i};
    [omega,omega_s]=omegas(pars,i);
    [freq_vec,Spec,Z,cfreq]=perfect_x(sol,omega);
    freq_vec=freq_vec/omega;
    freq=omega_s;
    freq2=freq_vec/freq;
    Spec=Spec(freq2<5);
    freq2=freq2(freq2<5);
    writetable(array2table([freq2',Spec'],'VariableNames',{'freqs','specs'}),[currDate,'/dat8_',num2str(i),'.txt']);
end

for i =1:length(sols)
    sol=sols{i};
    [omega,omega_s]=omegas(pars,i);
    [freq_vec,Spec,Z,cfreq]=perfect_x(sol,omega);
    cfreq=cfreq/omega;
    freq=omega_s;
    freq2=cfreq/freq;
    Z=Z(abs(freq2)<3);
    freq2=freq2(abs(freq2)<3);
    writetable(array2table([freq2',abs(Z)'],'VariableNames',{'freqs','specs'}),[currDate,'/dat9_',num2str(i),'.txt']);
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

