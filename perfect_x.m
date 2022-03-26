function [freq_vec,Spec,Z,cfreq,x,t_sample]=perfect_x(ODESol_struct,Omega)
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
