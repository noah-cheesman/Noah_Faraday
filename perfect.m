function [I_res,freq_vec,Spec,t_sample]=perfect(ODESol_struct,Omega)
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
