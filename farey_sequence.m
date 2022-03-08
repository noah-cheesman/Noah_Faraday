function [out]=farey_sequence(n)
% algorithm to find farey sequence of order n
  a=0;b=1;c=1;d=n;
  % Because we will only consider smallish n  (<10), we need not worry
  % about pre-allocation
  out=[];
  % add 0 to the sequence
  out=[out,a/b];
  while c <= n 
    k = ((n+b)-mod(n+b,d))/d;
    % temporary holder
    tempo = [c,d,c*k-a,d*k-b];
    % update values for
    a=tempo(1);b=tempo(2);c=tempo(3);d=tempo(4);
    % add next term in sequence
    out=[out,a/b];
  end
end