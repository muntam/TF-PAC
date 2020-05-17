function [tfd]=rid_rihaczek4(x,fbins)
% This function computes reduced interference Rihaczek distribution;
% Input: x: signal, fbins=required frequency bins
% Output: tfd = Generated reduced interference Rihaczek distribution

tbins = length(x);
amb = zeros(tbins);

for tau = 1:tbins
    amb(tau,:) = (conj(x) .* x([tau:tbins 1:tau-1]) );
end

ambTemp = [amb(:,tbins/2+1:tbins) amb(:,1:tbins/2)];
amb1 = [ambTemp(tbins/2+1:tbins,:); ambTemp(1:tbins/2,:)];

D=(-1:2/(tbins-1):1)'*(-1:2/(tbins-1):1);
L=D;
K=chwi_krn(D,L,0.01);
[s,d]=size(amb1);
df=K(1:s,1:d);
ambf = amb1 .* df;

A = zeros(fbins,tbins);
tbins=tbins-1;
if tbins ~= fbins
    for tt = 1:tbins
        A(:,tt) = datawrap(ambf(:,tt), fbins);
    end
else
    A = ambf; 
end

tfd = fft(A);

function K=chwi_krn(D,L,A)
%CHWI_KRN Choi-Williams kernel function.
%   K = CHWI_KRN(D,L,A) returns the values K of the Choi-Williams kernel
%   function evaluated at the doppler-values in matrix D and the lag-
%   values in matrix L. Matrices D and L must have the same size. The  
%   values in D should be in the range between -1 and +1 (with +1 being
%   the Nyquist frequency). The parameter A is optional and controls the
%   "diagonal bandwidth" of the kernel. Matrix K is of the same size as
%   the matrices D and L. Parameter A defaults to 10 if omitted.

%   Copyright (c) 1998 by Robert M. Nickel
%   $Revision: 1.1.1.1 $
%   $Date: 2001/03/05 09:09:36 $

if nargin<3; A=[]; end
if isempty(A); A=10; end
K=exp((-1/(A*A))*(D.*D.*L.*L));
end
end
