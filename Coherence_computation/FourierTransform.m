function zM = FourierTransform(xV,nsur)
% zM = FourierTransform(xV,nsur)
% FOURIERTRANSFORM generates a given number of Fourier Transform (FT) 
% surrogates for the given time series 'xV'. FT suggests that the phases of
% the Fourier transform of 'xV' are randomized and the inverse Fourier 
% transform gives the surrogate time series.
% INPUT
% - xV  : the given time series
% - nsur: the number of surrogate time series (default is 1)
% OUTPUT
% - zM  : the n x nsur matrix of 'nsur' FT surrogate time series
%========================================================================
%     <FourierTransform.m>, v 1.0 2010/02/11 22:09:14  Kugiumtzis & Tsimpiris
%     This is part of the MATS-Toolkit http://eeganalysis.web.auth.gr/

%========================================================================
% Copyright (C) 2010 by Dimitris Kugiumtzis and Alkiviadis Tsimpiris 
%                       <dkugiu@gen.auth.gr>

%========================================================================
% Version: 1.0

% The FreeBSD Copyright:	
% Copyright 1992-2012 The FreeBSD Project. All rights reserved.	

% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:	

% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.	
% Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.	

% "THIS SOFTWARE IS PROVIDED BY THE FREEBSD PROJECT ``AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE FREEBSD PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
% OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
% IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."

% The views and conclusions contained in the software and documentation are those of the authors and should not
% be interpreted as representing official policies, either expressed or implied, of the FreeBSD Project.

%=========================================================================
% Reference : D. Kugiumtzis and A. Tsimpiris, "Measures of Analysis of Time Series (MATS): 
% 	          A Matlab  Toolkit for Computation of Multiple Measures on Time Series Data Bases",
%             Journal of Statistical Software, Vol. 33, Issue 5, 2010
% Link      : http://eeganalysis.web.auth.gr/
%========================================================================= 

 
if nargin==1
    nsur=1;
end
n=length(xV);
if rem(n,2) == 0
    n2 = n/2;
else
    n2 = (n-1)/2;
end
nfreq = n2;  
fxV = fft(xV,2*nfreq);
magnV = abs(fxV);
% The magnitudes of the Fourier transform
mV = [magnV(1:nfreq+1)' flipud(magnV(2:nfreq))']';
fiV = angle(fxV);
zM = NaN*ones(n,nsur);
for isur=1:nsur
    % The random phases
    rfiV = rand(nfreq-1,1) * 2 * pi;
    nfiV = [0; rfiV; fiV(nfreq+1); -flipud(rfiV)];
    fzV = mV .* exp(nfiV .* i); 
    % The inverse transform
    zM(:,isur)=real(ifft(fzV,n));
end
