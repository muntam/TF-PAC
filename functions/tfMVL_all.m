function [tf_MVL,tf_PLV,tf_MI]=tfMVL_all(x,high_freq,low_freq, Fs)
% This function computes the phase amplitude coupling using TF-MVL, TF-PLV and TF-MI methods.
% Input:   x            : input signal 
%          high_freq    : Amplitude Frequency range 
%          low_freq     : Phase Frequency range 
%          Fs           : Sampling Frequency
% Output:  tf_MVL   : Computed PAC using TF-MVL method
%          tf_PLV   : Computed PAC using TF-PLV method
%          tf_MI    : Computed PAC using TF-MI method

% Written by: Tamanna T. K. Munia, January 2019

% These scripts have been optimised for the Windows operating systm  
% MATLAB version used 2018a.


%% Amplitude and Phase calculation
[tfd]=rid_rihaczek4(x,Fs);
W=tfd;
W2=W(2:end,:);
tfd_low=W2(low_freq:low_freq,:);
angle_low=angle(tfd_low);
Amp=((W2(high_freq:high_freq,:)));
Phase=((angle_low));
tf_MVL =(calc_MVL(Phase,Amp));
tf_PLV  = (calc_PLV(Phase,Amp));
tf_MI = (calc_MI(Phase,Amp));
 
end

%% Functions to compute PAC
% Using the metrics MLV(Canolty et al.,2006), PLV (Cohen 2008). and MI (Tort et al.,2010)
% Code adapted from sensory_PAC_master written by Robert Seymour - Aston Brain Centre. July 2017.

function [MVL] = calc_MVL(Phase,Amp)
         Amp=(Amp);     
         z1=(exp(1i*Phase));
         z=Amp.*(z1);% Generate complex valued signal
         MVL = abs((mean(z)));
end

function [PLV] = calc_PLV(Phase,Amp)
         Amp=(Amp);
         amp_phase = angle(((Amp))); % Phase of amplitude envelope
         PLV = abs(mean((exp(1i*(Phase-amp_phase)))));
end

function [MI] = calc_MI(Phase,Amp)
         nbin=18; 
         position=zeros(1,nbin); 
         winsize = 2*pi/nbin;
        for j=1:nbin
            position(j) = -pi+(j-1)*winsize;
        end
        
        % now we compute the mean amplitude in each phase:
        MeanAmp=zeros(1,nbin);
        for j=1:nbin
            I = find(Phase <  position(j)+winsize & Phase >=  position(j));
            MeanAmp(j)=mean(Amp(I));
        end       
        MI=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);
        MI=abs(MI);
end


