% function that computes the power spectrum of data
%inputs: 
% amplitude - array of amplitude data
% t - array of time data
% Nfft- length of data that the FFT is averaged over (higher number
%           gives a higher resolution
% len - length of psd
% fstart - start freq. for plot
% fstop - stop freq. for plot

function [P,f,len]=power_spec(amplitude,t,Nfft,plt,fstart,fstop);

dt=t(2)-t(1);
fs=1/dt;

amplitude=amplitude-mean(amplitude);

[P,f]= psd(amplitude,Nfft,fs,[]);

rr=1;
hh=1;
hh=length(f);

    for i =1:length(f)
        if f(i)>fstart 
            rr=i;
            break
        end
    end
      for j =rr:length(f)
        if f(j)>fstop
            hh=j;
            break
        end
    end
if plt ==1
    plot(f(rr:hh),P(rr:hh))
    xlabel('f [Hz]')
end

f=f(rr:hh)';
P=P(rr:hh)';

len=length(P);





