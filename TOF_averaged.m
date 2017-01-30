function [ ] = TOF_averaged( )
%TOF_averaged( )
%   Adds up the TOF for many shots, to see trap freq etc

bin_size = 0.0001;
tmax = 0.55;
tmin = 0.4;
binned_time_tot = zeros(floor((tmax-tmin)/bin_size),1);

for file = 1:50
    s = num2str(file);
    filename_call = ['z:\release\9_8_trap_freq_AL2',s];
    three_ch_out = dld_read_5channels_reconst_multi(filename_call,1,1,0);
    [binned_time,T] = dld_tof_plotter(three_ch_out,tmin,tmax,888,bin_size,0);
    binned_time_tot = binned_time_tot + binned_time;
    
end

plot(binned_time_tot)

%%%%%% FFT

%fo = 4;   %frequency of the cos wave
Fs = ceil(1/bin_size); %sampling rate
%Ts = 1/Fs; %sampling time interval

t = [bin_size:bin_size:tmax-tmin];
%t = 0:Ts:1-Ts; %sampling period
n = length(t); %number of samples
%y = 2*sin(2*pi*fo*t); %the cosine curve
y = binned_time_tot;

%plot the cosine curve in the time domain
cosPlot = figure;
plot(t,y)
xlabel('time (seconds)')
ylabel('y(t)')
title('Normalised input')
%set(cosPlot,'Position',[500,500,500,300])
grid

[YfreqDomain,frequencyRange] = positiveFFT(y,Fs);
positiveFFTplot = figure;
plot_length = length(frequencyRange);
plot(frequencyRange(2:plot_length),abs(YfreqDomain(2:plot_length)));
%set(positiveFFT,'Position',[500,500,500,300])
xlabel('Freq (Hz)')
ylabel('Amplitude')
title('Using the positiveFFT function')
grid
%axis([0,20,0,1.5])


%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P,f,len]=power_spec(binned_time_tot,t,1e4,1,0,5000);

end