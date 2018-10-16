% TOF fitting function
% here we just fit a gaussian to the data and find the FWHM
% then as in Yavin et. al we assine the FWHM = 2v0/g
%%%%%%%%%%%%%%%inputs%%%%%%%%%%%%%%%%%%%%%%%%
% t - time array (1-D) (in seconds)
% d - data array (1-D)
% plot_initial_data - flag (either 1 or 0) if 1 plots initial data
% plot_final_results - flag (either 1 or 0) if 1 plots fitted data

%%%%%%%%%%%%%%%outputs%%%%%%%%%%%%%%%%%%%%%%%%
% Tfw - Temperature based on FW1_e =2v0/g, this was the guess used
% t - same as raw time data
% rw - normalised raw data
% y - is a curve produced with fitting parameters
% background - fitted background
% T - fitted Temperature (in Kelvin)
%
function[Tfw,t,rw,y,T]=tof_gauss_dld_a(t,d,plot_i_data,plot_f_res,fig_number)

% time_intensity_array=dlmread('pulse_TOF_2d.txt',',');
% 
% t = time_intensity_array(:,1);
% d =  time_intensity_array(:,2);

%close all
[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;

size(t)
size(d)

if plot_i_data ==1
    figure(fig_number)
plot(t,d)
end

d=d/max(d); % normalise data

%%%%%Lets find the peak in the data%%%%%%
[a,pk] = max(d);
t2=t-t(pk); % centre peak
%%%%%%estimate half width for temperature guess%%%%
ll=min(d); % minimum of data

fw1=0;
for i =1:length(t)
    if d(i)> exp(-1)
        fw1=t2(i);
        j=i;
        break
    end
end

fw2=0;
for i =j+3:length(t)
    if d(i)< exp(-1)
         fw2=t2(i);
        break
    end
end

% g=9.8;
% Hemass=4e-27;
% kb=1.38e-23;

FW1_e=(abs(fw1)+abs(fw2)); % approximate FWHM
v0=g*FW1_e/2;
Tfw=(0.5*Hemass/kb)*v0^2;

B=1;
t0=t(pk);
w=FW1_e;
bg=0;
fp(1)=B;
fp(2)=t0;
fp(3)=w;
fp(4)=bg;
fitt=fminsearch('guassfun_a',fp,[],t,d); 
y=fitt(1)*exp(-((2*(t-fitt(2)))/fitt(3)).^2)+fitt(4);

w1_e=fitt(3);
v0=g*w1_e/2;
T=Hemass*v0^2/(2*kb);
if plot_f_res == 1
plot(t,d,'b.',t,y,'r')
xlabel('Time')
end

rw=d;