%program to fit a cloud of atoms below the critical temperature.
% it fits the thermal part to an Integrated density, i.e. as we see it
% on our electron multiplier.  This expression is integrated in two
% dimensions, just as our detector does.
% The fitting procedure is as follows, the user estimates the TF radius,
% i.e the half-width of the bottom of the visible condensate.  The program
% then fits the thermal tail and subtracts this from the rest of the data,
% this subtracted curve, the condensate if you like, is then fit.  This
% whole procedure is then itterated so that a best fit for both curves is
% found.  Thus the TF radius is a fitted value, using the users first guess
% only.

%%%%%ONLY GOOD IN TF APPROXIMATION%%%%%%%%%%%

%%%%%%%%%%%%%%%inputs%%%%%%%%%%%%%%%%%%%%%%%%
% z - fall distance in metres
% tau - fall time in seconds
% Thomas Fermi radius - in seconds
% fr - radial trapping frequency in Hertz (should be the higher one)
% fa - axial trapping frequency in Hertz
% nmax - number of thermal states to sum over (200)
% t - time array (1-D) (in seconds) NOT USED
% d - data array (1-D) NOTUSED
% plot_i_data - flag (either 1 or 0) if 1 plots initial data
% plot_f_res - flag (either 1 or 0) if 1 plots fitted data
% time_min - minimum time in seconds
% time_max - minimum time in seconds
% bin_time - time bin in seconds
% three_channel_output_sorted - dld data in form of t,x,y matrix

%%%%%%%%%%%%%%%outputs%%%%%%%%%%%%%%%%%%%%%%%%
% t - same as raw time data
% rw - normalised raw data
% y - is a curve produced with fitting parameters
% thermID - fitted thermal profile
% condID - fitted condensate profile
% T - fitted Temperature (in Kelvin)
% Tv0 - Temperature based on FWHM = 2v0/g (only approx below Tc)
% TFr - Thomas Fermi Radius in seconds
% thermfrac - Thermal fraction obtained from fit
% condfrac - Condensate fraction obtained from fit
% mu - chemical potential in joules obtained from the fit
% N01 - Number of atoms in condensate obtained from the fit
% N02 - Number of atoms in condensate obtained from using the critcal
%       temperature relationship
% Ntherm - Number of thermal atoms obtained from using the critcal
%       temperature relationship
function[t,rw,y,thermID,condID,T,Tv0,TFr,thermfrac,condfrac,mu,N01,N02,Ntherm]=...
                                TC_and_below_a(z,tau,TF_radius,fr,fa,nmax,plot_i_data,plot_f_res,t_min,t_max,bin_size,three_channel_output_sorted) ; 

close all

[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;

t = t_min:bin_size:t_max-bin_size;

binned_time = zeros(floor((t_max-t_min)/bin_size),1);

count6a = size(three_channel_output_sorted);
count6 = count6a(1);

max_bin = floor(t_max/bin_size);
min_bin = ceil(t_min/bin_size);

for count7=1:count6
    bin_for_I_plot = ceil((three_channel_output_sorted(count7,1))/bin_size)-min_bin;
    if bin_for_I_plot > 0 && bin_for_I_plot < max_bin - min_bin
        d(bin_for_I_plot) = binned_time(bin_for_I_plot) + 1;
    end
    
end

min_index = min(length(t),length(d));
t = t(1:min_index);
d = d(1:min_index);

if plot_i_data ==1
    figure(150)
    plot(t,d)
    rw = 0;
    y = 0;
    thermID = 0;
    condID = 0;
    T = 0;
    Tv0 = 0;
    TFr = 0;
    thermfrac = 0;
    condfrac = 0;
    mu = 0;
    N01 = 0;
    N02 = 0;
    Ntherm = 0;
%     return
end

d=d/max(d); % normalise data

%%%%%Lets find the peak in the data%%%%%%
[a,pk] = max(d);
t=t-t(pk); % centre peak
%%%%%%estimate half width for temperature guess%%%%
ll=min(d); % minimum of data

%%%%%%fit parameters%%%%
bg=0;
T=1e-6;
Amptherm=.5;
Ampcond=.5;
t0=0;
fp(1)=bg;
fp(2)=T;
fp(3)=Amptherm;
fp(4)=Ampcond;
fp(5)=TF_radius;
fp(6)=t0;
vel=g*tau;
wr=2*pi*fr;
wa=2*pi*fa;
vel=g*tau;

OPTIONS = OPTIMSET('MaxFunEvals',20000,'MaxIter',20000);


fitt=fminsearch('TC_and_below_tof',fp,OPTIONS,t,d,vel,tau,wr,nmax)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%plot out results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Thermal profile%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0=sqrt(((2*kb*fitt(2))/(Hemass*wr^2))*(1+(wr^2)*tau^2));
n=1:1:nmax;
t=t+fitt(6);
zz=vel*t; %convert flight into distance

Thermtot=0;
for i=1:length(t)
    ID(i)=fitt(3)*sum(((exp(-(zz(i)^2/z0^2))).^n)./n.^(5/2));
end
thermID=ID;

[Tfw,tt,rww,yy,Tgf]=tof_gauss(t,thermID,0,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Condensate profile%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(t)
    if t(i)>-fitt(5) && t(i)<fitt(5)
ID(i)=ID(i)+(fitt(4)*((1-zz(i)^2/(fitt(5)*vel)^2)^2));
    end
end
condID=ID-thermID;
y=ID+fitt(1);

if plot_f_res ==1
plot(t,d,'ob',t,y,'-r')
xlabel('Time (seconds)')
end
rw=d;
tot=sum(condID)+sum(thermID);
T=fitt(2);
Tv0=Tgf;%use 1/e fitted width
TFr =fitt(5);
thermfrac=sum(thermID)/tot;
condfrac=sum(condID)/tot;
mu=0.5*Hemass*wr^2*(TFr*vel)^2/(1+(wr^2*tau^2));
mu=mu/(2*pi*hbar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%working out N0 from Chem. pot%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wbar=(wa*wr^2)^(1/3);
aho=sqrt(hbar/(Hemass*wbar));
a=7.512e-9;

cp=0.5*Hemass*(wr^2/(1+wr^2*tau^2))*(TFr*vel)^2
N01=((2*cp/(hbar*wbar))^(5/2))*aho/(15*a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%working out thermal fraction from critical temp%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ntherm=1.202*(kb*T/(hbar*wbar))^3
total=Ntherm/thermfrac;
N02=total*condfrac
a=((2*cp/(hbar*wbar))^(5/2))*aho/(15*N02)

































