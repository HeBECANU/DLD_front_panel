function err= TC_and_below_tof_a(p,s,c,vel,tau,wr,nmax)

bg=p(1);
T=p(2);
Amptherm=p(3);
Ampcond=p(4);
TFradius=p(5);
t0=p(6);


[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Thermal profile%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0=sqrt(((2*kb*T)/(Hemass*wr^2))*(1+(wr^2)*tau^2));

n=1:1:nmax;
s=s+t0;
zz=vel*s; %convert flight into distance

for i=1:length(c)
ID(i)=Amptherm*sum(((exp(-(zz(i)^2/z0^2))).^n)./n.^(5/2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Condensate profile%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(c)
    if s(i)>-TFradius && s(i)<TFradius
ID(i)=ID(i)+(Ampcond*((1-zz(i)^2/(TFradius*vel)^2)^2));
    end
end

y=ID+bg;

yy=y-c; % remove ' for labview!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


err=sum(yy.^2);

