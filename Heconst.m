%%%% function returns all the parameters relevant to metastable helium
function[Hemass,Hegamma,Helambda,Helife,HeIs,Hemu,hbar,kb,Hek,g]=Heconst(); 
     

Na=6e23;
Hemass=4/Na/1000;
Hegamma=1.62e6*2*pi;
Helambda=1083.33e-9;
Helife=1/(Hegamma*2*pi);
HeIs=0.17;
Hemu=2.8e6;
hbar=1.05457e-34;
kb=1.38066e-23;
Hek=(2*pi)/Helambda;
g=9.8;