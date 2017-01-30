function err= guassfun_b(p,s,c)

B=p(1);
xc=p(2);
w=p(3);
bg=p(4);

v=B*exp(-((2*(s-xc))/w).^2)+bg;
v=v-c;
err=sum(v.^2);