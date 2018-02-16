alpha =1.;
x = 1.;
v = 0.;

Nstep = 10;
period = 2*pi*alpha;

dt = period/(Nstep)
N=Nstep*1+1;
ti= zeros(1,N);
for(i=1:N)
  ti(i) = (i-1)*dt;
endfor
tie= zeros(1,N*10);
for(i=1:N*10)
  tie(i) = (i-1)*dt/10.;
endfor
xi = zeros(1,N);
vi = zeros(1,N);
xiee = cos(alpha*tie);
viee = -alpha*sin(alpha*tie);
xie = cos(alpha*ti);
vie = -alpha*sin(alpha*ti);
Ee =  alpha*xie.*xie+vie.*vie;

% explicit (Euler) fi 12
fi = 1.
x = 1.;
v = 0.;
xi(1)  = x;
vi(1)  = v;
for(i=2:N)
  vp = v-dt*alpha*x;
  x = x+dt*v;
  v = vp;
  xi(i) = x;
  vi(i)  = v;  
endfor
E =  alpha*xi.*xi+vi.*vi;

figure(fi)
plot(tie,xiee,ti,xi,"-o")
hold off
figure(fi+1)
plot(ti, alpha*(xi.*xi), ti, vi.*vi , ti , E)
hold off
printf("explicite Euler, figures %d, %d \n", fi, fi+1)
printf("error x: %d, error v: %d , errorE: %d \n", abs(xie(N)-xi(N)), abs(vie(N)-vi(N)), E(N)-Ee(N) )

