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

% pred corr fi 12
fi = 1
x = 1.;
v = 0.;
xi(1)  = x;
vi(1)  = v;
for(i=2:N)
  xp = x + dt *v;
  fx = -alpha*x;
  fxp = -alpha*xp;
  vp = v +dt*fx + dt*0.5*(fxp-fx);
  printf("pred %d \n", vp) 
  cond = true;
  while(cond)
    xpp = xp;
    xp = x +dt*v+dt*0.5*(vp-v);
    vpp = vp;
    fxp = -alpha*xp;
    vp = v +dt*fx+dt*0.5*(fxp-fx);
    printf("corr x %d %d \n", xp, abs( xpp-xp))
    printf("corr  v %d %d \n", vp, abs( vpp-vp))
    cond = abs( xpp-xp) > 1.e-6; 
  end
  x = xp;
  v = vp;
  xi(i) = x;
  vi(i)  = v; 
 %pause(-1) 
endfor
E =  alpha*xi.*xi+vi.*vi;

figure(fi)
plot(tie,xiee,ti,xi,"-o")
hold off
figure(fi+1)
plot(ti, alpha*(xi.*xi), ti, vi.*vi , ti , E)
hold off
printf("pred corr, figures %d, %d \n", fi, fi+1)
printf("error x: %d, error v: %d , errorE: %d \n", abs(xie(N)-xi(N)), abs(vie(N)-vi(N)), E(N)-Ee(N) )