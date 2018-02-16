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



% explicit (Euler)Cromer f34
fi = fi+2
x = 1.;
v = 0.;
xi(1)  = x;
vi(1)  = v;
for(i=2:N)
  v = v-dt*alpha*x;
  x = x+dt*v;
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
printf("explicite Euler cromer , figures %d, %d \n", fi, fi+1)
printf("error x: %d, error v: %d , errorE: %d \n", abs(xie(N)-xi(N)), abs(vie(N)-vi(N)), E(N)-Ee(N) )


% midpoint method f56
fi = fi+2;
x = 1.;
v = 0.;
xi(1)  = x;
vi(1)  = v;
for(i=2:N)
  vp = v-dt*alpha*x;
  x = x+0.5*dt*(vp+v);
  v= vp;
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
printf("midpoint , figures %d, %d \n", fi, fi+1)
printf("error x: %d, error v: %d , errorE: %d \n", abs(xie(N)-xi(N)), abs(vie(N)-vi(N)), E(N)-Ee(N) )

% half step method f78 second order in time for x, v et E, but E increas (slightly)
fi = fi+2;
x = 1.;
v = 0.;
vhalf = v -0.5*alpha*x*dt;
xi(1)  = x;
vi(1)  = v;
for(i=2:N)
  vhalfp = vhalf - alpha*x*dt;
  xp     = x + vhalfp*dt;
  
  vhalfpp = vhalfp - alpha*xp*dt;
  vhalf = vhalfp;
  x= xp;
  v = 0.5*(vhalfp+vhalfpp);
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
printf("half step , figures %d, %d \n", fi, fi+1)
printf("error x: %d, error v: %d , errorE: %d \n", abs(xie(N)-xi(N)), abs(vie(N)-vi(N)), E(N)-Ee(N))


%implicit
fi = fi+2;
x=1.;
v =0.;
xi(1)  = x;
vi(1)  = v;
for(i=2:N)
  xp = (x +dt*v)/(1.+dt*dt*alpha);
  v = (xp -x)/dt;
  x = xp ;
  xi(i) = x;
  vi(i)  = v;  
  %v = (xp -x)/dt
  %x = xp
endfor
E =  alpha*xi.*xi+vi.*vi;
figure(fi)
plot(tie,xiee,ti,xi,"-o")
hold off
figure(fi+1)
plot(ti, alpha*(xi.*xi), ti, vi.*vi , ti , E)
hold off

printf("implicit , figures %d, %d \n", fi, fi+1)
printf("error x: %d, error v: %d , errorE: %d, \n", abs(xie(N)-xi(N)), abs(vie(N)-vi(N)), E(N)-Ee(N))

%crank nikolson
fi = fi+2
x=1.;
v =0.;
xi(1)  = x;
vi(1)  = v;
beta = 0.5;
for(i=2:N)
  xp = (x +dt*v - dt*dt*beta*alpha*x)/(1.+(1.-beta)*dt*dt*alpha);
  v = (xp -x)/dt;
  x = xp ;
  xi(i) = x;
  vi(i)  = v;  
  %v = (xp -x)/dt
  %x = xp
endfor
E =  alpha*xi.*xi+vi.*vi;
figure(fi)
plot(tie,xiee,ti,xi,"-o")
title("implicit")
hold off
figure(fi+1)

plot(ti, alpha*(xi.*xi), ti, vi.*vi , ti , E)
title("implicit")

hold off
printf("implicit , figures %d, %d \n", fi, fi+1)
printf("error x: %d, error v: %d , errorE: %d \n", abs(xie(N)-xi(N)), abs(vie(N)-vi(N)), E(N)-Ee(N))

%verlet
fi =fi+2
x=1.;
v =0.;
xi(1)  = x;
vi(1)  = v;
for(i=2:N)
  xp = x + v*dt-0.5*alpha*x*dt*dt;
  v = v+ (-alpha*x- alpha*xp)*0.5*dt;
  x = xp ;
  xi(i) = x;
  vi(i)  = v;  
  %v = (xp -x)/dt
  %x = xp
endfor
E =  alpha*xi.*xi+vi.*vi;
figure(fi)
plot(tie,xiee,ti,xi,"-o")
title("velocity verlet");
hold off
figure(fi+1)
plot(ti, alpha*(xi.*xi), ti, vi.*vi , ti , E)
title("velocity verlet");

printf("velocity verlet , figures %d, %d \n", fi, fi+1)
printf("error x: %d, error v: %d , errorE: %d \n", abs(xie(N)-xi(N)), abs(vie(N)-vi(N)), E(N)-Ee(N))

% all firt order for x and v and E 
%but verlet second order for E, but less precise for x (error 2 time bigger)

%Euler Richardson
fi =fi+2
x=1.;
v =0.;
xi(1)  = x;
vi(1)  = v;
for(i=2:N)
  xmid = x +0.5*dt*alpha*v;
  vmid = v -0.5*dt*alpha*x;
  amid = -alpha*xmid;
  v    = v + amid*dt;
  x    = x + vmid*dt;
  xi(i) = x;
  vi(i)  = v;  
  %v = (xp -x)/dt
  %x = xp
endfor
E =  alpha*xi.*xi+vi.*vi;
figure(fi)
plot(tie,xiee,ti,xi,"-o")
title("euler richardson")
hold off
figure(fi+1)
title("euler richardson")
plot(ti, alpha*(xi.*xi), ti, vi.*vi , ti , E)
printf("euler richardson , figures %d, %d \n", fi, fi+1)
printf("error x: %d, error v: %d , errorE: %d \n", abs(xie(N)-xi(N)), abs(vie(N)-vi(N)), E(N)-Ee(N))

%Beeman
fi =fi+2
x=1.;
v =0.;
xi(1)  = x;
vi(1)  = v;
% VErlet for the first step
a = -alpha*x;
xp = x + v*dt+0.5*a*dt*dt;
ap = -alpha*xp;
vp = v +0.5*(ap+a)*dt;
am = a;
a = ap;
x = xp;
v = vp;

xi(2)  = x;
vi(2)  = v;

for(i=3:N)
  xp = x + v*dt+(4*a -am)*dt*dt/6.;
  ap = -alpha*xp;
  vp = v+ (2*ap+5*a-am)*dt/6.;
  am = a;
  a = ap;
  v =vp;
  x =xp;
  xi(i) = x;
  vi(i)  = v;  
  %v = (xp -x)/dt
  %x = xp
endfor
E =  alpha*xi.*xi+vi.*vi;
figure(fi)
plot(tie,xiee,ti,xi,"-o")
hold off
figure(fi+1)
plot(ti, alpha*(xi.*xi), ti, vi.*vi , ti , E)
printf("Beeman , figures %d, %d \n", fi, fi+1)
printf("error x: %d, error v: %d , errorE: %d \n", abs(xie(N)-xi(N)), abs(vie(N)-vi(N)), E(N)-Ee(N))

% all firt order for x and v and E 
%but verlet second order for E, but less precise for x (error 2 time bigger)
hold off
