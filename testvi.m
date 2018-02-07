alpha =1.;
x = 1.;
v = 0.;
dt = 1.03/64;
N = 64000
ti= zeros(1,N);
xi = zeros(1,N);
vi = zeros(1,N);

xi(1)  = x;
vi(1)  = v;
ti(1)   =dt;
% explicit
for(i=2:N)
  v = v-dt*alpha*x;
  x = x+dt*v;
  ti(i) = i*dt;
  xi(i) = x;
  vi(i)  = v;  
  %v = (xp -x)/dt
  %x = xp
endfor
xie = cos(alpha*ti);
vie = -alpha*sin(alpha*ti);
figure(1)
plot(ti,xie,ti,xi)
hold off
figure(2)
plot(ti, alpha*(xi.*xi), ti, vi.*vi , ti , alpha*xi.*xi+vi.*vi)
hold off
E =  alpha*xi.*xi+vi.*vi;
E(N)-E(1)

%implicit
x=1.;
v =0.;
for(i=2:N)
  xp = (x +dt*v)/(1.+dt*dt*alpha);
  v = (xp -x)/dt;
  x = xp ;
  ti(i) = i*dt;
  xi(i) = x;
  vi(i)  = v;  
  %v = (xp -x)/dt
  %x = xp
endfor
figure(3)
plot(ti,xie,ti,xi)
hold off
figure(4)
plot(ti, alpha*(xi.*xi), ti, vi.*vi , ti , alpha*xi.*xi+vi.*vi)
hold off
E =  alpha*xi.*xi+vi.*vi;
E(N)-E(1)

%crank nikolson
x=1.;
v =0.;
beta = 0.5;
for(i=2:N)
  xp = (x +dt*v - dt*dt*beta*alpha*x)/(1.+(1.-beta)*dt*dt*alpha);
  v = (xp -x)/dt;
  x = xp ;
  ti(i) = i*dt;
  xi(i) = x;
  vi(i)  = v;  
  %v = (xp -x)/dt
  %x = xp
endfor
figure(5)
plot(ti,xie,ti,xi)
hold off
figure(6)
plot(ti, alpha*(xi.*xi), ti, vi.*vi , ti , alpha*xi.*xi+vi.*vi)
hold off
E =  alpha*xi.*xi+vi.*vi;
E(N)-E(1)
%verlet
x=1.;
v =0.;
for(i=2:N)
  xp = x + v*dt-0.5*alpha*x*dt*dt;
  v = v+ (-alpha*x- alpha*xp)*0.5*dt;
  x = xp ;
  ti(i) = i*dt;
  xi(i) = x;
  vi(i)  = v;  
  %v = (xp -x)/dt
  %x = xp
endfor
figure(7)
plot(ti,xie,ti,xi)
hold off
figure(8)
plot(ti, alpha*(xi.*xi), ti, vi.*vi , ti , alpha*xi.*xi+vi.*vi)
E =  alpha*xi.*xi+vi.*vi;
E(N)-E(1)

hold off
