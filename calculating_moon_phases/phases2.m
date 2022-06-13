format long;
phases=csvread("phases.csv");
t = phases(:,1);
phase=phases(:,2);
poly4=polyfit(t,phase,4) ## Or I could use the value of D below?
phase2 = phase - polyval(poly4,t);
mp = @(T) polyval([-6.797e-8,14.348e-6,89.97e-4,477198.86763133,134.96341138]/180*pi,T/36525)
D = @(T) polyval([-0.884e-8,1.832e-6,-16.3e-4,445267.11151675,297.8502042]/180*pi,T/36525)
F = @(T) polyval ([0.116e-8,-0.284e-6,-34.029e-4,483202.01752731,93.27209932]/180*pi,T/36525)
M = @(T) polyval ([0.041e-6,-1.536e-4,35999.05029094,357.52910918]/180*pi,T/36525)
dt=t(2)-t(1)
n = rows(t)
freq=linspace(0,1/dt,n);

fnc1 = @(p) sumsq(phase2-p(1)*sin(mp(t)))
p1=[0.109764];
[xmin1,fval1]=fminsearch(fnc1,p1)
phase2 = phase2 - xmin1 * sin(mp(t));

fnc1 = @(p) sumsq(phase2-p(1)*sin(M(t)))
p1=[-0.036652];
[xmin2,fval2]=fminsearch(fnc1,p1)
phase2 = phase2 - xmin2 * sin(M(t));

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*D(t)-mp(t)))
p1=[0.022235 ];
[xmin3,fval3]=fminsearch(fnc1,p1)
phase2 = phase2 - xmin3*sin(2*D(t)-mp(t));


fnc1 = @(p) sumsq(phase2-(ifelse(mod(D(t),pi) <= p(1)/2 | mod(D(t),pi) >= pi-p(1)/2, ...
  p(2)*mod(D(t)+p(1)/2,pi)/p(1), p(2)*(1-mod(D(t)-p(1)/2,pi)/(pi-p(1))))-p(2)/2)-p(3)*sin(2*D(t)))
p1=[1.6191e-01,   2.5418e-02,   6.5894e-03];
[xmin4,fval4]=fminsearch(fnc1,p1)
phase2 = phase2 - (ifelse(mod(D(t),pi) <= xmin4(1)/2 | mod(D(t),pi) >= pi-xmin4(1)/2, ...
  xmin4(2)*mod(D(t)+xmin4(1)/2,pi)/xmin4(1), xmin4(2)*(1-mod(D(t)-xmin4(1)/2,pi)/ ...
  (pi-xmin4(1))))-xmin4(2)/2) - xmin4(3) * sin(2*D(t));

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*mp(t)))
p1=[0.003735   ];
[xmin5,fval5]=fminsearch(fnc1,p1)
phase2 = phase2 - xmin5 * sin(2*mp(t));


fnc1 = @(p) sumsq(phase2-p(1)*sin(D(t)))
p1=[0.002  ];
[xmin10,fval10]=fminsearch(fnc1,p1)
phase2 = phase2 - xmin10 * sin(D(t));

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*(F(t)-D(t))))
p1=[0.0017  ];
[xmin11,fval11]=fminsearch(fnc1,p1)
phase2 = phase2 - xmin11 * sin(2*(F(t)-D(t)));


pfft=fft(phase2);
[maxv,maxi]=max(abs(pfft))
freq(maxi)*360*36525

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*(F(t)-3*D(t))))
p1=[0.0016  ];
[xmin12,fval12]=fminsearch(fnc1,p1)
phase2 = phase2 - xmin12 * sin(2*(F(t)-3*D(t)));

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*F(t)))
p1=[-0.0017  ];
[xmin13,fval13]=fminsearch(fnc1,p1)
phase2 = phase2 - xmin13 * sin(2*F(t));

pfft=fft(phase2);

actualphase = (1-cos(phase))/2;
myangle=(polyval(poly4,t) + xmin1 * sin(mp(t)) + xmin2 * sin(M(t)) + xmin3*sin(2*D(t)-mp(t)) + xmin4(3) * sin(2*D(t)) ...
  + (ifelse(mod(D(t),pi) <= xmin4(1)/2 | mod(D(t),pi) >= pi-xmin4(1)/2,xmin4(2)*mod(D(t)+xmin4(1)/2,pi)/xmin4(1), ...
    xmin4(2)*(1-mod(D(t)-xmin4(1)/2,pi)/(pi-xmin4(1))))-xmin4(2)/2)
  + xmin5 * sin(2*mp(t)) + xmin10 * sin(D(t)) + xmin11*sin(2*(F(t)-D(t))) ...
  + xmin12 * sin(2*(F(t)-3*D(t))) + xmin13 * sin(2*F(t)));

max(abs(phase-myangle))*180/pi
mean(abs(phase-myangle))*180/pi
median(abs(phase-myangle))*180/pi
sqrt(meansq(phase-myangle))*180/pi
std(abs(phase-myangle))*180/pi
quantile(abs(phase-myangle)*180/pi)
quantile(abs(phase-myangle)*180/pi,[0:0.1:1])'

myphase=(1-cos(myangle))/2;


max(abs(actualphase-myphase))
mean(abs(actualphase-myphase))
median(abs(actualphase-myphase))
sqrt(meansq(actualphase-myphase))
std(abs(actualphase-myphase))
quantile(abs(actualphase-myphase))
quantile(abs(actualphase-myphase),[0:0.1:1])'

