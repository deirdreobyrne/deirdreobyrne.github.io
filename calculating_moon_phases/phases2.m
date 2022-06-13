format long;
phases=csvread("phases.csv");
t = phases(:,1);
phase=phases(:,2);
poly4=polyfit(t,phase,4); ## Or I could use the value of D below?
poly4
# -2.522610894287007e-21   2.224416831388318e-16  -5.558948467500520e-12   2.127687346562486e-01   5.198337681802380e+00
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
#xmin1 = 0.109702964843750
#fval1 = 601.7590645893829
phase2 = phase2 - xmin1 * sin(mp(t));

fnc1 = @(p) sumsq(phase2-p(1)*sin(M(t)))
p1=[-0.036652];
[xmin2,fval2]=fminsearch(fnc1,p1)
#xmin2 = -3.659096484375000e-02
#fval2 = 250.6708194814297
phase2 = phase2 - xmin2 * sin(M(t));

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*D(t)-mp(t)))
p1=[0.022235 ];
[xmin3,fval3]=fminsearch(fnc1,p1)
#xmin3 = 2.183827148437500e-02
#fval3 = 125.6149640393175
phase2 = phase2 - xmin3*sin(2*D(t)-mp(t));

## RE-DONE!

fnc1 = @(p) sumsq(phase2-(ifelse(mod(D(t),pi) <= p(1)/2 | mod(D(t),pi) >= pi-p(1)/2, ...
  p(2)*mod(D(t)+p(1)/2,pi)/p(1), p(2)*(1-mod(D(t)-p(1)/2,pi)/(pi-p(1))))-p(2)/2)-p(3)*sin(2*D(t)))
p1=[0.22431, 0.03871, 0.011484 ];
[xmin4,fval4]=fminsearch(fnc1,p1)
# OLD -
#xmin4 = 1.508507421875000e-02
#fval4 = 65.99357998471909
# NEW -
#xmin4 = 1.619147787470724e-01   2.541848302958411e-02   6.589447562519009e-03
#fval4 = 56.68453435575019
phase2 = phase2 - (ifelse(mod(D(t),pi) <= xmin4(1)/2 | mod(D(t),pi) >= pi-xmin4(1)/2, ...
  xmin4(2)*mod(D(t)+xmin4(1)/2,pi)/xmin4(1), xmin4(2)*(1-mod(D(t)-xmin4(1)/2,pi)/ ...
  (pi-xmin4(1))))-xmin4(2)/2) - xmin4(3) * sin(2*D(t));

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*mp(t)))
p1=[0.003735   ];
[xmin5,fval5]=fminsearch(fnc1,p1)
#xmin5 = 3.704482421875000e-03
#fval5 = 53.06244608031772
phase2 = phase2 - xmin5 * sin(2*mp(t));

#fnc1 = @(p) sumsq(phase2-p(1)*sin(4*D(t)))
#p1=[0.003  ];
#[xmin6,fval6]=fminsearch(fnc1,p1)
##xmin6 = 3.244140625000000e-03
##fval6 = 59.62133016885514
#phase2 = phase2 - xmin6 * sin(4*D(t));
#
#fnc1 = @(p) sumsq(phase2-p(1)*sin(8*D(t)))
#p1=[0.002  ];
#[xmin7,fval7]=fminsearch(fnc1,p1)
##xmin7 = 2.518798828125000e-03
##fval7 = 57.97312530292286
#phase2 = phase2 - xmin7 * sin(8*D(t));
#
#fnc1 = @(p) sumsq(phase2-p(1)*sin(6*D(t)))
#p1=[0.001  ];
#[xmin8,fval8]=fminsearch(fnc1,p1)
##xmin8 = 2.831054687500000e-03
##fval8 = 55.84992477494759
#phase2 = phase2 - xmin8 * sin(6*D(t));
#
#fnc1 = @(p) sumsq(phase2-p(1)*sin(12*D(t)))
#p1=[0.001  ];
#[xmin9,fval9]=fminsearch(fnc1,p1)
##xmin9 = 1.823974609375000e-03
##fval9 = 54.97223532979434
#phase2 = phase2 - xmin9 * sin(12*D(t));

fnc1 = @(p) sumsq(phase2-p(1)*sin(D(t)))
p1=[0.001  ];
[xmin10,fval10]=fminsearch(fnc1,p1)
#xmin10 = 1.946044921875000e-03
#fval10 = 52.05837548024523
phase2 = phase2 - xmin10 * sin(D(t));

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*(F(t)-D(t))))
p1=[0.001  ];
[xmin11,fval11]=fminsearch(fnc1,p1)
#xmin11 = 1.732421875000000e-03
#fval11 = 51.26660815571190
phase2 = phase2 - xmin11 * sin(2*(F(t)-D(t)));

#fnc1 = @(p) sumsq(phase2-p(1)*sin(10*D(t)))
#p1=[0.001  ];
#[xmin12,fval12]=fminsearch(fnc1,p1)
##xmin12 = 2.159667968750000e-03
##fval12 = 51.94662778380171
#phase2 = phase2 - xmin12 * sin(10*D(t));

pfft=fft(phase2);
[maxv,maxi]=max(abs(pfft))
freq(maxi)*360*36525
# Going with 2F-6D

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*(F(t)-3*D(t))))
p1=[0.001  ];
[xmin12,fval12]=fminsearch(fnc1,p1)
#xmin12 = 1.579833984375000e-03
#fval12 = 50.62356039167265
phase2 = phase2 - xmin12 * sin(2*(F(t)-3*D(t)));

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*F(t)))
p1=[0.001  ];
[xmin13,fval13]=fminsearch(fnc1,p1)
#xmin13 = -1.716064453125000e-03
#fval13 = 49.86356679109232
phase2 = phase2 - xmin13 * sin(2*F(t));

pfft=fft(phase2);

actualphase = (1-cos(phase))/2;
myphase=(1-cos(polyval(poly4,t) + xmin1 * sin(mp(t)) + xmin2 * sin(M(t)) + xmin3*sin(2*D(t)-mp(t)) + xmin4(3) * sin(2*D(t)) ...
  + (ifelse(mod(D(t),pi) <= xmin4(1)/2 | mod(D(t),pi) >= pi-xmin4(1)/2,xmin4(2)*mod(D(t)+xmin4(1)/2,pi)/xmin4(1), ...
    xmin4(2)*(1-mod(D(t)-xmin4(1)/2,pi)/(pi-xmin4(1))))-xmin4(2)/2)
  + xmin5 * sin(2*mp(t)) + xmin10 * sin(D(t)) + xmin11*sin(2*(F(t)-D(t))) ...
  + xmin12 * sin(2*(F(t)-3*D(t))) + xmin13 * sin(2*F(t))   ))/2;

##
## DO AGAIN - RESULTS NOT AS GOOD?!
##

max(abs(actualphase-myphase))
#ans = 3.347400398570399e-03
mean(abs(actualphase-myphase))
#ans = 6.534936952047384e-04
median(abs(actualphase-myphase))
#ans = 5.514659669792477e-04
sqrt(meansq(abs(actualphase-myphase)))
#ans = 8.201897924756784e-04
std(abs(actualphase-myphase))
#ans = 4.956387339220904e-04
quantile(abs(actualphase-myphase))
# 2.038354984801316e-09   2.568847409428199e-04   5.514659669792477e-04   9.428177831294071e-04   3.347400398570399e-03
quantile(abs(actualphase-myphase),[0:0.1:1])'
#   2.038354984801316e-09
#   1.001372553499913e-04
#   2.037594189659209e-04
#   3.116419006228777e-04
#   4.261685020933381e-04
#   5.514659669792477e-04
#   6.895670473660465e-04
#   8.503985835515870e-04
#   1.051758413577020e-03
#   1.362115283503818e-03
#   3.347400398570399e-03

