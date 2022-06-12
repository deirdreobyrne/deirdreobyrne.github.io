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
F = @(T) polyval ([0.116e-8,-0.284e-6,-34.029e-4,483202.01752731,93.27209932]/180*pi,T/36525);
M = @(T) polyval ([0.041e-6,-1.536e-4,35999.05029094,357.52910918]/180*pi,T/36525);
dt=t(2)-t(1)
n = rows(t)
freq=linspace(0,1/dt,n);

fnc1 = @(p) sumsq(phase2-p(1)*sin(mp(t)));
p1=[0.109764];
[xmin1,fval1]=fminsearch(fnc1,p1)
# xmin1 = 0.109702964843750
# fval1 = 601.7590645893829
phase2 = phase2 - xmin1 * sin(mp(t));

fnc1 = @(p) sumsq(phase2-p(1)*sin(M(t)));
p1=[-0.036652];
[xmin2,fval2]=fminsearch(fnc1,p1)
# xmin2 = -3.659096484375000e-02
# fval2 = 250.6708194814297
phase2 = phase2 - xmin2 * sin(M(t));

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*D(t)-mp(t)));
p1=[0.022235 ];
[xmin3,fval3]=fminsearch(fnc1,p1)
# xmin3 = 2.183827148437500e-02
#fval3 = 125.6149640393175
phase2 = phase2 - xmin3*sin(2*D(t)-mp(t));

### csvwrite("phase2-pre-D.csv",phase2);

#
# 2/pi*abs(mod(x-pi/2,2*pi)-pi)-1 gives me a triangle wave from 0 -- 1 -- 0 -- -1 -- 0
# pi-mod(x,2*pi) gives me a triangle wave which starts at pi, goes to -pi, and immediately repeats
#

#
# CONCLUSION - I got rid of all terms except those involving 2n*D up to 12D, and what I got was an utter mess.
# I have to go through all the 2n*D terms!!!
#





## Looks like there is a left sawtooth component in D!
## Left sawtooth p(2)*(pi-mod(D(t),2*pi)) gave xmin4 = 1.358536317940849e-02   1.489650600516604e-03 fval4 = 62.73527662318310
# Now try a more probable triangle wave p(2)*(2/pi*abs(mod(D(t)-pi/2,2*pi)-pi)-1)
## Gives xmin4 = 1.508457477088387e-02   2.375603692198599e-03 fval4 = 65.00526145508256
# xmin4(2) is pretty close to the ratio of the moon's orbit semi-major axis to the earth's!!!
#   -- BUT I think my triangle goes from -0.5 to 0.5, so it would be half of the ratio
# BUT it has increased the magnitude of the harmonic terms.
fnc1 = @(p) sumsq(phase2-p(1)*sin(2*D(t)));
p1=[0.011484 ];
[xmin4,fval4]=fminsearch(fnc1,p1)
# xmin4 = 1.508507421875000e-02
# fval4 = 65.99357998471909
phase2 = phase2 - xmin4 * sin(2*D(t));

pfft=fft(phase2);
[maxv,maxi]=max(abs(pfft))
#maxv = 922.1859014132265
#maxi = 4758
freq(maxi)*360*36525
# ans = 954435.9177320818

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*mp(t)));
p1=[0.003735   ];
[xmin5,fval5]=fminsearch(fnc1,p1)
#xmin5 = 3.704482421875000e-03
#fval5 = 62.37158817387459
phase2 = phase2 - xmin5 * sin(2*mp(t));

pfft=fft(phase2);
[maxv,maxi]=max(abs(pfft))
# maxv = 881.5598922056477
# maxi = 8878
freq(maxi)*360*36525
# 1781065.302019696

fnc1 = @(p) sumsq(phase2-p(1)*sin(4*D(t)));
p1=[0.003  ];
[xmin6,fval6]=fminsearch(fnc1,p1)
# xmin6 = 3.244140625000000e-03
# fval6 = 59.62133016885514
phase2 = phase2 - xmin6 * sin(4*D(t));

pfft=fft(phase2);
[maxv,maxi]=max(abs(pfft))
#maxv = 687.4306761403453
#maxi = 17755
freq(maxi)*360*36525
# 3562130.604039391 = 2*1781065.302019696 = 4*890532.6510098478 = 6*593688.4340065651 = 8*445266.3255049239
fnc1 = @(p) sumsq(phase2-p(1)*sin(8*D(t)));
p1=[0.002  ];
[xmin7,fval7]=fminsearch(fnc1,p1)
#xmin7 = 2.518798828125000e-03
#fval7 = 57.97312530292286
phase2 = phase2 - xmin7 * sin(8*D(t));
pfft=fft(phase2);
[maxv,maxi]=max(abs(pfft))
# maxv = 543.6182275472694
# maxi = 13317

freq(maxi)*360*36525
# 6D
fnc1 = @(p) sumsq(phase2-p(1)*sin(6*D(t)));
p1=[0.001  ];
[xmin8,fval8]=fminsearch(fnc1,p1)
# xmin8 = 2.342773437500000e-03
# fval8 = 54.84151215213033
phase2 = phase2 - xmin8 * sin(6*D(t));

pfft=fft(phase2);
[maxv,maxi]=max(abs(pfft))
# ౳maxv = 447.0957005046549
# maxi = 379
freq(maxi)*360*36525
# ans = 75841.23962638783
# 2D-2F?

fnc1 = @(p) sumsq(phase2-p(1)*sin(2*(D(t)-F(t))));
p1=[-0.0015  ];
[xmin9,fval9]=fminsearch(fnc1,p1)
#xmin9 = -1.744140625000000e-03
#fval9 = 54.04997806701135
phase2 = phase2 - xmin9 * sin(2*(D(t)-F(t)));

pfft=fft(phase2);
[maxv,maxi]=max(abs(pfft))
#maxv = 444.2175210652240
#maxi = 26632
freq(maxi)*360*36525
#ans = 5343195.906059086 12D!

##
## STOP HERE!!!
###

fnc1 = @(p) sumsq(phase2-p(1)*sin(D(t)));
p1=[0.00192   ];
[xmin10,fval10]=fminsearch(fnc1,p1)
#xmin9 = 1.965776367187500e-03
#fval9 = 54.04997806701135 n/a
phase2 = phase2 - xmin10 * sin(D(t));

pfft=fft(phase2);
[maxv,maxi]=max(abs(pfft))


actualphase = (1-cos(phase))/2;
myphase=(1-cos(polyval(poly4,t) + xmin1 * sin(mp(t)) + xmin2 * sin(M(t)) + xmin3*sin(2*D(t)-mp(t)) + xmin4 * sin(2*D(t)) ...
  + xmin5 * sin(2*mp(t)) + xmin6 * sin(4*D(t))+ xmin7 * sin(8*D(t))+ xmin8 * sin(6*D(t))))/2;

max(abs(actualphase-myphase))
# 4.183684770429741e-03 -- 0.0042
mean(abs(actualphase-myphase))
# ans = 9.231399169979842e-04 -- 0.0009
median(abs(actualphase-myphase))
# ans = 7.966864956935327e-04
sqrt(meansq(abs(actualphase-myphase)))
# ans = 1.147918894841079e-03
std(abs(actualphase-myphase))
# ans = 6.822985935121816e-04
quantile(abs(actualphase-myphase),[0:0.1:1])'
#ans =
#
#   4.934626263164432e-09
#   1.427066555683476e-04
#   2.925694626539521e-04
#   4.501522089345623e-04
#   6.170882512702925e-04
#   7.966864956935327e-04
#   9.941175769826562e-04
#   1.218003832573989e-03
#   1.488449995713265e-03
#   1.863680528871447e-03
#   4.183684770429741e-03
quantile(abs(actualphase-myphase))
#ans =
#
#   4.934626263164432e-09   3.713515682289281e-04   7.966864956935327e-04   1.345582460542727e-03   4.183684770429741e-03

myphase=(1-cos(polyval(poly4,t) + xmin1 * sin(mp(t)) + xmin2 * sin(M(t)) + xmin3*sin(2*D(t)-mp(t)) + xmin4 * sin(2*D(t)) ...
  + xmin5 * sin(2*mp(t))))/2;
max(abs(actualphase-myphase))
# 4.400805384995765e-03
mean(abs(actualphase-myphase))
# ans = 9.837375116916378e-04
median(abs(actualphase-myphase))
# ans = 8.513377903127117e-04
sqrt(meansq(abs(actualphase-myphase)))
# ans = 1.226054633530540e-03
std(abs(actualphase-myphase))
# ans = 7.317591774780291e-04
quantile(abs(actualphase-myphase),[0:0.1:1])'
#   2.578901536764988e-09
#   1.413141016769481e-04
#   2.991023606531030e-04
#   4.689526137639177e-04
#   6.514166687180857e-04
#   8.513377903127117e-04
#   1.067065126844591e-03
#   1.310459170393508e-03
#   1.607521498874615e-03
#   1.980552194768048e-03
#   4.400805384995765e-03


