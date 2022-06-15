pkg load optim;
format long;

#ang=csvread("i.csv");
t = (ang([1:524288],1)-ang(1,1))*86400;
dt = 45.0 + 50.0*t/86400/36525;
t = t + dt / 86400;
phaseangle=ang([1:524288],2);
dt=t(2)-t(1);
n = rows(t);
freq=linspace(0,1/dt,n);
degjcy = freq*360/86400/36525;

phase = (1 - cos(phaseangle))/2;

args = [
   4.455933781445361e+00 ...
   3.789022554950839e+05 ...
   4.847300486062351e+00 ...
   4.060747530841359e+05 ...
   3.711741520539088e+00 ...
   3.741949504350046e+05 ...
   6.245093325555559e+00 ...
   5.022681132791225e+06 ...
   1.098294220753609e-01 ...
  -3.647737653400068e-02 ...
   2.181414668659628e-02 ...
   1.343381074086885e-02 ...
   3.665903202450143e-03 ...
   1.958316814329805e-03 ...
   1.106039239677234e-03 ...
   1.066672928781793e-03 ...
  -1.240209781239226e-03 ...
   9.823080994024243e-04 ...
   1.168396449167793e-03 ...
   8.193661471503118e-04 ...
  -6.220874755401757e-04 ...
  -5.298160075353295e-04 ...
  -1.004011787385146e-04 ...
   2.736281167893605e-04 ...
  -6.578343139873168e-05 ...
   2.068112392990309e-04 ...
   1.890063960917113e-04 ...
   2.455817964409304e-06 ...
  -4.495544450489619e-06 ...
   8.222196325073695e-07 ...
   3.190588938867953e-05 ...
   ];

dpargs = [
  1e-3, 1e-6, 1e-3, 1e-6, 1e-3, 1e-6, 1e-3, 1e-5, ...
  1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, ...
  1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, ...
  1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1 ...
  1e-1 ...
  ];

#fitfunc = @(p) sumsq((1-cos( ...
func = @(t,p) (1-cos( ...
    mod(polyval([1/p(4),p(3)],t),2*pi) ...
  + p(9)*sin(polyval([1/p(2),p(1)],t)) ... # mp
  + p(10)*sin(polyval([1/p(8),p(7)],t)) ... # M
  + p(11)*sin(2*polyval([1/p(4),p(3)],t)-polyval([1/p(2),p(1)],t)) ... # 2D-mp
  + p(12)*sin(2*polyval([1/p(4),p(3)],t)) ... # 2D
  + p(13)*sin(2*polyval([1/p(2),p(1)],t)) ... # 2mp
  + p(14)*sin(polyval([1/p(4),p(3)],t)) ... # D
  + p(15)*sin(2*(polyval([1/p(6),p(5)],t)-polyval([1/p(4),p(3)],t))) ... # 2(F-D)
  + p(16)*sin(2*(polyval([1/p(4),p(3)],t)-polyval([1/p(2),p(1)],t))) ... # 2(D-mp)
  + p(17)*sin(2*polyval([1/p(6),p(5)],t)) ... # 2F
  + p(18)*sin(2*polyval([1/p(4),p(3)],t)-polyval([1/p(8),p(7)],t)-polyval([1/p(2),p(1)],t)) ... # 2D - M - mp
  + p(19)*sin(2*polyval([1/p(4),p(3)],t)+polyval([1/p(2),p(1)],t)) ... # 2D + mp
  + p(20)*sin(2*polyval([1/p(4),p(3)],t)-polyval([1/p(8),p(7)],t)) ... # 2D - M
  + p(21)*sin(polyval([1/p(8),p(7)],t)-polyval([1/p(2),p(1)],t)) ... # M - mp
  + p(22)*sin(polyval([1/p(8),p(7)],t)+polyval([1/p(2),p(1)],t)) ... # M + mp
  + p(23)*sin(polyval([1/p(2),p(1)],t)+2*polyval([1/p(6),p(5)],t)) ... # mp + 2F
  + p(24)*sin(polyval([1/p(2),p(1)],t)-2*polyval([1/p(6),p(5)],t)) ... # mp - 2F
  + p(25)*sin(4*polyval([1/p(4),p(3)],t)-polyval([1/p(2),p(1)],t)) ... # 4D - mp
  + p(26)*sin(3*polyval([1/p(2),p(1)],t)) ... # 3mp
  + p(27)*sin(4*polyval([1/p(4),p(3)],t)-2*polyval([1/p(2),p(1)],t)) ... # 4D-2mp
  + p(28)*sin(polyval([1/p(6),p(5)],t)) ... # F
  + p(29)*sin(polyval([1/p(2),p(1)],t)+polyval([1/p(6),p(5)],t)) ... # mp+F
  + p(30)*sin(polyval([1/p(2),p(1)],t)-polyval([1/p(6),p(5)],t)) ... # mp-F
  + p(31)*sin(2*polyval([1/p(4),p(3)],t)-polyval([1/p(8),p(7)],t)+polyval([1/p(2),p(1)],t)) ... # 2D - M + mp
  ))/2;


#fitfunc(args)

[myphase,newargs,cvg,iter]=leasqr(t, phase, args, func, 0.00000001, 20, ones(size(phase)), dpargs);
cvg
iter
newargs

max(abs(myphase-phase))


#
#myphase = (1-cos( ...
#    mod(polyval([result(4),result(3)],t),2*pi) ...
#  + result(9)*sin(polyval([result(2),result(1)],t)) ... # mp
#  + result(10)*sin(polyval([result(8),result(7)],t)) ... # M
#  + result(11)*sin(2*polyval([result(4),result(3)],t)-polyval([result(2),result(1)],t)) ... # 2D-mp
#  + result(12)*sin(2*polyval([result(4),result(3)],t)) ... # 2D
#  + result(13)*sin(2*polyval([result(2),result(1)],t)) ... # 2mp
#  + result(14)*sin(polyval([result(4),result(3)],t)) ... # D
#  + result(15)*sin(2*(polyval([result(6),result(5)],t)-polyval([result(4),result(3)],t))) ... # 2(F-D)
#  ))/2;
#
#dphase = myphase - phase;
#dfft = fft(dphase);
#
#dphase = abs(dphase);
#max(dphase) # ans = 2.812761265502894e-03
#mean(dphase) # ans = 6.428006984862622e-04
#median(dphase) # ans = 5.042315038243095e-04
#sqrt(meansq(dphase)) # ans = 8.277289335442166e-04
#std(dphase) # ans = 5.214815127517401e-04
#quantile(dphase)'
#quantile(dphase,[0:0.1:1])'
#
# Quantiles
#
#   8.034322651617742e-10
#   2.310393715379833e-04
#   5.042315038243095e-04
#   9.337538482949437e-04
#   2.812761265502894e-03
#
#   8.034322651617742e-10
#   8.867288309312160e-05
#   1.827553005553539e-04
#   2.805617935229987e-04
#   3.857476484745770e-04
#   5.042315038243095e-04
#   6.437388032005665e-04
#   8.220775543042680e-04
#   1.066496447293541e-03
#   1.444741621761308e-03
#   2.812761265502894e-03


# My FFT gives that the next terms might be -
#
# 1,335,840  3D
# 369,375    3D - 2F
# 445216 - a residual of D??
# 521057     D - 2F??????
# 822616 -- 3D - M - mp
# 381417    3D - 2mp
# 32000     D - mp
# 1,829,700
#






