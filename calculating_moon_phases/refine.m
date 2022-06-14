format long;

ang=csvread("i.csv");
t = (ang([1:524288],1)-ang(1,1))/36525;
phaseangle=ang([1:524288],2);

phase = (1 - cos(phaseangle))/2;

args = [
   4.455894194873053e+00 ... # mp @ 1970
   deg2rad(477198.86763133) ...
   4.847069099776490e+00 ... # D @ 1970
   deg2rad(445267.11151675) ...
   3.711731242352926e+00 ... # F @ 1970
   deg2rad(483202.01752731) ...
   6.245032550629560e+00 ... # M @ 1970
   deg2rad(35999.05029094) ...
   1.214833984375000e-01 ... # mp
  -3.652992968750000e-02 ... # M
   1.660861328125000e-02 ... # 2D-mp
   1.216943359375000e-02 ... # 2D
   4.445800781250000e-03 ... # 2mp
   1.938964843750000e-03 ... # D
   1.683105468750000e-03 ... # 2(F-D)
  ];

fitfunc = @(p) sumsq((1-cos( ...
    mod(polyval([p(4),p(3)],t),2*pi) ...
  + p(9)*sin(polyval([p(2),p(1)],t)) ... # mp
  + p(10)*sin(polyval([p(8),p(7)],t)) ... # M
  + p(11)*sin(2*polyval([p(4),p(3)],t)-polyval([p(2),p(1)],t)) ... # 2D-mp
  + p(12)*sin(2*polyval([p(4),p(3)],t)) ... # 2D
  + p(13)*sin(2*polyval([p(2),p(1)],t)) ... # 2mp
  + p(14)*sin(polyval([p(4),p(3)],t)) ... # D
  + p(15)*sin(2*(polyval([p(6),p(5)],t)-polyval([p(4),p(3)],t))) ... # 2(F-D)
  ))/2 - phase);

fitfunc(args)

[result,fval]=fminsearch(fitfunc,args);
result'
fval

myphase = (1-cos( ...
    mod(polyval([result(4),result(3)],t),2*pi) ...
  + result(9)*sin(polyval([result(2),result(1)],t)) ... # mp
  + result(10)*sin(polyval([result(8),result(7)],t)) ... # M
  + result(11)*sin(2*polyval([result(4),result(3)],t)-polyval([result(2),result(1)],t)) ... # 2D-mp
  + result(12)*sin(2*polyval([result(4),result(3)],t)) ... # 2D
  + result(13)*sin(2*polyval([result(2),result(1)],t)) ... # 2mp
  + result(14)*sin(polyval([result(4),result(3)],t)) ... # D
  + result(15)*sin(2*(polyval([result(6),result(5)],t)-polyval([result(4),result(3)],t))) ... # 2(F-D)
  ))/2;

dphase = myphase - phase;
dfft = fft(dphase);

dphase = abs(dphase);
max(dphase) # ans = 2.812761265502894e-03
mean(dphase) # ans = 6.428006984862622e-04
median(dphase) # ans = 5.042315038243095e-04
sqrt(meansq(dphase)) ans = 8.277289335442166e-04
std(dphase) # ans = 5.214815127517401e-04
quantile(dphase)'
quantile(dphase,[0:0.1:1])'

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
# 1,829,700
#






