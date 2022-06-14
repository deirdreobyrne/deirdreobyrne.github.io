format long;
ang=csvread("i.csv");
t = ang(:,1);
phaseangle=ang(:,2);
mp = @(T) polyval([-6.797e-8,14.348e-6,89.97e-4,477198.86763133,134.96341138]/180*pi,T/36525)
D = @(T) polyval([-0.884e-8,1.832e-6,-16.3e-4,445267.11151675,297.8502042]/180*pi,T/36525)
F = @(T) polyval ([0.116e-8,-0.284e-6,-34.029e-4,483202.01752731,93.27209932]/180*pi,T/36525)
M = @(T) polyval ([0.041e-6,-1.536e-4,35999.05029094,357.52910918]/180*pi,T/36525)
dt=t(2)-t(1)
n = rows(t)
freq=linspace(0,1/dt,n);
degjcy = freq*360*36525;

phase = (1 - cos(phaseangle))/2;

myphaseangle = mod(D(t), 2*pi);
myphase = (1 - cos(myphaseangle))/2;

sumsq(myphase - phase) ans = 543.1683364341201
max(abs(myphase-phase)) # ans = 8.601478848068173e-02

fitfunc = @(p) sumsq((1-cos(myphaseangle+p(1)*sin(mp(t))))/2 - phase);
testp = [0.121483398437500];
[arg(1),fval]=fminsearch(fitfunc,testp) # fval 62.51789991186224
myphaseangle = myphaseangle + arg(1)*sin(mp(t));

myphase = (1 - cos(myphaseangle))/2;
dphase = myphase - phase;
max(abs(dphase)) # ans = 3.296245371903672e-02
dfft = fft(dphase);
[val,index]=max(abs(dfft)) # val = 2245.248559630474
degjcy(index) # ans = 409301.928142410
# plot(degjcy,abs(dfft));

fitfunc = @(p) sumsq((1-cos(myphaseangle+p(1)*sin(M(t))))/2 - phase);
testp=[-3.652992968750000e-02];
[arg(2),fval]=fminsearch(fitfunc,testp) # fval 19.07154369275808
myphaseangle = myphaseangle + arg(2)*sin(M(t));

myphase = (1 - cos(myphaseangle))/2;
dphase = myphase - phase;
max(abs(dphase)) # ans = 1.748781084608919e-02
dfft = fft(dphase);
[val,index]=max(abs(dfft))
degjcy(index) # ans = 858530.8580987131
# plot(degjcy,abs(dfft));

fitfunc = @(p) sumsq((1-cos(myphaseangle+p(1)*sin(2*D(t)-mp(t))))/2 - phase);
testp=[0.001655]; # actual fit is 1.660861328125000e-02 -- I think a bit better is out there, though
[arg(3),fval]=fminsearch(fitfunc,testp) # fval = 10.08844159752731
myphaseangle = myphaseangle + arg(3)*sin(2*D(t)-mp(t));

myphase = (1 - cos(myphaseangle))/2;
dphase = myphase - phase;
max(abs(dphase)) # ans = 1.059048426308473e-02
dfft = fft(dphase);
[val,index]=max(abs(dfft)) # val = 819.3529219301031
degjcy(index) # ans = 445216.1659549064
# plot(degjcy,abs(dfft));

fitfunc = @(p) sumsq((1-cos(myphaseangle+p(1)*sin(2*D(t))))/2 - phase);
testp=[1.216943359375000e-02];
[arg(4),fval]=fminsearch(fitfunc,testp) # fval = 5.188965085172714
myphaseangle = myphaseangle + arg(4)*sin(2*D(t));

myphase = (1 - cos(myphaseangle))/2;
dphase = myphase - phase;
max(abs(dphase)) # ans = 8.799220108316130e-03
dfft = fft(dphase);
[val,index]=max(abs(dfft)) # val = 713.7318493715557
degjcy(index) # ans = 922534.4439209822
# plot(degjcy,abs(dfft));

fitfunc = @(p) sumsq((1-cos(myphaseangle+p(1)*sin(2*mp(t))))/2 - phase);
testp=[4.445800781250000e-03];
[arg(5),fval]=fminsearch(fitfunc,testp) # fval = 4.530590347658622
myphaseangle = myphaseangle + arg(5)*sin(2*mp(t));

myphase = (1 - cos(myphaseangle))/2;
dphase = myphase - phase;
max(abs(dphase)) # ans = 7.058738796966169e-03
dfft = fft(dphase);
[val,index]=max(abs(dfft)) # val = 701.1142505396080
degjcy(index) # ans = 922534.4439209822
# plot(degjcy,abs(dfft));

fitfunc = @(p) sumsq((1-cos(myphaseangle+p(1)*sin(D(t))))/2 - phase);
testp=[1.938964843750000e-03];
[arg(6),fval]=fminsearch(fitfunc,testp) # fval = 4.341681598598860
myphaseangle = myphaseangle + arg(6)*sin(D(t));

myphase = (1 - cos(myphaseangle))/2;
dphase = myphase - phase;
max(abs(dphase)) # ans = 6.514760489687688e-03
dfft = fft(dphase);
[val,index]=max(abs(dfft)) # val = 701.3615254759620
degjcy(index) # ans = 922534.4439209822 <-- D + mp -- makes no difference!
# plot(degjcy,abs(dfft));

# Tried
# D + mp
# F + mp


fitfunc = @(p) sumsq((1-cos(myphaseangle+p(1)*sin(2*(F(t)-D(t)))))/2 - phase);
testp=[1.683105468750000e-03];
[arg(7),fval]=fminsearch(fitfunc,testp) # fval = 4.246048869823010
myphaseangle = myphaseangle + arg(7)*sin(2*(F(t)-D(t)));

myphase = (1 - cos(myphaseangle))/2;
dphase = myphase - phase;
max(abs(dphase)) # ans = 6.199360409463117e-03
dfft = fft(dphase);
[val,index]=max(abs(dfft))
degjcy(index) # ans = 922534.4439209822 <-- D + mp -- makes no difference!
# plot(degjcy,abs(dfft));

dphase = abs(dphase);
mean(dphase) # ans = 2.422249336288884e-03
median(dphase) # ans = 2.351272490478062e-03
sqrt(meansq(dphase)) # ans = 2.845820708813536e-03
std(dphase) # ans = 1.493789782945632e-03
quantile(dphase)'
quantile(dphase,[0:0.1:1])'


# Quantiles
#
#   1.534401194991375e-08
#   1.120117296322859e-03
#   2.351272490478062e-03
#   3.644469004007000e-03
#   6.199360409463117e-03
#
# and
#
#   1.534401194991375e-08
#   4.315384808369127e-04
#   8.840477452019269e-04
#   1.360622956519747e-03
#   1.841673525083976e-03
#   2.351272490478062e-03
#   2.890898819356420e-03
#   3.392637984892095e-03
#   3.903925012936816e-03
#   4.468324609582398e-03
#   6.199360409463117e-03
#
