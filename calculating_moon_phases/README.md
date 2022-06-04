# Calculating moon phases

## Introduction

There seems to be a proliferation of a poor algorithm for calculating the moon's phase on the Internet. This algorithm, for instance, assumes the moon's phase
progression is smooth, whereas in fact the progression of the phase is slow near to new and full moon, and quicker nearer to the quarters.

## A better algorithm

My source for this algorithm is Jean Meeus' excellent 1991 book "Astronomical Algorithms" published by Willmann-Bell of Richmond, Virginia, USA. Specifically
the algorithm is based on formula 48.4 on page 346.

## Derivation of the algorithm

The time argument for Meeus' formulae is _T_, the number of Julian centuries (of 36,525 days each) since 2000 Jan 1 12:00:00 UT. Most computer systems have some
sort of _millis()_ function, which is the number of milliseconds since 1970 Jan 1 00:00:00. Hence 
```
T = (millis() - 946728000000) / 3155760000000
```
Meeus works in degrees, whereas computers work in radians. So we convert Meeus' formulae 47.2 - 47.4 into radians and, since our main formula 48.4 is of relatively low
accuracy, we drop unnecessary terms
```
D = (7771.37714483372 * T) + 5.19846652984
M = (628.30195516723 * T) + 6.24006012726
MP = (8328.69142475915 * T) + 2.35555563685
```
These values will be much larger than _2pi_ for modern times but, since most math libraries will accept such large values, we can proceed.

Meeus' formula 48.4 calculates the quantity _i_ - the selenocentric elongation of the Earth from the Sun - which is a quantity whose values _decreases_ with the
passage of time. Instead we are going to calculate a quantity _E_ - the geocentric elongation of the Moon from the Sun - where _E = pi - i_.
```
E = D + (0.109764 * sin (MP)) - (0.036652 * sin(M)) + (0.022235 * sin(2*D - MP)) + (0.011484 * sin(2*D)) + (0.003735 * sin(2*MP)) + (0.00192 * sin(D))
```
From this, we can calculate the fraction of the moon's diameter which is illuminated
```
FRAC = (1 - cos(E)) / 2
```
We note that
```
if (sin(E) >= 0) {
  # Moon is waxing
} else {
  # Moon is waning
}
```

## A worked example

Our timestamp is 2022 June 4 09:31:10 UT

The value of millis() is 1654335070000, hence
```
T = 0.224227149720
D = 1747.752213114103
M = 147.122416697834
MP = 1869.874294706978
E =  1747.7522131 - 0.0643708 - 0.0186059 - 0.0220047 + 0.0101724 + 0.0035484 + 0.0016427 = 1747.7522131 - 0.0896179 = 1747.6625952
FRAC = 0.2039
```
We note that _sin(E)_ is positive, hence the moon is waxing.

<!--

I'm not sure if Stellarium calculates the geocentric phase angle or the topocentric one.


### Comparison to Stellarium

The Stellarium program gives a phase angle of _125d 25' 52.4" = 125.43122d_ for this instant, and a value of _FRAC_ of _21.0%_.
Our value of _E_ gives a phase angle of _126.3093d_.
-->

<!-- 

The better performance of this formula is probably a coincidence!


### Comparison to simpler formula

I have a simpler formula in my (current) moon phase bangle widget -

```
k = (millis - 947168438000) / 2551442878;
mp = 3.5179607 + 6.73377583059 * k;
m =  0.0445652 + 0.50798430438 * k;
d =  6.28318530718 * k;
t = d + 0.109764 * Math.sin(mp) - 0.036652 * Math.sin(m);
k = (1 - Math.cos(t))/2;
```
This actually gives a slightly better result!!!
```
k = 277.16341921568976626722661827117
mp = 1869.8742941382957222058568853447 OK
m = 140.8392319098644910991419028743 = 147.12241721704407757606718964086 OK
d = 1741.4691233037928187393266814888 = 1747.7523086109724052162519682554 A bit off, but OK

t = d - 0.06437076 - 0.01860583 = 1.02679322 - 0.06437076 - 0.01860583 = 0.94381663 (54.0767 deg -- 125.9233)
k = 0.20665
```
-->

## Rate of change of the phase

It is sensible to cache the calculations for the phase. Hence we need to know how quickly the phase changes, so we can get an idea as to when the cache should
be considered stale.

Using calculus, we can calculate that the maximum rate of change of the phase is a bit slower than 1 part in 800,000 per second.
