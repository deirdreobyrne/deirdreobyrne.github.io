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
These values will be much larger than _2 &#03c0; _ for modern times but, since most math libraries will accept such large values, we can proceed.

Meeus' formula 48.4 calculates the quantity _i_ - the selenocentric elongation of the Earth from the Sun - which is a quantity whose values _decreases_ with the
passage of time. Instead we are going to calculate a quantity _E_ - the geocentric elongation of the Moon from the Sun - where _E = pi - i_.
