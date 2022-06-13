format long;

printf("\nLoading data for phase angle calculation\n");
sun=csvread("SUN.csv");
moon=csvread("MOON.csv");

ang=[zeros(524288,2)];
#psi=[zeros(1,524288)];

printf("Starting...\n");

for i = 1 : 524288
sd = deg2rad(sun(i,3));
sr=deg2rad(sun(i,2));
md=deg2rad(moon(i,3));
mr=deg2rad(moon(i,2));
#psi(i)=acos(sin(sd)*sin(md)+cos(sd)*cos(md)*cos(sr-mr));
psi=acos(sin(sd)*sin(md)+cos(sd)*cos(md)*cos(sr-mr));
r = sun(i,4)*149.59787066;
rm = moon(i,4);
ang(i,1)=sun(i,1)-2451545;
ang(i,2)=psi + atan2(rm*sin(psi), r-rm*cos(psi));
if (mod(i,5000) == 0)
  printf("(%d / 524288)\n",i);
endif;
end;

printf("\nDONE!\n");

csvwrite("i1.csv",ang);



