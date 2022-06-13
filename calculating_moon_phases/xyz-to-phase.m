format long;

printf("\nLoading data for phase angle calculation\n");
sun=csvread("SUN.csv");
moon=csvread("MOON.csv");

ang=[zeros(524288,2)];

printf("Starting...\n");

for i = 1 : 524288
ang(i,1)=sun(i,1)-2451545;
s=sun(i,[2:4]);
m=-moon(i,[2:4]);
s=s+m;
x = dot(s,m);
n = cross(s,m);
n = n / norm(n);
if (n(3)<0)
n = -n;
endif;
y = det([n; s; m]);
ang(i,2)=mod(atan2(y,x),2*pi);
if (mod(i,5000) == 0)
  printf("(%d / 524288)\n",i);
endif;
end;

printf("\nDONE!\n");

csvwrite("i1.csv",ang);



