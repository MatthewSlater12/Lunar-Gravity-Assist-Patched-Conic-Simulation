function X = COE2RV(coe, mu)
%establising variables from classic elements
a = coe(1);
e = coe(2);
inc = coe(3);
RAAN = coe(4);
argp = coe(5);
Ta = coe(6);

%Calculating JP matrix values
a1 = (cos(RAAN)*cos(argp)) - (sin(RAAN)*cos(inc)*sin(RAAN));
b1 = (-cos(RAAN)*sin(argp)) - (sin(RAAN)*cos(inc)*cos(argp));
c1 = sin(argp)*sin(inc);

a2 = (sin(RAAN)*cos(argp)) + (cos(RAAN)*cos(inc)*sin(argp));
b2 = (-sin(RAAN)*sin(argp)) + (cos(RAAN)*cos(inc)*cos(argp));
c2 = -cos(RAAN)*sin(inc);

a3 = sin(inc)*sin(argp);
b3 = sin(inc)*cos(argp);
c3 = cos(inc);

%Calculating r and h, Inisalising JP matrix
r = (a*(1-e^2))/(1+e*cos(Ta));
h = sqrt(mu*a*(1-e^2));
jp = [a1, b1, c1; a2, b2, c2; a3, b3, c3];

%Inisalising position and velocity translation matrix
postran = [r*cos(Ta); r*sin(Ta); 0];
veltran = [-(mu/h)*sin(Ta); (mu/h)*(cos(Ta)+e); 0];

%finding position and velocity vectors
posvec = jp*postran;
velvec = jp*veltran;

%Reshapeing Vectors into a linier array
posvec = reshape(posvec, [1,3]);
velvec = reshape(velvec, [1,3]);

X = [posvec velvec];
