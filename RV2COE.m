function X = RV2COE(r,v,mu)

r2d = 180/pi;

MagR = norm(r);
MagV = norm(v);

E = 0.5*(MagV)^2-(mu/MagR);
a = -(mu/(2*E));

h = cross(r,v);
h = h.';

e = ((1/mu)*(cross(v,h))-(r/MagR));
Mage = norm(e);

MagH = norm(h);

Khat = [0; 0; 1];
Jhat = [0; 1; 0];
Ihat = [1; 0; 0];

i = acos(dot(h,Khat)/MagH);

Nhat = (cross(Khat,h)/(norm(cross(Khat,h))));

RAAN = atan2(dot(Nhat,Jhat),dot(Nhat,Ihat));

NhaT = (cross(h,Nhat)/norm(cross(h,Nhat)));
ArgP = atan2(dot(e,NhaT),dot(e,Nhat));

iehat = e/Mage;
iphat = cross(h,e)/norm(cross(h,e));

Ta = atan2(dot(r,iphat),dot(r,iehat));

i = i*r2d;
RAAN = RAAN*r2d;
ArgP = ArgP*r2d;
Ta = Ta*r2d;

X = [a, Mage, i , RAAN, ArgP, Ta];