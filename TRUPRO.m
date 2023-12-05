function Tati = TRUPRO(a,mu,e)

%Initalising Array
Mti = zeros(1,5000);
Eti = zeros(1,5000);

%Finding orbital period and establishing time vector 
P = 2*pi/mu^(1/2)*a^(3/2);
t = linspace(0,P,5000);

%Establising M0 and Setting the intial value in the M hold array
M0 = 0;
Mti(1) = M0;

%Calculating mean motion
N = (mu/a^3)^(1/2);

%Iterative index loop finding the values of M over orbit time 
i = 2;
while i ~= 5001
    Mti(i) = M0 + N*t(i);
    i = i + 1;
end

%Iterative index loop calculating the value of E over time using the M over
%time array
ie = 1;
while ie ~= 5001
    Eti(ie) = Kepler(e, Mti(ie), 10^-10);
    ie = ie + 1;
end

%Finding the value of T over time
Tati = 2*atan2(sqrt(1+e)*tan(Eti/2), sqrt(1-e));