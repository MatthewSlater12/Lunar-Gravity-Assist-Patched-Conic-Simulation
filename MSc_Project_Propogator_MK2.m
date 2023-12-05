clear; clc; close all

%Simulation Setup For each site set SiteSelect to either 1 for sutherland,
%2 for Cornwall and 3 for SaxaVord
SiteSelect = 3;

if SiteSelect == 1
    SIMV = [300, 300, 83, 77, 147, 3951, 1049];
    STitle = 'Sutherland';
elseif SiteSelect == 2
    SIMV = [230, 230, 70, 33, 50, 3852, 1150];
    STitle = 'Cornwall';
elseif SiteSelect == 3
    SIMV = [300, 300, 62, 3, 6, 3791, 1210];
    STitle = 'SaxaVord';
else
    error('Please select valid site selection value')
end

%%Orbit Calculation
%Establising Iterated Arrays
T_coe = zeros(6,5000);
T_X = zeros(6,5000);
I_coe = zeros(6,5000);
I_X = zeros(6,5000);
LU_coe = zeros(6,5000);
LU_X = zeros(6,5000);
coe_FO1 = zeros(6,5000);
X_FO1 = zeros(6,5000);
coe_FO2 = zeros(6,5000);
X_FO2 = zeros(6,5000);
coe_G = zeros(6,5000);
X_G = zeros(6,5000);
coe_G2 = zeros(6,5000);
X_G2 = zeros(6,5000);
%Degree to radian converstion
d2r = pi/180;
Re = 6378.137;

%Standard Gravititational Parameters
mu = 398600.4418;
muM = 4900;

%Initial Orbit Conditions
I_Ap = SIMV(1);
I_Pe = SIMV(2);
I_inc = SIMV(3)*d2r;
I_RAAN = 0*d2r;
I_argp = 0*d2r;
I_a = ((I_Ap+Re)+(I_Pe+Re))/2;
I_e = ((I_Ap+Re)/I_a)-1;
I_Tati = TRUPRO(I_a,mu,I_e);

%Transfer Orbit Conditions
T_Ap = 400000;
T_Pe = SIMV(2);
T_inc = SIMV(3)*d2r;
T_RAAN = 0*d2r;
T_argp = 0*d2r;
T_a = ((T_Ap+Re)+(T_Pe+Re))/2;
T_e = ((T_Ap+Re)/T_a)-1;
T_Tati = TRUPRO(T_a,mu,T_e);

%Convering the classical orbital eliments over time to cartisian
%coordinates array

for tt = 1:length(I_Tati)
    I_coe(:,tt) = [I_a,I_e,I_inc,I_RAAN,I_argp,I_Tati(tt)];
    I_X(:,tt) = COE2RV(I_coe(:,tt),mu);
end


for tt = 1:length(T_Tati)
    T_coe(:,tt) = [T_a,T_e,T_inc,T_RAAN,T_argp,T_Tati(tt)];
    T_X(:,tt) = COE2RV(T_coe(:,tt),mu);
end

%Adding Luner Orbit
%Circular orbit
LU_a = 400000+Re;
LU_e = 0;
LU_inc = 28.58*d2r;
LU_RAAN = 0;
LU_argP = 0;

LU_Ta = TRUPRO(LU_a,mu,LU_e);
for tt = 1:length(LU_Ta)
    LU_coe(:,tt) = [LU_a,LU_e,LU_inc,LU_RAAN,LU_argP,LU_Ta(tt)];
    LU_X(:,tt) = COE2RV(LU_coe(:,tt),mu);
end

%Closest Approch 

LU_Closest_Coe = [LU_a,LU_e,LU_inc,LU_RAAN,LU_argP,180*d2r];
SAT_Closest_Coe = [T_a,T_e,T_inc,T_RAAN,T_argp,180*d2r];

LU_Closest_PV = COE2RV(LU_Closest_Coe,mu);
SAT_Closest_PV = COE2RV(SAT_Closest_Coe,mu);


%% Vinf Map 

%Velocity Values

VgaVec = LU_Closest_PV(1,4:6);
VscVec = SAT_Closest_PV(1,4:6);
gaPos = LU_Closest_PV(1,1:3);
VinfVec = VscVec - VgaVec ;


Vga = norm(VgaVec);
Vsc = norm(VscVec);
Vinf = norm(VinfVec);

Renc = norm(SAT_Closest_PV(1,1:3));

vc = sqrt(mu/Renc);

GaFPA = 0;

P1 = (SAT_Closest_PV(1,1:3)/Renc);
P3 = (cross(P1,LU_Closest_PV(1,4:6))/norm(cross(P1,LU_Closest_PV(1,4:6))));
P2 = cross(P3,P1);

PI = [P1; P2; P3];

ScFPA = sqrt((T_a/Renc)*(1-T_e^2)/(2-(Renc/T_a)));
ScFPA = acosd(ScFPA);


VinfP = PI*VinfVec.';

CheckP = acosd(VinfP(2)/Vinf);
CheckC = atan2d(-VinfP(3),VinfP(1));

InBPump = acosd(((Vsc^2)-(Vinf^2)-(Vga^2))/(2*(Vinf*Vga)));
InBKrank = acosd((1/Vinf*sind(InBPump)*cosd(ScFPA))*((Vsc*sin(GaFPA))-Vga*sind(ScFPA)-Vinf*cosd(InBPump)*sind(GaFPA)));

%Orbit Size Equation

%Pumpcrank plus

PumpNew = linspace(0,180,2000);
CrankNew = linspace(0,360,2000);

[p,c] = meshgrid(PumpNew,CrankNew);

aNew = zeros(size(p));
eNew = zeros(size(p));
iNew = zeros(size(p));
ApNew = zeros(size(p));
PeNew = zeros(size(p));

for i = 1:length(CrankNew)
    for i2 = 1:length(PumpNew)
        VinfNewVec = (Vinf*sind(p(i,i2))*cosd(c(i,i2))*P1)+(Vinf*cosd(p(i,i2))*P2)-(Vinf*sind(p(i,i2))*sind(c(i,i2))*P3);
        VscEciNew = VinfNewVec+VgaVec;
        COENEW = RV2COE(gaPos,VscEciNew,mu);
        aNew(i,i2) = COENEW(1);
        eNew(i,i2) = COENEW(2);
        iNew(i,i2) = COENEW(3);
        ApNew(i,i2) = aNew(i,i2)*(1+eNew(i,i2));
        PeNew(i,i2) = aNew(i,i2)*(1-eNew(i,i2));
    end
end

[x1,y1] = contour(c,p,iNew,[0.2,0]);

[x2,y2] = contour(c,p,PeNew,[42164,42164]);

Int = InterX(x1,x2);

figure;
contour(c,p,iNew,[0.2,0],'b','LineWidth',0.5,'DisplayName','Orbits of 0.2 Inclination');
hold on 
contour(c,p,PeNew,[42164,42164],'r','LineWidth',1,'DisplayName','Orbits of Pe at GEO Ring');
scatter(InBKrank,InBPump,'DisplayName','Intercept Orbit');
scatter(Int(1,SIMV(4)),Int(2,SIMV(4)),'k','LineWidth',2,'DisplayName','Anti-Planet Encounter');
scatter(Int(1,SIMV(5)),Int(2,SIMV(5)),'g','LineWidth',2,'DisplayName','Planet Encounter');
title('Pump Crank Angle Graph For Lunar Intercept', STitle)
xlabel('Crank Angle (\circ)')
ylabel('Pump Angle (\circ)')
legend
hold off

FO_1_Pump = Int(2,SIMV(4));
FO_1_Crank = Int(1,SIMV(4));

FO_1_Vinf = (Vinf*sind(FO_1_Pump)*cosd(FO_1_Crank)*P1)+(Vinf*cosd(FO_1_Pump)*P2)-(Vinf*sind(FO_1_Pump)*sind(FO_1_Crank)*P3);

FO_1_Vsc = FO_1_Vinf + VgaVec;
FO_1_COE = RV2COE(gaPos,FO_1_Vsc,mu);

FO_1_Bending = acosd((dot(VinfVec,FO_1_Vinf))/(Vinf^2));

FO_1_Radius = ((muM*cscd(FO_1_Bending/2))/(Vinf^2))-(muM/(Vinf^2));


FO_2_Pump = Int(2,SIMV(5));
FO_2_Crank = Int(1,SIMV(5));

FO_2_Vinf = (Vinf*sind(FO_2_Pump)*cosd(FO_2_Crank)*P1)+(Vinf*cosd(FO_2_Pump)*P2)-(Vinf*sind(FO_2_Pump)*sind(FO_2_Crank)*P3);

FO_2_Vsc = FO_2_Vinf + VgaVec;
FO_2_COE = RV2COE(gaPos,FO_2_Vsc,mu);

FO_2_Bending = acosd((dot(VinfVec,FO_2_Vinf))/(Vinf^2));

FO_2_Radius = ((muM*cscd(FO_2_Bending/2))/(Vinf^2))-(muM/(Vinf^2));

FO_1_TaTi = TRUPRO(FO_1_COE(1,1),mu,FO_1_COE(1,2));

for tt = 1:length(FO_1_TaTi)
    coe_FO1(:,tt) = [FO_1_COE(1,1),FO_1_COE(1,2),FO_1_COE(1,3)*d2r,FO_1_COE(1,4)*d2r,FO_1_COE(1,5)*d2r,FO_1_TaTi(tt)];
    X_FO1(:,tt) = COE2RV(coe_FO1(:,tt),mu);
end

FO_2_TaTi = TRUPRO(FO_2_COE(1,1),mu,FO_2_COE(1,2));

for tt = 1:length(FO_2_TaTi)
    coe_FO2(:,tt) = [FO_2_COE(1,1),FO_2_COE(1,2),FO_2_COE(1,3)*d2r,FO_2_COE(1,4)*d2r,FO_2_COE(1,5)*d2r,FO_2_TaTi(tt)];
    X_FO2(:,tt) = COE2RV(coe_FO2(:,tt),mu);
end


%GeoOrbits
G_Ap = 42164;
G_Pe = 42164;
G_inc = 0*d2r;
G_RAAN = FO_1_COE(1,4)*d2r;
G_argp = FO_1_COE(1,5)*d2r;
G_a = ((G_Ap)+(G_Pe))/2;
G_e = ((G_Ap)/G_a)-1;
G_Tati = TRUPRO(G_a,mu,G_e);

for tt = 1:length(G_Tati)
    coe_G(:,tt) = [G_a,G_e,G_inc,G_RAAN,G_argp,G_Tati(tt)];
    X_G(:,tt) = COE2RV(coe_G(:,tt),mu);
end

SDV_1 = norm(T_X(4:6,1)-I_X(4:6,1));
SDV_2 = norm(X_FO1(4:6,1)-X_G(4:6,1));
SDV_T = SDV_1 + SDV_2;

G_Ap2 = 42164;
G_Pe2 = 42164;
G_inc2 = 0*d2r;
G_RAAN2 = FO_2_COE(1,4)*d2r;
G_argp2 = FO_2_COE(1,5)*d2r;
G_a2 = ((G_Ap2)+(G_Pe2))/2;
G_e2 = ((G_Ap2)/G_a2)-1;
G_Tati2 = TRUPRO(G_a2,mu,G_e2);

for tt = 1:length(G_Tati2)
    coe_G2(:,tt) = [G_a2,G_e2,G_inc2,G_RAAN2,G_argp2,G_Tati2(tt)];
    X_G2(:,tt) = COE2RV(coe_G2(:,tt),mu);
end

LDV_1 = norm(T_X(4:6,1)-I_X(4:6,1));
LDV_2 = norm(X_FO2(4:6,1)-X_G2(4:6,1));
LDV_T = LDV_1 + LDV_2;

%Orbit Time Equations

    T_OrbP = 2*pi*sqrt((T_a)^3/(mu));

    R_MMS = sqrt(mu/(FO_1_COE(1)^3));

    R_MML = sqrt(mu/(FO_2_COE(1)^3));
    
    MDT_S = abs(T2M(FO_1_COE(6),FO_1_COE(2)));

    MDT_L = (180*d2r)+((180*d2r)-T2M(FO_2_COE(6),FO_2_COE(2)));

    Transfer_TimeS = (T_OrbP/2) + ((MDT_S/R_MMS));
    Transfer_TimeL = (T_OrbP/2) + ((MDT_L/R_MML));


%% 3D plot 


[Xe,Ye,Ze] = sphere(50);

% Create figure
figure('units','normalized','outerposition',[0 0 1 1])

% Plot Earth
surf(Re*Xe, Re*Ye, Re*Ze, 'EdgeColor', 'none', 'FaceColor', 'c','DisplayName','Earth Surface'); hold on; axis equal;
xlabel('X, ECI (Km)','FontSize',20);
ylabel('Y, ECI (Km)','FontSize',20);
zlabel('Z, ECI (Km)','FontSize',20);


% Plot ECI axes
% quiver3(0, 0, 0, 1, 0, 0, 1e4, 'k', 'LineWidth', 2);    % I-axis
% quiver3(0, 0, 0, 0, 1, 0, 1e4, 'k', 'LineWidth', 2);    % J-axis
% quiver3(0, 0, 0, 0, 0, 1, 1e4, 'k', 'LineWidth', 2);    % K-axis
% text(1e4, 0, 0, 'I', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'k')
% text(0, 1e4, 0, 'J', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'k')
% text(0, 0, 1e4, 'K', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'k')
title('Lunar Gravity Assist from Launch Site',STitle,'FontSize',20)

% Plot Trajectory
plot3(I_X(1,:), I_X(2,:), I_X(3,:), 'r--', 'LineWidth', 2,'DisplayName','Inital Orbit');
plot3(I_X(1,1), I_X(2,1), I_X(3,1), 'ok', 'MarkerFaceColor', 'b','DisplayName','TLI');
text(I_X(1,1), I_X(2,1), I_X(3,1), '\leftarrow Trans-Lunar Injection Burn','DisplayName','TLI Label','FontSize',20);

plot3(T_X(1,1:2500), T_X(2,1:2500), T_X(3,1:2500), 'k', 'LineWidth', 2,'DisplayName','Transfer Orbit');

plot3(SAT_Closest_PV(1), SAT_Closest_PV(2), SAT_Closest_PV(3),'o', 'LineWidth', 2,'DisplayName','Lunar Encounter');
text(SAT_Closest_PV(1), SAT_Closest_PV(2), SAT_Closest_PV(3), '\leftarrow Lunar Encounter','FontSize',20)
% 
plot3(LU_X(1,:), LU_X(2,:), LU_X(3,:), 'r', 'LineWidth', 2,'DisplayName','Lunar Orbit');

plot3(X_FO1(1,SIMV(6):5000), X_FO1(2,SIMV(6):5000), X_FO1(3,SIMV(6):5000), 'b', 'LineWidth', 2,'DisplayName','Anti-Planet Orbit');


plot3(X_FO2(1,SIMV(7):5000), X_FO2(2,SIMV(7):5000), X_FO2(3,SIMV(7):5000), 'g', 'LineWidth', 2,'DisplayName','Planet Orbit');

plot3(X_G(1,:), X_G(2,:), X_G(3,:), 'm', 'LineWidth', 2,'DisplayName','GEO Orbit');
plot3(X_G(1,1), X_G(2,1), X_G(3,1), 'ok', 'MarkerFaceColor', 'k','DisplayName','GEO Circularisation Burn (Anit-Planet Return)');
text(X_G(1,1), X_G(2,1), X_G(3,1),'\leftarrow GEO Circularisation Burn (Anti-Planet)','FontSize',20);



plot3(X_G2(1,:), X_G2(2,:), X_G2(3,:), 'm', 'LineWidth', 2,'DisplayName','GEO Orbit');
plot3(X_G2(1,1), X_G2(2,1), X_G2(3,1), 'ok', 'MarkerFaceColor', 'k','DisplayName','GEO Circularisation Burn (Planet Return)');
text(X_G2(1,1), X_G2(2,1), X_G2(3,1),'\leftarrow GEO Circularisation Burn (Planet)','FontSize',20);

legend('FontSize',15)