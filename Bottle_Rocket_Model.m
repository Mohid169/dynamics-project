% 3/6/2019 March-Launch interations
clc;
clear;
P_intPSI=input('Please input an initial pressure in psi: \n');

% Constants 
P_int = 6894.75729*P_intPSI; % Converting original reading to pascal 
P_ambient = 101185.27; %Atmospheric pressure in Pascal
P_tot = P_int + P_ambient;

mt = 1.7697;
mh20launch = 0.552; 
mbottle = 0.0646;
rho = 1000; 

d_bot = 0.08549; 
d_neck = 0.02149; 
A = pi*(d_neck/2)^2; %Area of Neck 

gamma=1.4;
h=0;


C = P_tot *((mt-mh20launch)/rho)^gamma; %adiabatic expansion constant gam = 1.4; 
t_thrust = 1/(sqrt(2*C*rho^(gamma+1))*A)*-1/(gamma/2+1)*((mt-mh20launch)^(gamma/2+1)-(mt)^(gamma/2+1));
dt=t_thrust/1000; % Integration time step
t=0;
u = sqrt(2*P_int/rho); % Velocity when we are seperated from pressure vessel
v0=0; %initial velocity
v=0;

% I kept making mistakes, so I just broke up the formula
c1 = -1*gamma*C*A/rho;
c2 = rho/(2*C);
c3 = P_ambient/C;
c4 = (gamma+1)/gamma;
g=-9.81;
m=mh20launch ;% lets consider the total mass to start 

% Thrust Force present
while (m > 0) 
v = u*log((mh20launch+mbottle)/(m+mbottle)) - g*dt ;
 u = u + (c1*(c2*u^2+c3)^c4)*dt;
h = h + v*dt; 
m = m - rho*A*u*dt; 
t = t + dt; 
end 


% After the water has run out (this implies there is no thrust force remaining)
g = -9.81;
maxH=h;
while (h >= 0) 
v = v + (g + -0.5*1.225*0.82*pi*d_bot^2*v^2*v/abs(v))*dt;
t = t + dt; h= h + v*dt; 
if (h > maxH) %This string of code I used to compare the current height with its previous index, so when it passes the zenith it will stop iteratiing
    maxH = h;
end 
end

max_height = maxH*3.28084; %convert to feet ;
flight_time = t;
%