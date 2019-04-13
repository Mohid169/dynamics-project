% 4/12/2019 March-Launch interations
clc;
clear;
P_intPSI=input('Please input an initial pressure in psi: \n');
mh20launch=input('Please  type the mass at launch: \n');

%% Constants 
P_int = 6894.75729*P_intPSI; % Converting original reading to pascal 
P_ambient = 101185.27; %Atmospheric pressure in Pascal
P_tot = P_int + P_ambient;


mt = 1.7697;
mbottle = 0.0646;
rho = 1000; 

d_bot = 0.08549; 
d_neck = 0.02149; 
A = pi*(d_neck/2)^2; %Area of Neck 

gamma=1.4;
h=0;


%% Defining launch state 
C = P_tot *((mt-mh20launch)/rho)^gamma; %adiabatic expansion constant gam = 1.4; 
t_thrust = 1/(sqrt(2*C*rho^(gamma+1))*A)*-1/(gamma/2+1)*((mt-mh20launch)^(gamma/2+1)-(mt)^(gamma/2+1));
dt=t_thrust/1000; % Integration time step
t=0;
u = sqrt(2*P_int/rho); % Velocity when we are seperated from pressure vessel
v0=0; %initial velocity
v=0;

%% Code for water still in bottle
c1 = -1*gamma*C*A/rho;
c2 = rho/(2*C);
c3 = P_ambient/C;
c4 = (gamma+1)/gamma;
g=-9.81;

m(1)=mh20launch ;% lets consider the total mass to start
i=1;
while (m > 0) 
 
if m==mh20launch
    v = u*log((mh20launch+mbottle)/(m+mbottle)) + g*dt;
else 
    v = u*log((mh20launch+mbottle)/(m+mbottle)) + (g-0.5*1.225*0.82*pi*d_bot^2*v^2*v/abs(v))*dt;
end
 u = u + (c1*(c2*u^2+c3)^c4)*dt;
h = h + v*dt; 
m = m - rho*A*u*dt; 
t = t + dt; 

end 

t_thrust=t; 

%% After the water has run out (this implies there is no thrust force remaining)
g = -9.81;
maxH=h;
disp(maxH);

myCoordList=[];
while (h >= 0) 
v = v + (g + -0.5*1.225*0.82*pi*d_bot^2*v^2*v/abs(v))*dt;
t = t + dt; h= h + v*dt;

if (h > maxH) %This string of code I used to compare the current height with its previous index, so when it passes the zenith it will stop iteratiing
    maxH = h;
end 
plot(t,v,'k');
plot(t,h,'r');
end


%% Reporting results
max_height = maxH*3.28084; %convert to feet ;
flight_time = t;
disp(max_height);
max_thrust=P_int*A;
disp(max_thrust);
disp(t_thrust);

