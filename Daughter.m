% Solution of 1-D Transport Equation by Bright Takyi (10878402)
% Transport Equation; s_t = -v*s_x + D*s_xx - lambda
% Initial and Boundary conditions;
% s(x,0)=s0, s(0,t)=s1, ds/dx=0
%***********************************************************************************
% Initializing Parameters
 clear
 close all
 clc
L = 10000; % lenght of river [m]
T = 2*3600 ; % simulation time [sec]
%m = 500; % number of subintervals of length of river
%n = 360; %1728; % number of subintervals of simulation time
%f = 0.1; % fraction of parent that form the daughter
load decayp; % fraction of parent that forms daughter
f = decayp;
load s1
sf = f.*s1;
dt = 1;
dx = 20;
n = T/dt;
m = L/dx;
s2max = 100;
%dx = L/m; % space step-size
%dt = 1;
%dt = T/n; % time step-size
%oneday = 24*3600;
x = linspace(0,L,m+1); % space discretization
t = linspace(0,T,n+1); % time discretization
rx = dt/(2.*dx);
rxx = dt/(dx^2);
halflife = 3*60;
%halflife = 20.8*24*3600; % halflife of U-230
decay = log(2)/halflife; % decay rate of U-230
%decay = 0;
e = 1-decay*dt;
%v = 1000.*ones(m+1,1);
v = 3.*ones(m+1,1); % wave celerity [m/s]                                        
%D1 = xlsread('C:\Users\G\Desktop\Diffusivity.xls','A4:B23');
%D=repmat(D1',m+1,1);% diffusion coefficient [m^3/s]
D = 0.003.*ones(m+1,1); % diffusion coefficient [m^3/s]
%D = 5.*ones(m+1,1);
% Setting up Q matrix, initial and boundary conditions
sd = zeros(m+1,n+1);  
% Boundary condition
sd(1,1) = 0;
sd(1,2:end)= sf(1,2:end); 
 
% initial condition
sd(:,1) = 0;
%sd(:,2:n+1) = sf(:,1:n);
sd(1,1) = 0;
so = zeros(m+1,n+1);
s2 = zeros(m+1,n+1);
%sa = 10; % initial concentration of radon-220 [mol/m^3]
%sb = 100;% boundary condition at X = 0
%sd(1,:) = 0; % initial fraction of parent that formed daughter
%sd(2:end,1) = 0; % boundary condition 
s2(:,1) = sd(:,1);
s2(1,:) = sd(1,:);
 
 
% For loop 
for k = 2:n+1   % temporal loop
    a = zeros(m+1,1);
    b = zeros(m+1,1);
    d = zeros(m+1,1);
    h = zeros(m+1,1);
    A = zeros(m,m);
    %F = zeros(m+1,1);
    for i = 2:m+1   % spatial loop
          h(i) = e*s2(i,k-1)+sd(i,k-1).*dt;
          a(i) = -0.5*(v(i-1)*rx+v(i)*rx)+0.5*(D(i-1)*rx+D(i)*rx);
          if i==m+1
              b(i) = 1+0.5*(v(i-1)*rx+v(i)*rx)+0.5*(D(i-1)+D(i)+D(i-1)+D(i))*rxx;
              d(i) = -0.5*(D(i-1)+D(i))*rxx;
          else   % for i < m+1
              b(i) = 1+0.5*(v(i+1)*rx+v(i)*rx)+0.5*(D(i+1)+D(i)+D(i-1)+D(i))*rxx;
 
              d(i) = -0.5*(D(i+1)+D(i))*rxx;
 
          end
          
        if i==2
            h(i) = h(i)-a(i)*s2(1,k);
            %F(i) = s(i,k-1)-a(i)*sa;
            A(i-1,i-1) = b(i);
            A(i-1,i) = d(i);
        elseif i==m+1
            A(i-1,i-1) = b(i);
            A(i-1,i-2) = a(i)+d(i);
           % h(i) = s(i,k-1);
        else
            A(i-1,i-2) = a(i);
            A(i-1,i-1) = b(i);
            A(i-1,i) = d(i);
            %h(i) = s(i,k-1);
        
        end
    end
    
    
    
    s2(2:end,k) = (A\h(2:end,1));
 
end
for j = 2:n+1
    so(:,j)= s2(:,j)./s2(1,j);
end
s2 = s2./s2max; % relative concentrations to plot 
 
figure
%plot(t(1:7201),s2(1,1:7201),'LineWidth',3);
plot(t(1:7201),s2(16,1:7201),'LineWidth',3);
hold on
plot(t(1:7201),s2(46,1:7201),'LineWidth',3);
plot(t(1:7201),s2(91,1:7201),'LineWidth',3);
plot(t(1:7201),s2(271,1:7201),'LineWidth',3);
%plot(t(1:7201),s2(51,1:7201),'LineWidth',3);
xlabel('Time [s]')
ylabel('Relative Concentration [s2/s2max]')
legend('x = 300 m','x = 900 m', 'x = 1800 m', 'x = 5400 m','location','East')
%title('Graph of U-230 for various distance')
 
 
figure
plot(x,s2(:,101),'LineWidth',3);
hold on 
plot(x,s2(:,301),'LineWidth',3);
plot(x,s2(:,601),'LineWidth',3);
plot(x,s2(:,1801),'LineWidth',3);
%plot(x,s2(:,361),'LineWidth',3);
xlabel('Distance from initial section [m]')
ylabel('Relative Concentration [s2/s2max]')
%title('Graph of U-230 for various days')
legend('After 100 s','After 300 s', 'After 600 s','After 1800 s') 
 
 
%figure
%mesh(x,t,s','LineWidth',3)
%xlabel('Distance [m]')
%ylabel('Time [s]')
%zlabel('Radionuclide Concentration (s)[Bq]')
%title('A 3D Graph of Pa-230 concentration')


