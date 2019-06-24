%Shooting method
%Based on Michalke 1965 & Koochesfahani_Frieler 1989

clear all;close all
%% Initial conditions

h = 0.001; %Step Size
r = 0.38; %Velocity ratio U2/U1
yoff = 8; %"Infinity" 
imax = 1+(yoff/h); 
tol = 1e-6;
lam_u = 0.45; %(U1-U2)/(U1+U2)
W = 0.6; % Wake deficit
%% Input Velocity profile

[U1,U2] = Velocity_input(yoff,imax,h,W,lam_u,r);

% Plot velocity profile
figure(1)
y = linspace(-yoff,yoff,2*imax);

U = flipud([U1(:,1);flipud(U2(:,1))]);
plot(U,y,'linewidth',2)
xlabel('U')
ylabel('Y')
set(gca,'linewidth',1,'fontsize',18)


%% Density Profile

rho = zeros(imax,2);

rho(:,1) = 2;
rho(:,2) = 0;

%% Main Routine

a = 0;
tic
for beta = 0.03:0.005:0.9
    a = a+1;
    if a==1
        alpha = beta - 0.01*1i;   %alpha = beta - 0.01*1i; alpha2 = 10*(beta  - beta*1i); work for Koochesfahani vel input lam = 0.45 W = 0.4,0.6,0.8; Tumins r = 0.38, W = 0.414 alpha 1 = 0.1274 - 0.0541*1i;alpha2 = 0.18  - 0.082*1i;
        alpha2 = 10*(beta  - beta*1i);
    end
   

%% Shooting scheme for alpha

% First guess (y1_1 & y2_1 is phi at high and low speed sides)
y0(1) = -alpha;% Boundary Conditions
y1_1 = RungeKuttaGill(y0,U1,rho,imax,-h,alpha,beta);
y2_1= RungeKuttaGill(-y0,U2,rho,imax,h,alpha,beta); 

% % Second guess 
y0(1) = -alpha2;% Boundary Conditions
y1_2 = RungeKuttaGill(y0,U1,rho,imax,-h,alpha2,beta);
y2_2= RungeKuttaGill(-y0,U2,rho,imax,h,alpha2,beta); 


%% fzero

%  alpha_comp = fzero(@(alp) matchcond(alp,U1,U2,rho,imax,h,beta,yoff),alpha);
% 
%  alpha = alpha_comp;

%% Secant method

matchcond1 = y1_1(1,imax) - y2_1(1,imax);
matchcond2 = y1_2(1,imax) - y2_2(1,imax);

alpha_a = alpha;
while (abs(matchcond1) > 1e-6)

alpha1 = alpha;
alpha = alpha - ((matchcond1)*(alpha - alpha2)/((matchcond1) - (matchcond2)));
alpha2 = alpha1;

matchcond2 = matchcond1;
y1_2 = y1_1;
y2_2 = y2_1;


% New guess
y0(1) = -alpha;% Boundary Conditions
y1_1 = RungeKuttaGill(y0,U1,rho,imax,-h,alpha,beta);
y2_1= RungeKuttaGill(-y0,U2,rho,imax,h,alpha,beta); 

matchcond1 = y1_1(1,imax) - y2_1(1,imax);


end
alpha2 = 2*alpha_a;

% ci = beta*imag(alpha)/(abs(alpha))^2;

cr = beta/real(alpha);
% alpha_final(a) = alpha;
% beta_final(a) = beta;


%% Plotting
figure(3)
plot(beta,real(alpha),'*')
hold on;
ylim([0 1])
xlabel('\beta_r')
ylabel('\alpha_r')
set(gca,'linewidth',1,'fontsize',15)


figure(4)
plot(beta,-imag(alpha),'*')
hold on;
ylim([0 0.3])
xlabel('\beta_r')
ylabel('-\alpha_i')
set(gca,'linewidth',1,'fontsize',15)

end
toc
% figure(3)
% 
% plot(beta_final,real(alpha_final),'*')
% hold on;
% ylim([0 1])
% xlabel('\beta_r')
% ylabel('\alpha_r')
% set(gca,'linewidth',1,'fontsize',15)
% figure(4)
% 
% plot(beta_final,-imag(alpha_final),'*')
% hold on;
% ylim([0 0.25])
% xlabel('\beta_r')
% ylabel('-\alpha_i')
% set(gca,'linewidth',1,'fontsize',15)
