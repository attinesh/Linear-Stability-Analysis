function y = RungeKuttaGill(y0,U,rho,imax,h,alpha,beta)
%% Initialize 
i = 1;
y(1) = y0(1);



%% Runge-Kutta

while i < imax 

%Forward Euler half-step predictor

y_2_1(i) = y(i) + (0.5*h*f(y,U,rho,i,alpha,beta));


%Backward Euler half-step corrector

y_2_2(i) = y(i) +h*(((-0.5+1/sqrt(2))*f(y,U,rho,i,alpha,beta))+((1-1/sqrt(2))*f(y_2_1,U,rho,i,alpha,beta)));


%Mid-point rule full-step predictor 

y_1_1(i) = y(i) + h*(((-1/sqrt(2))*f(y_2_1,U,rho,i,alpha,beta))+((1+1/sqrt(2))*f(y_2_2,U,rho,i,alpha,beta)));


%Simpson's full-step corrector

y(i+1) = y(i) + 0.16666*h*(f(y,U,rho,i,alpha,beta)+((2-sqrt(2))*f(y_2_1,U,rho,i,alpha,beta))+((2+sqrt(2))*f(y_2_2,U,rho,i,alpha,beta))+f(y_1_1,U,rho,i,alpha,beta));


i = i+1;

end 
end