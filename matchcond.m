function fun = matchcond(alp,U1,U2,rho,imax,h,beta,yoff)

y0(1) = exp(-alp*yoff);% Boundary Conditions
y1_1 = RungeKuttaGill(y0,U1,rho,imax,-h,alp,beta);
y2_1= RungeKuttaGill(y0,U2,rho,imax,h,alp,beta); 

fun = real(y1_1(1,imax) - y2_1(1,imax));
end