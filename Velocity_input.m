function [U1,U2] = Velocity_input(yoff,imax,h,W,lam_u,r)
% Lower

U1 = zeros(imax,3);
for i = 1:imax

% U1(i,1) = 0.5*((1+r)+((1-r)*tanh(yoff-((i-1)*h)))); %Tumin's suggestion
% U1(i,1) = 0.5*(1 + tanh(yoff-((i-1)*h))); % Michalke
U1(i,1) = (1 + lam_u*tanh(yoff-((i-1)*h)))-W.*exp(-0.6931.*(yoff-((i-1)*h)).^2); %Koochesfahani
% U1(i,1) = 0.5*((1+r)+((1-r)*tanh(yoff-((i-1)*h)))) -W.*exp(-0.6931.*(yoff-((i-1)*h)).^2); %Tumin's suggestion + Wake like
end 

for i = 2:imax-1

U1(i,2) = (-U1(i-1,1)+U1(i+1,1))/(2*h);
U1(i,3) = (U1(i-1,1)-2*U1(i,1)+U1(i+1,1))/(h^2);

end  

U1(1,2) = (-3*U1(1,1)+4*U1(2,1)-U1(3,1))/(2*h);
U1(1,3) = (U1(1,1)-2*U1(2,1)+U1(3,1))/(h^2);
U1(imax,2) = (U1(imax-2,1)-4*U1(imax-1,1)+3*U1(imax,1))/(2*h);
U1(imax,3) = (U1(imax-2,1)-2*U1(imax-1,1)+U1(imax,1))/(h^2);



% Upper
U2 = zeros(imax,3);
for i = 1:imax

% U2(i,1) = 0.5*((1+r)+((1-r)*tanh(-yoff+((i-1)*h)))); %Tumin's suggestion

% U2(i,1) = 0.5*(1+tanh(-yoff + ((i-1)*h)));% Michalke
U2(i,1) = (1 + lam_u*tanh(-yoff+((i-1)*h)))-W.*exp(-0.6931.*(-yoff+((i-1)*h)).^2); %Koochesfahani

% U2(i,1) = 0.5*((1+r)+((1-r)*tanh(-yoff+((i-1)*h)))) -W.*exp(-0.6931.*(-yoff+((i-1)*h)).^2); %Tumin's suggestion + Wake like
end 

for i = 2:imax-1

U2(i,2) = (-U2(i-1,1)+U2(i+1,1))/(2*h);
U2(i,3) = (U2(i-1,1)-2*U2(i,1)+U2(i+1,1))/(h^2);

end  

U2(1,2) = (-3*U2(1,1)+4*U2(2,1)-U2(3,1))/(2*h);
U2(1,3) = (U2(1,1)-2*U2(2,1)+U2(3,1))/(h^2);
U2(imax,2) = (U2(imax-2,1)-4*U2(imax-1,1)+3*U2(imax,1))/(2*h);
U2(imax,3) = (U2(imax-2,1)-2*U2(imax-1,1)+U2(imax,1))/(h^2);



end