function f = f(y,U,rho,i,alpha,beta)
    
% if i == 1
%     f = 0;
% else
    f = (alpha^2 + (U(i,3)+(rho(i,2)*U(i,2)/rho(i,1)))/(U(i,1)- beta/alpha)) - (rho(i,2)/rho(i,1))*y(i) - y(i)^2 ;
%       f = (alpha^2 + U(i,3)/(U(i,1)-beta/alpha)) - y(i)^2; 
% end
end