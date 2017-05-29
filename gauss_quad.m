function [ GQweight, GQpoint ] = gauss_quad(n)

% Program to get the quadrature points
% and weight for Gauss-Legendre Quadrature
% Rule

syms y
% Input n: Quad pt rule
% Calculating the Pn(x)
% Legendre Polynomial
% Using recursive relationship
% P(order of polynomial, value of x)
% P(0,x)=1; P(1,x)=0;
% (i+1)*P(i+1,x)=(2*i+1)*x*P(i,x)-i*P(i-1,x)
m=n-1;
P0=1;
P1=y;
for i=1:1:m
    Pn=((2.0*i+1)*y*P1-i*P0)/(i+1.0);
    P0=P1;
    P1=Pn;
end
if n==1
    Pn=P1;
end
Pn=expand(Pn);
GQpoint=solve(vpa(Pn,32));
GQpoint=sort(GQpoint);
% Finding the weights
% Formula for weights is given at
% http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
% Equation (13)
for k=1:1:n
    P0=1;
    P1=y;
    m=n;
    % Calculating P(n+1,x)
    for i=1:1:m
        Pn=((2.0*i+1)*y*P1-i*P0)/(i+1.0);
        P0=P1;
        P1=Pn;
    end
    Pn=P1;
    GQweight(k)=vpa(2*(1-GQpoint(k)^2)/(n+1)^2/ ...
                                   subs(Pn,y,GQpoint(k))^2,32);
end
   % fprintf('Quad point rule for n=%g \n',n)
%disp(' ')
%disp('Abscissas')
%disp(GQpoint)
%disp('  ')
%disp('Weights')
%disp(GQweight)
end
