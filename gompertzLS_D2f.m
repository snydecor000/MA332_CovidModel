function D2fvalue = gompertzLS_D2f(parameters)
%GOMPERTZLS_D2F evaluates the Gompertz Least Squares Hessian Matrix given the
%Gompertz coefficients.  
%
%   This simple function loads in the Coronavirus .mat data files, and
%   then uses the error function, the first partial derivatives of the error 
%   function, and the second partial derivatives of the error function to
%   solve for the Least Squares Gompertz Hessian given the 3 input
%   parameters
%   
%   Cory Snyder 10/30/2020

% Rename parameters to a
a = parameters;

% Load in function values
load('TotalConfirmedCasesinUS.mat');
y = TotalConfirmedCasesinUS;
load('Day.mat');
x = Day;

% Function Definitions pulled directly from GompertsLSDerivatives.mw Maple
% Worksheet (Located in this repository)
error = @(i) y(i)-a(1)*exp(-log(a(1)/a(2))*exp(-a(3)*x(i)));
derrorda1 = @(i) (exp(-a(3)*x(i))-1)*(a(1)/a(2))^(-exp(-a(3)*x(i)));
derrorda2 = @(i) -a(1)*exp(-a(3)*x(i))*(a(1)/a(2))^(-exp(-a(3)*x(i)))/a(2);
derrorda3 = @(i) -a(1)*log(a(1)/a(2))*x(i)*exp(-a(3)*x(i))*(a(1)/a(2))^(-exp(-a(3)*x(i)));
d2errorda1da1 = @(i) -exp(-a(3)*x(i))*(a(1)/a(2))^(-exp(-a(3)*x(i)))*(exp(-a(3)*x(i))-1)/a(1);
d2errorda1da2 = @(i) exp(-a(3)*x(i))*(a(1)/a(2))^(-exp(-a(3)*x(i)))*(exp(-a(3)*x(i))-1)/a(2);
d2errorda1da3 = @(i) exp(-a(3)*x(i))*x(i)*(-1+(exp(-a(3)*x(i))-1)*log(a(1)/a(2)))*(a(1)/a(2))^(-exp(-a(3)*x(i)));
d2errorda2da2 = @(i) -a(1)*exp(-a(3)*x(i))*(a(1)/a(2))^(-exp(-a(3)*x(i)))*(exp(-a(3)*x(i))-1)/a(2)^2;
d2errorda2da3 = @(i) -a(1)*x(i)*exp(-a(3)*x(i))*(a(1)/a(2))^(-exp(-a(3)*x(i)))*(log(a(1)/a(2))*exp(-a(3)*x(i))-1)/a(2);
d2errorda3da3 = @(i) -(a(1)/a(2))^(-exp(-a(3)*x(i)))*log(a(1)/a(2))*a(1)*x(i)^2*exp(-a(3)*x(i))*(log(a(1)/a(2))*exp(-a(3)*x(i))-1);     
D2fvalue = zeros(3,3);

% Assemble the Hessian
for i = 1:length(x)
    D2fvalue(1,1) = D2fvalue(1,1) + 2*(derrorda1(i))^2 + 2*error(i)*d2errorda1da1(i);
    D2fvalue(1,2) = D2fvalue(1,2) + 2*derrorda1(i)*derrorda2(i) + 2*error(i)*d2errorda1da2(i);
    D2fvalue(1,3) = D2fvalue(1,3) + 2*derrorda1(i)*derrorda3(i) + 2*error(i)*d2errorda1da3(i);
    D2fvalue(2,2) = D2fvalue(2,2) + 2*(derrorda2(i))^2 + 2*error(i)*d2errorda2da2(i);
    D2fvalue(2,3) = D2fvalue(2,3) + 2*derrorda2(i)*derrorda3(i) + 2*error(i)*d2errorda2da3(i);
    D2fvalue(3,3) = D2fvalue(3,3) + 2*(derrorda3(i))^2 + 2*error(i)*d2errorda3da3(i);
    % Hessian matrix is symetric because it doesn't matter which order you
    % take partial derivatives in
    D2fvalue(2,1) = D2fvalue(1,2);
    D2fvalue(3,2) = D2fvalue(2,3);
    D2fvalue(3,1) = D2fvalue(1,3);
end


