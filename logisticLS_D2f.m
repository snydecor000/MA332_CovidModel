function [D2fvalue] = logisticLS_D2f(parameters)
%LOGISTICLS_D2F calculates the sum of the Least Squares errors for a logistic
%model 
%   

% Rename parameters to a
a = parameters;

% Load in function values
load('TotalConfirmedCasesinUS.mat');
y = TotalConfirmedCasesinUS;
load('Day.mat');
x = Day;

error = @(i) y(i) - a(1) / (1 + a(2)*exp(-a(3)*x(i)));
derrorda1 = @(i) -1 / (1 + a(2)*exp(-a(3)*x(i)));
derrorda2 = @(i) ( a(1)*exp(-a(3)*x(i))   ) /...
                 ( 1+a(2)*exp(-a(3)*x(i)) )^2;
derrorda3 = @(i) ( -a(1)*a(2)*x(i)*exp(-a(3)*x(i)) ) /...
                 (      1+a(2)*exp(-a(3)*x(i))     )^2;
d2errorda1da2 = @(i) (     exp(-a(3)*x(i))    ) /...
                     ( 1+a(2)*exp(-a(3)*x(i)) )^2;
d2errorda1da3 = @(i) ( -a(2)*x(i)*exp(-a(3)*x(i)) ) /...
                     (   1+a(2)*exp(-a(3)*x(i))   )^2;
d2errorda2da2 = @(i) ( -2*a(1)*(exp(-a(3)*x(i)))^2 ) /...
                     (    1+a(2)*exp(-a(3)*x(i))   )^3;
d2errorda2da3 = @(i) (2*a(1)*(exp(-a(3)*x(i)))^2*a(2)*x(i))/(1+a(2)*exp(-a(3)*x(i)))^3 -...
                     (a(1)*x(i)*exp(-a(3)*x(i)))           /(1+a(2)*exp(-a(3)*x(i)))^2;
d2errorda3da3 = @(i) (-2*a(1)*a(2)^2*x(i)^2*(exp(-a(3)*x(i)))^2)/(1+a(2)*exp(-a(3)*x(i)))^3 + ...
                     (a(1)*a(2)*x(i)^2*exp(-a(3)*x(i)))         /(1+a(2)*exp(-a(3)*x(i)))^2;     
D2fvalue = zeros(3,3);

for i = 1:length(x)
    D2fvalue(1,1) = D2fvalue(1,1) + 2*(derrorda1(i))^2;
    D2fvalue(1,2) = D2fvalue(1,2) + 2*derrorda1(i)*derrorda2(i) + 2*error(i)*d2errorda1da2(i);
    D2fvalue(1,3) = D2fvalue(1,3) + 2*derrorda1(i)*derrorda3(i) + 2*error(i)*d2errorda1da3(i);
    D2fvalue(2,2) = D2fvalue(2,2) + 2*(derrorda2(i))^2 + 2*error(i)*d2errorda2da2(i);
    D2fvalue(2,3) = D2fvalue(2,3) + 2*derrorda2(i)*derrorda3(i) + 2*error(i)*d2errorda2da3(i);
    D2fvalue(3,3) = D2fvalue(3,3) + 2*(derrorda3(i))^2 + 2*error(i)*d2errorda3da3(i);
    % This matrix is symetric 
    D2fvalue(2,1) = D2fvalue(1,2);
    D2fvalue(3,2) = D2fvalue(2,3);
    D2fvalue(3,1) = D2fvalue(1,3);
end

