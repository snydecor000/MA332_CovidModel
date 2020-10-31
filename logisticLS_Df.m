function [Dfvalue] = logisticLS_Df(parameters)
%LOGISTICLS_DF evaluates the Logistic Least Squares Gradient Matrix given the
%Logistic coefficients.  
%
%   This simple function loads in the Coronavirus .mat data files, and
%   then uses the error function and the partial derivatives of the error 
%   function to solve for the Least Squares Logistic Gradient Matrix 
%   given the 3 input parameters
%   
%   Cory Snyder 10/30/2020

% Rename parameters to a
a = parameters;

% Load in function values
load('TotalConfirmedCasesinUS.mat');
y = TotalConfirmedCasesinUS;
load('Day.mat');
x = Day;

% Function Definitions pulled directly from LogisticLSDerivatives.mw Maple
% Worksheet (Located in this repository)
error = @(i) y(i) - a(1)/(1 + a(2)*exp(-a(3)*x(i)));
derrorda1 = @(i) -1/(1 + a(2)*exp(-a(3)*x(i)));
derrorda2 = @(i) (a(1)*exp(-a(3)*x(i)))/(1+a(2)*exp(-a(3)*x(i)))^2;
derrorda3 = @(i) (-a(1)*a(2)*x(i)*exp(-a(3)*x(i)))/(1+a(2)*exp(-a(3)*x(i)))^2;

Dfvalue = zeros(3,1);

% Assemble the Gradient Matrix
for i = 1:length(x)
    Dfvalue(1) = Dfvalue(1) + 2*error(i)*derrorda1(i);
    Dfvalue(2) = Dfvalue(2) + 2*error(i)*derrorda2(i);
    Dfvalue(3) = Dfvalue(3) + 2*error(i)*derrorda3(i);
end

