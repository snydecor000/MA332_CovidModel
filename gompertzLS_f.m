function fvalue = gompertzLS_f(parameters)
%GOMPERTZLS_F evaluates the Gompertz Least Squares value given the
%Gompertz coefficients.  
%
%   This simple function loads in the Coronavirus .mat data files, and
%   then uses the error function to solve for the Least Squares Gompertz 
%   value given the 3 input parameters
%   
%   Cory Snyder 10/30/2020

% Rename parameters to a
a = parameters;

% Load in function values
load('TotalConfirmedCasesinUS.mat');
y = TotalConfirmedCasesinUS;
load('Day.mat');
x = Day;

% Function Definition pulled directly from GompertzLSDerivatives.mw Maple
% Worksheet (Located in this repository)
error = @(i) y(i)-a(1)*exp(-log(a(1)/a(2))*exp(-a(3)*x(i)));

fvalue = 0;

% Then sum the Least Squares errors
for i = 1:length(x)
    fvalue = fvalue + error(i)^2;
end


