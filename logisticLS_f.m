function fvalue = logisticLS_f(parameters)
%LOGISTICLS_F evaluates the Logistic Least Squares value given the
%Logtistic coefficients.  
%
%   This simple function loads in the Coronavirus .mat data files, and
%   then uses the error function to solve for the Least Squares Logistic 
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

% Function Definition pulled directly from LogisticLSDerivatives.mw Maple
% Worksheet (Located in this repository)
error = @(i) y(i) - a(1) / (1 + a(2)*exp(-a(3)*x(i)));
fvalue = 0;

% Then sum the Least Squares errors
for i = 1:length(x)
    fvalue = fvalue + error(i)^2;
end

