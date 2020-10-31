function fvalue = logisticLS_f(parameters)
%LOGISTICLS_F calculates the value of the Least Squares function for a 
%logistic model 
% 
%   fvalue = logisticLS_f(parameters)
%   
%   
% Cory Snyder 10/30/2020

% Rename parameters to a
a = parameters;

% Load in function values
load('TotalConfirmedCasesinUS.mat');
y = TotalConfirmedCasesinUS;
load('Day.mat');
x = Day;

% Define the function which represents error
error = @(i) y(i) - a(1) / (1 + a(2)*exp(-a(3)*x(i)));
fvalue = 0;

% Then sum the Least Squares errors
for i = 1:length(x)
    fvalue = fvalue + error(i)^2;
end

