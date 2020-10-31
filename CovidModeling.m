%   MA332_CovidModel
%   10/30/2020
%   
%   Cory Snyder
%% Polynomial Fit 
% Sets up a simple linear solver to fit the COVID-19 data to a 5th order 
% polynomial.  

% Load in variables
load('TotalConfirmedCasesinUS.mat');
y = TotalConfirmedCasesinUS;
load('Day.mat');
x = Day;

% Setup Matrix
A = zeros(length(x),6);
for i = 1:length(x)
    A(i,1) = 1;
    A(i,2) = i;
    A(i,3) = i^2;
    A(i,4) = i^3;
    A(i,5) = i^4;
    A(i,6) = i^5;
end

% Solve for the equation coeffs
alpha = linsolve(A'*A,A'*y);

% Assemble the Polynomial fucntion based on the coefficients
f = @(x) alpha(6)*x.^5 + alpha(5)*x.^4 + alpha(4)*x.^3 + alpha(3)*x.^2 + alpha(2)*x + alpha(1);
xx = 1:0.1:length(x);
% 60 Day Prediction
xx2 = x(end):0.1:x(end)+60;

% Plot it
figure(1);
hold on;
grid on;
plot(Day,TotalConfirmedCasesinUS, '-o');
plot(xx,f(xx),'LineWidth',2);
plot(xx2,f(xx2),'--','LineWidth',1.5);
title('Polynomial Fit: Total Confirmed Cases in the US');
xlabel('Days');
ylabel('Confirmed Cases');
legend('Original Data','5th Order Polynomial','2 Month Forecast','Location','southeast');
hold off;

%% Logistic Fit
% Load in variables
load('TotalConfirmedCasesinUS.mat');
y = TotalConfirmedCasesinUS;
load('Day.mat');
x = Day;

% Use Newton Optimization and Least Squares Logistic Function,Gradient,and
% Hessian to find suitable parameters
[abest,fbest,itr,status] = NewtonOpt(@logisticLS_f,@logisticLS_Df,@logisticLS_D2f,[9000000;100;0.025],1,1,100,3);
%[9000000;105;0.025]
% Assemble the Logistic model based on the coefficients from NewtonOpt
f = @(x) abest(1) ./ (1 + abest(2).*exp(-abest(3).*x));
xx = 1:0.1:length(Day);
% 60 Day Prediction
xx2 = Day(end):0.1:Day(end)+60;

% Plot it
figure(2);
hold on;
grid on;
plot(Day,TotalConfirmedCasesinUS, '-o');
plot(xx,f(xx),'LineWidth',1.5);
plot(xx2,f(xx2),'--','LineWidth',1.5);
title('Logistic Fit: Total Confirmed Cases in the US');
xlabel('Days');
ylabel('Confirmed Cases');
legend('Original Data','Logistic Curve','2 Month Forecast','Location','southeast');
hold off;

%% Gompertz Fit
% Load in variables
load('TotalConfirmedCasesinUS.mat');
y = TotalConfirmedCasesinUS;
load('Day.mat');
x = Day;

% Use Newton Optimization and Least Squares Logistic Function,Gradient,and
% Hessian to find suitable parameters
[abest,fbest,itr,status] = NewtonOpt(@gompertzLS_f,@gompertzLS_Df,@gompertzLS_D2f,[13000000;8000;0.01],1,1,2,3);

% Assemble the Gompertz model based on the coefficients from NewtonOpt
f = @(x) abest(1)*exp(-log(abest(1)/abest(2))*exp(-abest(3)*x));
xx = 1:0.1:length(Day);
% 60 Day Prediction
xx2 = Day(end):0.1:Day(end)+60;

% Plot it
figure(3);
hold on;
grid on;
plot(Day,TotalConfirmedCasesinUS, '-o');
plot(xx,f(xx),'LineWidth',1.5);
plot(xx2,f(xx2),'--','LineWidth',1.5);
title('Gompertz Fit: Total Confirmed Cases in the US');
xlabel('Days');
ylabel('Confirmed Cases');
legend('Original Data','Gompertz Curve','2 Month Forecast','Location','southeast');
hold off;