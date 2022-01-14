%% Antenna and radio sources project %% 
% Sergio Cuevas del Valle, Brandon Escamilla, MiSE 2021 

%% Conformal array design %% 
% This file introduces an optimization problem to design a conformal antenna array form by wires. Specifically, several 
% analytical solutions are introduced. 

set_graphics();

%% Cube dimensions 
% Electromagnetic variables 
c = 3e8; 
f = 20e9;
lambda = c/f;

I = 1; 

% Cube dimensions
a = 0.1; 
b = 0.1; 
l = 0.1; 

% Dimension length 
D = a/lambda;                   % Continuity is needed

%% Single wire/continuous current distribution (monopole)
% Straigth wire in the ZX/ZY plane (symmetry), current parallel to the Z direction (+)
theta = 0:1e-2:pi;
beta = 2*pi*b/lambda;
K = 0.999;
k = cos(pi/2*K)*beta; 

Phi = 30*pi*I^2/lambda^2*((1-cos(beta*cos(theta)-k)*l)./(beta*cos(theta)-k).^2).*(sin(theta)).^2;
Phi = 2*Phi; 

figure
polarplot(theta, sqrt(Phi))
title('Radiation pattern for ZX wire, Z current')

figure
plot(theta, sqrt(Phi))
title('Radiation pattern for ZX wire, Z current')

%% Rectangular array in the ZX planes
% Radiation pattern for a ZX rectangular array, current parallel to the - X direction
phi = 0:1e-2:pi; 

% Preallocation 
Phi = zeros(length(theta), length(phi));
for i = 1:length(theta)
    for j = 1:length(phi)
        F = sin(pi*a/lambda*sin(theta(i))*cos(-phi(j)))^2;
        F = F*sin(pi*b/lambda*cos(theta(i)))^2;
        F = F/(pi*a/lambda*sin(theta(i))*cos(-phi(j)))^2;
        F = F/(pi*b/lambda*cos(theta(i)))^2;
        Phi(i,j) = (15*pi*a^2*I^2/lambda^2)*F*(1-sin(theta(i))^2*cos(-phi(j))^2);
        M = [rad2deg(theta(i)) rad2deg(phi(j)) Phi(i,j)];
    end
end

Phi = 2*Phi;
Phi = Phi(2:end,2:end);
theta = theta(2:end);

% Results 
figure
polarplot(theta, sqrt(Phi(:,end)))

figure 
plot(theta, sqrt(Phi(:,end)))
title('Radiation pattern of a ZY rectangular array, +Y direction')
grid on;

% Radiation pattern for a ZX rectangular array, current parallel to the + X direction
Phi = zeros(length(theta), length(phi));
for i = 1:length(theta)
    for j = 1:length(phi)
        F = sin(pi*a/lambda*sin(theta(i))*cos(phi(j)))^2;
        F = F*sin(pi*b/lambda*cos(theta(i)))^2;
        F = F/(pi*a/lambda*sin(theta(i))*cos(phi(j)))^2;
        F = F/(pi*b/lambda*cos(theta(i)))^2;
        Phi(i,j) = (15*pi*a^2*I^2/lambda^2)*F*(1-sin(theta(i))^2*cos(phi(j))^2);
        M = [rad2deg(theta(i)) rad2deg(phi(j)) Phi(i,j)];
    end
end

figure
polarplot(theta, sqrt(Phi(:,end)))
title('Radiation pattern of a ZX rectangular array, -X direction')

figure
plot(theta, sqrt(Phi(:,end)))
title('Radiation pattern of a ZX rectangular array, -X direction')
grid on;

%% Radiation pattern for a ZX rectangular array, current parallel to the Z direction
% Radiation pattern for a ZX rectangular array, current parallel to the Z direction
phi = 0:1e-2:pi; 

% Preallocation 
Phi = zeros(length(theta), length(phi));
for i = 1:length(theta)
    for j = 1:length(phi)
        F = sin(pi*a/lambda*sin(theta(i))*cos(phi(j)))^2;
        F = F*sin(pi*b/lambda*cos(theta(i)))^2;
        F = F/(pi*a/lambda*sin(theta(i))*cos(phi(j)))^2;
        F = F/(pi*b/lambda*cos(theta(i)))^2;
        Phi(i,j) = (15*pi*a^2*I^2/lambda^2*F*cos(theta(i))^2);
        M = [rad2deg(theta(i)) rad2deg(phi(j)) Phi(i,j)];
    end
end

Phi = 2*Phi;
Phi = Phi(2:end,2:end);
theta = theta(2:end);

% Results 
figure
polarplot(theta, sqrt(Phi(:,end)))

figure 
plot(theta, sqrt(Phi(:,end)))
title('Radiation pattern of a ZY rectangular array, +Y direction')
grid on;

%% Rectangular array in the ZY planes
% ZY plane, Z direction
% Preallocation 
Phi = zeros(length(theta), length(phi));
for i = 1:length(theta)
    for j = 1:length(phi)
        F = sin(pi*a/lambda*sin(theta(i))*sin(phi(j)))^2;
        F = F*sin(pi*b/lambda*cos(theta(i)))^2;
        F = F/(pi*a/lambda*sin(theta(i))*sin(phi(j)))^2;
        F = F/(pi*b/lambda*cos(theta(i)))^2;
        Phi(i,j) = (15*pi*a^2*I^2/lambda^2*F*cos(theta(i))^2);
        M = [rad2deg(theta(i)) rad2deg(phi(j)) Phi(i,j)];
    end
end

Phi = 2*Phi;
Phi = Phi(2:end,2:end);
theta = theta(2:end);

% Results 
figure
polarplot(theta, sqrt(Phi(:,end)))

figure 
plot(theta, sqrt(Phi(:,end)))
title('Radiation pattern of a ZY rectangular array, Z direction')
grid on;


% Radiation pattern for a ZY rectangular array, current parallel to the X direction
% Preallocation 
Phi = zeros(length(theta), length(phi));
for i = 1:length(theta)
    for j = 1:length(phi)
        F = sin(pi*a/lambda*sin(theta(i))*sin(phi(j)))^2;
        F = F*sin(pi*b/lambda*cos(theta(i)))^2;
        F = F/(pi*a/lambda*sin(theta(i))*sin(phi(j)))^2;
        F = F/(pi*b/lambda*cos(theta(i)))^2;
        Phi(i,j) = (15*pi*a^2*I^2/lambda^2*F*(1-sin(theta(i))*cos(phi(j)))^2);
        M = [rad2deg(theta(i)) rad2deg(phi(j)) Phi(i,j)];
    end
end

Phi = 2*Phi;
Phi = Phi(2:end,2:end);
theta = theta(2:end);

% Results 
figure
polarplot(theta, sqrt(Phi(:,end)))

figure 
plot(theta, sqrt(Phi(:,end)))
title('Radiation pattern of a ZY rectangular array, +Y direction')
grid on;

% Preallocation 
Phi = zeros(length(theta), length(phi));
for i = 1:length(theta)
    for j = 1:length(phi)
        F = sin(pi*a/lambda*sin(theta(i))*sin(phi(j)))^2;
        F = F*sin(pi*b/lambda*cos(theta(i)))^2;
        F = F/(pi*a/lambda*sin(theta(i))*sin(phi(j)))^2;
        F = F/(pi*b/lambda*cos(theta(i)))^2;
        Phi(i,j) = (15*pi*a^2*I^2/lambda^2*F*(1-sin(theta(i))*sin(phi(j)))^2);
        M = [rad2deg(theta(i)) rad2deg(phi(j)) Phi(i,j)];
    end
end

Phi = 2*Phi;
Phi = Phi(2:end,2:end);
theta = theta(2:end);

% Results 
figure
polarplot(theta, sqrt(Phi(:,end)))

figure 
plot(theta, sqrt(Phi(:,end)))
title('Radiation pattern of a ZY rectangular array, +Y direction')
grid on;

%% Optimization of a rectangular array with finite numer of elements 
% Objective function 
f = pseudostep(theta, pi/2, 0.2);

% Linear constraints
A = []; 
c = [];
Aeq = []; 
beq = [];

% Upper and lower bounds
lb = [1e-3 1e-3 2 2];      % Lower bound for the initial point x coordinate
ub = [a b 100 100];        % Upper bound for the initial point y coordinate

% Optimization options 
options = optimset('Display', 'off');

% Initial guess optimization
nvars = 4; 

% Integer constraints (to be on the cube mesh)
intlcon = [3,4]; 

% Nonlinear constraints 
nonlcon = [];

% Design the wire and obtain its performance
array = ga(@(array)costfunc(lambda, I, f, theta, phi, array), nvars, A, c, Aeq, beq, lb, ub, nonlcon, intlcon, options); 

% Array performance
c = array(1); 
d = array(2); 
n = array(3); 
m = array(4); 

for i = 1:length(theta)
    for j = 1:length(phi)
        F = sin(pi*a/(lambda*(1-1/m))*sin(theta(i))*sin(phi(j)))^2;
        F = F*sin(pi*b/(lambda*(1-1/n))*cos(theta(i)))^2;
        F = F/(pi*a/(lambda*(m-1))*sin(theta(i))*sin(phi(j)))^2;
        F = F/(pi*b/(lambda*(1-1/n))*cos(theta(i)))^2;
        Phi(i,j) = (15*pi*a^2*I^2/lambda^2*F*(1-sin(theta(i))*cos(phi(j)))^2);
    end
end

Phi = 2*Phi; 

% Plotting 
figure
polarplot(theta, sqrt(Phi(:,end)))

figure 
plot(theta, sqrt(Phi(:,end)))
grid on; 

% Gain and power radiated
Phip = Phip(2:end,2:end);
Phi_max = max(max(Phi));
Phip = Phip/Phi_max; 
P = Phi_max*trapz(theta(size(P(:,1),1)), Phip(:,end).*sin(theta(size(P(:,1),1))));
P = trapz(phi,P);
G = 4*pi/P;

%% Auxiliary functions
% Function to fake a step function continuously 
function [u] = pseudostep(phi, dphi, W)
    % Constants
    delta = 1e6;        % Saturation exponent 
    tol = 1e-4;         % Resolution tolerance

    % Preallocation of the step 
    u = zeros(length(phi),1); 

    % Compute the step of width W center at dphi as sigmoid function 
    for i = 1:length(phi)
        if (phi(i)-(dphi-W) > tol) && (phi(i)-(dphi+W) < tol)
            u(i) = 1/(1+exp(-delta*(dphi+W-phi(i))));
        end
    end
end

% Optimization cost function 
function [residual] = costfunc(lambda, I, U, theta, phi, array)
    % Array constants
    a = array(1); 
    b = array(2); 
    n = array(3); 
    m = array(4); 

    % Array radiation pattern
    Phi = zeros(length(theta), length(phi));
    for i = 1:length(theta)
        for j = 1:length(phi)
            F = sin(pi*a/(lambda*(1-1/m))*sin(theta(i))*sin(phi(j)))^2;
            F = F*sin(pi*b/(lambda*(1-1/n))*cos(theta(i)))^2;
            F = F/(pi*a/(lambda*(m-1))*sin(theta(i))*sin(phi(j)))^2;
            F = F/(pi*b/(lambda*(1-1/n))*cos(theta(i)))^2;
            Phi(i,j) = (15*pi*a^2*I^2/lambda^2*F*(1-sin(theta(i))*cos(phi(j)))^2);
        end
    end

    f = sqrt(Phi);

    % Compute the residual
    residual = 0; 
    for i = 1:length(phi)
        residual = residual + sqrt(sum(dot(f(:,i)-U,f(:,i)-U,2))/length(theta));
    end
end

% Some cool graphics setup
function set_graphics()
    %Set graphical properties
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
    set(groot, 'defaultAxesFontSize', 11); 
    set(groot, 'defaultAxesGridAlpha', 0.3); 
    set(groot, 'defaultAxesLineWidth', 0.75);
    set(groot, 'defaultAxesXMinorTick', 'on');
    set(groot, 'defaultAxesYMinorTick', 'on');
    set(groot, 'defaultFigureRenderer', 'painters');
    set(groot, 'defaultLegendBox', 'off');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultLegendLocation', 'best');
    set(groot, 'defaultLineLineWidth', 1); 
    set(groot, 'defaultLineMarkerSize', 3);
    set(groot, 'defaultTextInterpreter','latex');
end

