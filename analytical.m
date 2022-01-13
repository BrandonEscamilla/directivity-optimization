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
theta = -pi/2:1e-2:pi/2;
beta = 2*pi*b/lambda;
K = 0.999;
k = cos(pi/2*K)*beta; 

Phi = 30*pi*I^2/lambda^2*((1-cos(beta*cos(theta)-k)*l)./(beta*cos(theta)-k).^2).*(sin(theta)).^2;
alpha = theta;

figure
polarplot(alpha, sqrt(Phi))
title('Radiation pattern for ZX wire, Z current')

% Straigth wire in the ZX plane, current parallel to the X direction (+)
phi = 0:1e-2:pi; 
Phi = 30*pi*I^2/lambda^2*((1-cos(beta*cos(theta)-k)*l)./(beta*cos(theta)-k).^2).*(cos(theta)).^2;
alpha = acos(sin(theta).*sin(phi));

figure
polarplot(alpha, sqrt(Phi))
title('Radiation pattern for ZX wire, X current')

%% Rectangular array in the ZX planes
% Angular domain
theta = -pi/2:1e-2:pi/2;

% Radiation pattern for a ZY rectangular array, current parallel to the + Y direction
phi = pi/2:1e-2:3*pi/2; 
phi = flip(phi);
F = sin(pi*a/lambda*sin(theta).*cos(phi)).^2;
F = F.*sin(pi*b/lambda*cos(theta)).^2;
F = F./(pi*a/lambda*sin(theta).*cos(phi)).^2;
F = F./(pi*b/lambda*cos(theta)).^2;
Phi_xz = (15*pi*a^2*I^2/lambda^2)*F.*(1-sin(theta).^2.*cos(phi).^2);
alpha = acos(sin(theta).*cos(phi));
Phi_xz = [Phi_xz Phi_xz];
alpha = [flip(alpha) alpha];

figure
polarplot(alpha, sqrt(Phi_xz))
title('Radiation pattern of a ZY rectangular array, +Y direction')
grid on;

% Radiation pattern for a ZY rectangular array, current parallel to the -Y direction
phi = 0:1e-2:pi; 
F = sin(pi*a/lambda*sin(theta).*cos(pi/2-phi)).^2;
F = F.*sin(pi*b/lambda*cos(theta)).^2;
F = F./(pi*a/lambda*sin(theta).*cos(pi/2-phi)).^2;
F = F./(pi*b/lambda*cos(theta)).^2;
Phi_xz = (15*pi*a^2*I^2/lambda^2)*F.*(1-sin(theta).^2.*cos(pi/2-phi).^2);
alpha = acos(sin(theta).*cos(pi/2-phi));
Phi_xz = [Phi_xz Phi_xz];
alpha = [flip(alpha) alpha];

figure
polarplot(alpha, sqrt(Phi_xz))
title('Radiation pattern of a ZX rectangular array, -X direction')
grid on;

% Radiation pattern for a ZX rectangular array, current parallel to the Z direction
phi = 0:1e-2:pi; 
F = sin(pi*a/lambda*sin(theta).*cos(phi)).^2;
F = F.*sin(pi*b/lambda*cos(theta)).^2;
F = F./(pi*a/lambda*sin(theta).*cos(phi)).^2;
F = F./(pi*b/lambda*cos(theta)).^2;
Phi_xz = (15*pi*a^2*I^2/lambda^2)*F.*(cos(theta).^2);
alpha = theta;
Phi_xz = [Phi_xz Phi_xz];
alpha = [flip(alpha) alpha];

figure
polarplot(alpha, sqrt(Phi_xz))
title('Radiation pattern of a ZX rectangular array, +Z direction')
grid on;

% Radiation pattern for a ZX rectangular array, current parallel to the -Z direction
phi = 0:1e-2:pi; 
F = sin(pi*a/lambda*sin(theta).*cos(phi)).^2;
F = F.*sin(pi*b/lambda*cos(-theta)).^2;
F = F./(pi*a/lambda*sin(theta).*cos(phi)).^2;
F = F./(pi*b/lambda*cos(-theta)).^2;
Phi_xz = (15*pi*a^2*I^2/lambda^2)*F.*(cos(-theta).^2);
alpha = theta;
Phi_xz = [Phi_xz Phi_xz];
alpha = [flip(alpha) alpha];

figure
polarplot(alpha, sqrt(Phi_xz))
title('Radiation pattern of a ZX rectangular array, +Z direction')
grid on;

%% Rectangular array in the ZY planes
% Angular domain
theta = -pi/2:1e-2:pi/2;

% Radiation pattern for a ZX rectangular array, current parallel to the + X direction
phi = 0:1e-2:pi; 
F = sin(pi*a/lambda*sin(theta).*sin(phi)).^2;
F = F.*sin(pi*b/lambda*cos(theta)).^2;
F = F./(pi*a/lambda*sin(theta).*sin(phi)).^2;
F = F./(pi*b/lambda*cos(theta)).^2;
Phi_xz = (15*pi*a^2*I^2/lambda^2)*F.*(1-sin(theta).^2.*sin(phi).^2);
alpha = acos(sin(theta).*sin(phi));
Phi_xz = [Phi_xz Phi_xz];
alpha = [flip(alpha) alpha];

figure
polarplot(alpha, sqrt(Phi_xz))
title('Radiation pattern of a ZY rectangular array, +Y direction')
grid on;

figure 
plot(alpha, sqrt(Phi_xz))
title('Radiation pattern of a ZY rectangular array, +Y direction')
grid on; 

% Radiation pattern for a ZY rectangular array, current parallel to the -Y direction
phi = pi/2:1e-2:3*pi/2;
F = sin(pi*a/lambda*sin(theta).*sin(pi/2-phi)).^2;
F = F.*sin(pi*b/lambda*cos(theta)).^2;
F = F./(pi*a/lambda*sin(theta).*sin(pi/2-phi)).^2;
F = F./(pi*b/lambda*cos(theta)).^2;
Phi_xz = (15*pi*a^2*I^2/lambda^2)*F.*(1-sin(theta).^2.*sin(pi/2-phi).^2);
alpha = acos(sin(theta).*sin(pi/2-phi));
Phi_xz = [Phi_xz Phi_xz];
alpha = [flip(alpha) alpha];

figure
polarplot(alpha, sqrt(Phi_xz))
title('Radiation pattern of a ZY rectangular array, -Y direction')
grid on;

% Radiation pattern for a ZX rectangular array, current parallel to the Z direction
phi = pi/2:1e-2:3*pi/2; 
F = sin(pi*a/lambda*sin(theta).*sin(phi)).^2;
F = F.*sin(pi*b/lambda*cos(theta)).^2;
F = F./(pi*a/lambda*sin(theta).*sin(phi)).^2;
F = F./(pi*b/lambda*cos(theta)).^2;
Phi_xz = (15*pi*a^2*I^2/lambda^2)*F.*(cos(theta).^2);
alpha = theta;
Phi_xz = [Phi_xz Phi_xz];
alpha = [flip(alpha) alpha];

figure
polarplot(alpha, sqrt(Phi_xz))
title('Radiation pattern of a ZX rectangular array, +Z direction')
grid on;

% Radiation pattern for a ZX rectangular array, current parallel to the -Z direction
phi = pi/2:1e-2:3*pi/2;  
F = sin(pi*a/lambda*sin(theta).*sin(phi)).^2;
F = F.*sin(pi*b/lambda*cos(-theta)).^2;
F = F./(pi*a/lambda*sin(theta).*sin(phi)).^2;
F = F./(pi*b/lambda*cos(-theta)).^2;
Phi_xz = (15*pi*a^2*I^2/lambda^2)*F.*(cos(-theta).^2);
alpha = theta;
Phi_xz = [Phi_xz Phi_xz];
alpha = [flip(alpha) alpha];

figure
polarplot(alpha, sqrt(Phi_xz))
title('Radiation pattern of a ZX rectangular array, +Z direction')
grid on;

%% Array of rectangular arrays in the XZ and ZY planes 
% For the XZ plane 
l = 0.1; 
phi = 0:1e-2:pi; 
Phi_xz = 4*sin(pi*l/lambda*sin(theta).*cos(phi)).^2.*Phi_xz;

% For ZY planes 
l = 0.1; 
phi = -pi/2:1e-2:pi/2; 
Phi_yz = 4*sin(pi*l/lambda*sin(theta).*sin(phi)).^2.*Phi_yz;

figure 
polarplot(alpha, sqrt(Phi_xz))
title('Radiation pattern')
grid on;

figure
polarplot(alpha, sqrt(Phi_yz))
title('Radiation pattern')
grid on;

% It is therefore a better option to avoid interaction between opposite
% surfaces, as directivity is lost. 

%% Fourier transform approach. Array synthesis for a current distribution (wire) on a given surface
% Generate the double sided frequency response 
theta = -pi/2:1e-3:pi/2; 
L = length(theta);

f = pseudostep(theta, 0, 0.2);
plot(theta,f)

% The inverse transform is 
I = ifft(f);
I = I.*conj(I)/L;
i = fft(I,L);

figure
plot(theta, I)
title('Required antenna synthesis PSD')
grid on; 

%% Optimization of a rectangular array with finite numer of elements 
% Objective function 
theta = -pi/2:1e-2:pi/2; 
phi = 0:1e-2:pi; 
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
a = array(1); 
b = array(2); 
n = array(3); 
m = array(4); 

phi = pi/2:1e-2:3*pi/2; 
phi = flip(phi);
F = sin(pi*a/(lambda*(1-1/m))*sin(theta).*cos(phi)).^2;
F = F.*sin(pi*b/(lambda*(1-1/n))*cos(theta)).^2;
F = F./(pi*a/((m-1)*lambda)*sin(theta).*cos(phi)).^2;
F = F./(pi*b/((n-1)*lambda)*cos(theta)).^2;
Phi_xz = (15*pi*a^2*I^2/lambda^2)*F.*(1-sin(theta).^2.*cos(phi).^2);
alpha = acos(sin(theta).*cos(phi));
Phi_xz = [Phi_xz Phi_xz];
alpha = [flip(alpha) alpha];

figure
polarplot(alpha, Phi_xz)

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
    F = sin(pi*a/(lambda*(1-1/m))*sin(theta).*cos(phi)).^2;
    F = F.*sin(pi*b/(lambda*(1-1/n))*cos(theta)).^2;
    F = F./(pi*a/((m-1)*lambda)*sin(theta).*cos(phi)).^2;
    F = F./(pi*b/((n-1)*lambda)*cos(theta)).^2;
    Phi = (15*pi*a^2*I^2/lambda^2)*F.*(1-sin(theta).^2.*cos(phi).^2);
    f = sqrt(Phi);

    % Compute the residual
    residual = sqrt(sum(dot(f-U,f-U,2))/length(theta));
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

