%% Antenna and radio sources project %% 
% Sergio Cuevas del Valle, Brandon Escamilla, MiSE 2021 

%% Conformal array design %% 
% This file introduces an optimization problem to design a conformal antenna array form by wires. Specifically, several 
% analytical solutions are introduced. 

%% Constant continuous array of non-radiative sources 
% Directivity domain
theta = -pi/2:1e-2:pi/2; 
phi = theta; 
thetad = 0;

% Wire characteristics
I = 1; 
c = 3e8; 
f = 20e9; 
lambda = c/f; 
ds = 1e-2; 
l = 1e-1; 

% Parallell to the x axis 
alpha = 1-sin(theta).^2; 
Phi0 = 15*pi*(I*ds/lambda)^2*alpha; 

% Space factor 
beta = 2*pi/lambda;
k = 0; 
S = 2*(abs(sin(0.5*l*(beta*cos(theta-thetad)-k))))./(beta*cos(theta-thetad)-k);

% Radition intensity 
Phi = Phi0; 

for i = 1:length(Phi0)
    Phi(i) = S(i)*Phi0(i);
end

% Results 
figure(1) 
plot(theta, Phi)
ylabel('Radiation pattern')
xlabel('Angular domain')
grid on; 

%% Rectangular array in the XY plane
% Dimensions of the array 
a = 1e-1; 
b = 1e-1; 

theta = -pi/2:1e-2:pi/2;
phi = theta; 

% Space factor
F = sin(pi*a/lambda*sin(theta).*cos(phi)).^2;
F = F.*sin(pi*b/lambda*sin(theta).*cos(phi)).^2;
F = F./(pi*a/lambda*sin(theta).*cos(phi)).^2;
F = F./(pi*a/lambda*sin(theta).*cos(phi)).^2;
Phi = (15*pi*a^2*I^2/lambda^2)*F.*(1-sin(theta).^2.*sin(phi).^2);

alpha = acos(sin(theta).*sin(phi));

% Results 
figure(2) 
plot(alpha, Phi)
ylabel('Radiation pattern')
xlabel('Angular domain')
grid on;

%% Fourier transform approach 
% Define the pattern of interest
k = 0; 
u = beta*cos(theta)-k; 
L = 40; 
nx = 2^10; 
dx = L/nx;
du = [0:nx-1]' * dx - L/2;

k = zeros(nx,1);
k(1:nx/2+1) = 2*[0:nx/2]/nx;
k(nx:-1:nx/2+2) = -k(2:nx/2);
k = k*pi/dx ;

theta = acos((1/beta)*(u-k));
f = pseudostep(theta, pi/2, 0.2);

Gi = ifft(f) / dx ;
Gi = ifftshift(f);

figure
plot(theta,Gi)


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

