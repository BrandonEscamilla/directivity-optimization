%% Antenna and radio sources project %% 
% Sergio Cuevas del Valle, Brandon Escamilla, MiSE 2021 

%% Conformal array design %% 
% This file introduces an optimization problem to design a conformal antenna array form by wires. Specifically, several 
% analytical solutions are introduced. 

%% Constant continuous array of non-radiative sources 
% % Directivity domain
% % theta = -pi/2:1e-2:pi/2; 
% % phi = theta; 
% % thetad = 0;
% 
% Wire characteristics
I = 1; 
c = 3e8; 
f = 20e9; 
lambda = c/f; 
l = 1e-1; 

% % Parallell to the x axis 
% % alpha = 1-sin(theta).^2; 
% % Phi0 = 15*pi*(I*ds/lambda)^2*alpha; 
% % 
% % Space factor 
% % beta = 2*pi/lambda;
% % k = 0; 
% % S = 2*(abs(sin(0.5*l*(beta*cos(theta-thetad)-k))))./(beta*cos(theta-thetad)-k);
% % 
% % Radition intensity 
% % Phi = Phi0; 
% % 
% % for i = 1:length(Phi0)
% %     Phi(i) = S(i)*Phi0(i);
% % end
% % 
% % Results 
% % figure(1) 
% % plot(theta, Phi)
% % ylabel('Radiation pattern')
% % xlabel('Angular domain')
% % grid on; 

%% Rectangular array in the XZ plane, current paralell to the X axis
% Dimensions of the array 
a = 1e-1;                       % Maximum 1U surface dimensions
b = 1e-1; 

% Dimension length 
D = a/lambda; 

theta = -pi/2:1e-2:pi/2;

% Space factor ZX
phi = 0:1e-2:pi; 
F = sin(pi*a/lambda*sin(theta).*cos(phi)).^2;
F = F.*sin(pi*b/lambda*cos(theta)).^2;
F = F./(pi*a/lambda*sin(theta).*cos(phi)).^2;
F = F./(pi*b/lambda*cos(theta)).^2;
Phi_xz = (15*pi*a^2*I^2/lambda^2)*F.*(1-sin(theta).^2.*cos(phi).^2);

alpha = acos(sin(theta).*cos(phi));

Phi_xz = [Phi_xz Phi_xz];
alpha = [-flip(alpha) alpha];

figure(2) 
polar(alpha, sqrt(Phi_xz))
ylabel('Radiation pattern')
xlabel('Angular domain')
grid on;


% Space factor ZY
phi = -pi/2:1e-2:pi/2; 
alpha = acos(sin(theta).*sin(phi));
F = sin(pi*a/lambda*sin(theta).*sin(phi)).^2;
F = F.*sin(pi*b/lambda*cos(theta)).^2;
F = F./(pi*a/lambda*sin(theta).*sin(phi)).^2;
F = F./(pi*b/lambda*cos(theta)).^2;
Phi_yz = (15*pi*b^2*I^2/lambda^2)*F.*(1-sin(theta).^2.*sin(phi).^2);

Phi_yz = [Phi_yz Phi_yz];
alpha = [-flip(alpha) alpha];

figure(3)
polar(alpha, sqrt(Phi_yz))
ylabel('Radiation pattern')
xlabel('Angular domain')
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

figure(4) 
polar(alpha, sqrt(Phi_xz))
ylabel('Radiation pattern')
xlabel('Angular domain')
grid on;

figure(5) 
polar(alpha, sqrt(Phi_yz))
ylabel('Radiation pattern')
xlabel('Angular domain')
grid on;

% It is therefore a better option to avoid interaction between opposite
% surfaces, as directivity is lost. 

%% Single wire on an edge (parallel to the Z axis)
% Straigth wire in the Z direction and ZX plane
theta = -pi/2:1e-2:pi/2;
phi = 0:1e-2:pi; 
beta = 2*pi*b/lambda;
K = 0.999;
k = cos(pi/2*K)*beta; 

Phi_zx = 30*pi*I^2/lambda^2*((1-cos(beta*cos(theta)-k)*l)./(beta*cos(theta)-k).^2).*(cos(theta)).^2;
alpha = theta;

figure(6) 
polar(alpha, sqrt(Phi_zx))
ylabel('Radiation pattern')
xlabel('Angular domain')

% Straight wire in the Z direction and ZY plane
theta = -pi/2:1e-2:pi/2;
phi = pi/2:1e-2:-pi/2; 
beta = 2*pi*b/lambda;
K = 0.999;
k = cos(pi/2*K)*beta; 

Phi_zx = 30*pi*I^2/lambda^2*((1-cos(beta*cos(theta)-k)*l)./(beta*cos(theta)-k).^2).*(cos(theta)).^2;
alpha = theta;

figure(7) 
polar(alpha, sqrt(Phi_zx))
ylabel('Radiation pattern')
xlabel('Angular domain')

% This directly connects to the Fourier transform of the required radiation
% pattern 

%% Fourier transform approach. Array synthesis for a current distribution (wire) on a given surface
% Generate the double sided frequency response 
theta = -pi/2:1e-3:pi/2; 
L = length(theta);
n = 2^nextpow2(L);

f = pseudostep(theta, 0, 0.2);
plot(theta,f)

% The inverse transform is 
I = ifft(f,'symmetric');
i = fft(I,n);


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

