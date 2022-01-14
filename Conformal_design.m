%% Conformal array design %% 
% Sergio Cuevas del Valle, MiSE 2021 

%% This file introduces an optimization problem to design a conformal antenna array form by up to 12 radiative elements 

set_graphics();

%% Initial conditions and constants
% Generate the cube mesh
d = -2:2;
[X, Y, Z] = meshgrid(d);
x = X(:);
y = Y(:);
z = Z(:);

% Delaunay triangulation
DT = delaunayTriangulation(x,y,z);
[Tfb, Xfb] = freeBoundary(DT);
T = triangulation(Tfb, Xfb);

% Calculate adjacency matrix
mesh.AdjencyMatrix = adjency_matrix(T);
mesh.Points = T.Points;
mesh.ConnectivityList = T.ConnectivityList;

figure('Name','Cube 3D mesh','NumberTitle','off')
view(3)
hold on 
trisurf(T,'FaceColor',[0.8 0.8 1.0]);
alpha(0.1)
hold off
xlabel('x','FontSize',13);
ylabel('y','FontSize',13);
zlabel('z','FontSize',13);
title('Cube 3D mesh','FontSize',13)

%% Optimization problem
% Design constraints 
L = 400;                % Total length of the wire 

% Optimization design variables
I = 1; 
lambda = 3e8/20e9;
params = [L I lambda];    

% Optimization 
[wire, Phi, dPhi] = antenna_synthesis(mesh, params);

%% Results 
figure('Name','Cube 3D mesh','NumberTitle','off')
view(3)
hold on 
trisurf(T,'FaceColor',[0.8 0.8 1.0]);
alpha(0.1)
scatter3(wire(2:end,1), wire(2:end,2), wire(2:end,3), 'b','filled')
scatter(P1(1), P1(2), P1(3));
scatter(P2(1), P2(2), P2(3));
scatter(P3(1,:), P3(2,:), P3(3,:));
scatter(P4(1,:), P4(2,:), P4(3,:));
hold off
xlabel('x','FontSize',13);
ylabel('y','FontSize',13);
zlabel('z','FontSize',13);
title('Cube 3D mesh','FontSize',13)

%% Auxiliary functions 
% Compute the adjency matrix of the mesh 
function [AdjMat] = adjency_matrix(T)
    N = size(T.Points,1);
    AdjMat = false(N);
    for i = 1:size(T.Points,1)
      for j = 1:size(T.Points,1)
          if (i ~= j)
              flag = any(T.ConnectivityList(T.ConnectivityList(:,1) == i & T.ConnectivityList(:,2) == j));
              flag = flag || any(T.ConnectivityList(T.ConnectivityList(:,1) == i & T.ConnectivityList(:,3) == j));
              flag = flag || any(T.ConnectivityList(T.ConnectivityList(:,2) == i & T.ConnectivityList(:,3) == j));
              flag = flag || any(T.ConnectivityList(T.ConnectivityList(:,2) == i & T.ConnectivityList(:,1) == j));
              flag = flag || any(T.ConnectivityList(T.ConnectivityList(:,3) == i & T.ConnectivityList(:,1) == j));
              flag = flag || any(T.ConnectivityList(T.ConnectivityList(:,3) == i & T.ConnectivityList(:,2) == j));
              AdjMat(i,j) = flag;
          else
            AdjMat(i,j) = false;
          end
      end
    end
end

% Function to optimize the problem 
function [wire, Phi, dPhi] = antenna_synthesis(mesh, params, phi)
    % Linear constraints
    A = []; 
    b = []; 
    Aeq = []; 
    beq = [];
    
    % Upper and lower bounds
    lb = [1 0 0];                         % Lower bound for the initial point x coordinate
    ub = [size(mesh.Points,1) 0.1 10];    % Upper bound for the initial point y coordinate

    % Optimization options 
    options = optimset('Display', 'off');

    % Initial guess optimization
    nvars = 3; 

    % Integer constraints (to be on the cube mesh)
    intlcon = 1; 

    % Nonlinear constraints 
    nonlcon = [];

    % Design the wire and obtain its performance
    PopSize = 10;
    MaxGenerations = 1; 
    options = optimoptions(@ga,'PopulationSize', PopSize, 'MaxGenerations', MaxGenerations);
    wire = ga(@(vars)costfunc(mesh, params, vars), nvars, A, b, Aeq, beq, lb, ub, nonlcon, intlcon, options); 
    dPhi = wire(end-1);

    theta = -pi/2:1e-2:pi/2; 
    F = pseudostep(theta, pi/2, wire(end));
    [wire, Phi] = wirebredth(mesh, wire(1), params, wire(end), phi, theta, F);
end

% Cost function 
function [cost] = costfunc(mesh, params, vars)
    % Constants 
    L = params(1); 
    
    % Angle domain 
    phi = 0:1e-2:pi; 
    theta = -pi/2:1e-2:pi/2; 

    % Objective function (pseudo unit step radiation pattern)
    dTheta = vars(end);
    F = pseudostep(theta, pi/2, dTheta);

    % Generate the wire and assess its performance  
    [r, f] = wirebredth(mesh, vars(1), params, vars(end), phi, theta, F);

    % Cost function
    cost = sqrt(sum(dot(f-F,f-F,2))/length(theta));
    s = darclength(r);
    cost = cost+(s-L);
end

% Heuristic search of the desired wire
function [wire, Phi] = wirebredth(mesh, source_index, params, k, phi, theta, F)
    % Generate the graph and the nodes lists
    open_list.ID = 1:size(mesh.Points,1);
    close_list.ID = [];
    path.ID = [];
    path.Points = [];
    path.Cost = []; 

    current = source_index; 

    % Set up the optimization loop
    GoOn = true;                         % Convergence flag 

    % Constants 
    tol = 1e-3; 

    % Explore the complete connectivity matrix of the source 
    while (GoOn)
        % Update the graph lists
        path.ID = [path.ID; current];
        path.Points = mesh.Points(path.ID,:);
        open_list.ID = open_list.ID(open_list.ID ~= current);
        close_list.ID = [close_list.ID; current];

        % Check the generated radiative pattern error
        f = radiation_pattern(params, path.Points, k, phi, theta);
        cost = sqrt(sum(dot(f(:,1)-F,f(:,1)-F,2))/length(theta));
        path.Cost = [path.Cost; cost]; 

        if (cost > tol)
            % Explore the connectivity matrix 
            adjecentNodes = mesh.AdjencyMatrix(current,:);
            adjecentNodes = find(adjecentNodes == 1);
            adjecentPoints = mesh.Points(adjecentNodes.',:); 

            f_cost = zeros(length(adjecentNodes));
            for i = 1:size(adjecentNodes,2)
                % Check if it has been closed already or if it hasn't been opened
                close_node = any(close_list.ID == adjecentNodes(i));
 
                % Directivity cost
                if (~close_node)
                    f = radiation_pattern(params, [path.Points; adjecentPoints(i,:)], k, phi, theta);
                    f_cost(i) = sqrt(sum(dot(f(:,1)-F,f(:,1)-F,2))/length(theta));
                else
                    f_cost(i) = Inf; 
                end
            end

            % Select the minimum cost node
            [~, index] = sort(f_cost);
            current = adjecentNodes(index(1));

            % Convergence 
            if (isempty(open_list.ID) || size(open_list.ID,1) < 10)
                GoOn = false;
                path.Cost = [path.Cost; cost]; 
            end
        else
            GoOn = false; 
        end
    end

    % Final output
    Phi = radiation_pattern(params, path.Points, k, phi, theta);
    wire = path.Points;
end

% Compute the radiation pattern through the electrical field
function [Phi] = radiation_pattern(params, r, k, phi, theta)
    % Constants 
    I = params(end-1);
    lambda = params(end);
    beta = 2*pi/lambda; 

    % Preallocation 
    Phi = zeros(length(theta), length(phi)); 

    % Position computations 
    s = darclength(r);
    r = r(end,:);

    % Main computation
    for i = 1:length(theta)
        for j = 1:length(phi)
            % Compute the electric field
            r = r/norm(r); 
            u = [cos(phi(j))*sin(theta(i)); sin(phi(j))*sin(theta(i)); cos(theta(i))];
            alpha = acos(dot(r,u));
            Et = compound_trapz(s, electric_field(beta, k, alpha, s));
            Et = Et*1i*60*pi*(I/lambda)*sin(alpha);

            % Compute the radiation pattern 
            Phi(i,j) = (1/(240*pi))*Et.*conj(Et);
        end
    end
end

% Compute the integral of the electric field
function [S] = electric_field(beta, k, alpha, s)
    S = exp(1i*(beta*cos(alpha)-k)*s);
end

% Compute the discrete arc length of a 3D-curve 
function [s] = darclength(r)
    % Preallocation of the arclength 
    si = zeros(size(r,1)-1,1); 

    % Compute the point-wise arclength 
    for i = 2:size(r,1)-1
        si(i) = norm(r(i,:)-r(i-1,:));
    end

    % Compute the accumulated arclength
    s = si;
    for i = 1:size(si,1)
        s(i) = sum(si(1:i)); 
    end

    if (isempty(s))
        s = 0; 
    end
end

% Compute the integral of a given function using the compound trapezoidal rule
function [I] = compound_trapz(s, func)
    % Constants 
    N = length(s);              % Number of intervals to analise

    % Compound trapezoidal rule
    fi = func; 

    % Complete integral 
    I = (fi(1)+fi(end))/2 + sum(fi(2:end-1));
    I = (s(end)-s(1))/N*I;
end

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

% Function to cube the taxicab distance between any two points on the
% cube's surface. Not working completely. DANGER: recursive functions in
% here!
function [d, P3, P4] = cube_distance(P1, P2, L)
    % Check if they are on the same face computing how many components are different
    index(1,:) = abs(P1) == L;
    index(2,:) = abs(P2) == L;

    I = eye(3);

    % Handle edge points
    if (sum(index(1,:)) > 1)
        a = find(index(1,:) == 1); 
        index(1,:) = false(1,3); 
        index(1,a(1)) = true;
    end

    if (sum(index(2,:)) > 1)
        a = find(index(2,:) == 1); 
        index(2,:) = false(1,3); 
        index(2,a(1)) = true;
    end

    N(:,1) = sign(P1(index(1,:)).').*I(:,index(1,:));
    N(:,2) = sign(P2(index(2,:)).').*I(:,index(2,:));

    if ((dot(N(:,1),N(:,2)) == 0) || (dot(N(:,1),N(:,2)) == -1))
        % Compute the intersection between the edges and the unfolding points
        plane = (abs(P1) ~= L);

        P3 = repmat(P1,1,4);
        k = 1; 
        for i = 1:length(plane)
            if (plane(i))
                P3(i,k) = L; 
                P3(i,k+2) = -L;
                k = k+1;
            end
        end

        % Check if the edge shares surface with the destination
        new_index = zeros(1,size(P3,2)); 

        for i = 1:size(P3,2)
            Paux = P3(:,i);
            index(3,:) = abs(Paux) == L;

            Naux = sign(Paux(index(3,:)).').*I(:,index(3,:));
            if (size(Naux,2) ~= 3)
                Naux = [Naux zeros(3,1)];
            end


            if all(Naux(:,1) == Naux(:,2)) || all(Naux(:,2) == N(:,2)) || all(Naux(:,3) == N(:,2))
                new_index(i) = 1; 
            else
                new_index(i) = 0; 
            end
        end

        new_index = logical(new_index);

        if (any(new_index))
            P4 = P3(:,new_index);
            d = sum(abs(P2-P4)) + sum(abs(P4-P1));
            d = min(d);
        else
            d = zeros(1,size(P3,2));
            P4 = P3; 
            P4(~plane,:) = P4(~plane,:)-sign(P4(~plane,:))*1;

            for i = 1:length(d)
                d(i) = sum(abs(P4(:,i)-P3(:,i))) + cube_distance(P4(:,i), P2, L) + sum(abs(P3(:,i)-P1));
            end

            [d, index] = sort(d);
            P4 = P4(:,index(1));
            d = d(1);
        end
    else
        d = sum(abs(P2-P1));        % Taxicab distance betweeen two points on the same surface
        P3 = P1; 
        P4 = P2;
    end
end

% Data regression using Chebyshev polynomials of the first kind
function [cn, W] = chebyshev_coefficients(x, y, order)
    %Sanity check on the data dimensions
    if (size(y,1) == 1)
        y = y.';
    end
    
    %Constants 
    a = x(1);                                        %Initial point in the domain
    b = x(end);                                      %Final point in the domain
    
    %Project the data on the interval [-1, 1]
    u = (2*x-(b+a))/(b-a);                           %New argument
    
    %Main computation 
    Tn = zeros(length(u), order);                    %Preallocation of the Chebyshev polynomial
    for i = 1:length(u) 
        Tn(i,:) = chebyshev('first', order, u(i));   %Compute the Chebyshev polynomials
    end
    
    %Sampling distribution characteristics
    W = diag(ones(1,length(y)));
    
    %Chebyshev coefficients by least-squares
    cn = (Tn.'*W*Tn)^(-1)*(Tn.'*W*y);
    cn = cn.';
end

% Generate a Chebyshev basis
function [Pn] = chebyshev(kind, order, u)
    %Preallocation of the polynomials 
    Pn = zeros(order,1); 

    %Main computation 
    switch (kind)
        case 'first'
            Pn(1) = 1;                    %Initialization of the Chebyshev polynomials of the first kind
            Pn(2) = u;                    %Initialization of the Chebyshev polynomials of the first kind
        case 'second'
            Pn(1) = 1;                    %Initialization of the Chebyshev polynomials of the second kind
            Pn(2) = 2*u;                  %Initialization of the Chebyshev polynomials of the second kind
        otherwise
            error('No valid kind of polynomials was selected'); 
    end
  
    for i = 2:order-1
        Pn(i+1) = 2*u*Pn(i)-Pn(i-1);      %Chebyshev polynomials
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