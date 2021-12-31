%% Conformal array design %% 
% Sergio Cuevas del Valle, MiSE 2021 

%% This file introduces an optimization problem to design a conformal antenna array form by up to 12 radiative elements 

%% Initial conditions and constants
% Generate the cube mesh
d = 1:5;
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

%% Optimization problem
% Design constraints 
L = 400;                % Total length of the wire 
dphi = pi/2;            % Directivity pulse center 
W = 1;                  % Width of the pulse 
method = 'fmincon';     % Solver to be used

num_nodes = 50;         % Number of nodes to optimize 
gamma = sqrt(2);        % Characteristic tesselation length
tol = 1;      % Continuity tolerance   

% Optimization design variables
params = [L dphi W tol];    

% Optimization 
[wire, cost, source, destination] = Astar_algorithm(mesh, params);

%% Results 
figure('Name','Cube 3D mesh','NumberTitle','off')
view(3)
hold on 
trisurf(T,'FaceColor',[0.8 0.8 1.0]);
alpha(0.1)
scatter3(source(1), source(2), source(3), 'k', 'filled')
scatter3(wire(2:end,1), wire(2:end,2), wire(2:end,3), 'b','filled')
scatter3(destination(1), destination(2), destination(3), 'r','filled')
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

% A^* algorithm  to obtain the needed space factor 
function [path, cost, source, destination] = Astar_algorithm(mesh, params)
    % Generate the graph and the nodes lists
    source = randi([1 size(mesh.Points,1)]);
    destination = randi([1 size(mesh.Points,1)]);
    close_list.ID = [];
    close_list.Cost = [];
    close_list.G = [];
    path.ID = [];
    path.Cost = [];

    % Initial cost estimation
    open_list.ID = source;
    open_list.Cost = Inf*ones(1,length(open_list)); 
    open_list.Cost(1) = 0;
    open_list.G = Inf;

    % Set up the optimization loop
    GoOn = true;                         % Convergence flag 

    % Explore the complete connectivity matrix of the source 
    while (GoOn)
        % Select the minimum cost node
        [current_cost, index] = sort(open_list.Cost);         
        current = open_list.ID(index(1));              % Current node
        current_cost = current_cost(1);
        current_G = open_list.G(index(1));

        path.ID = [path.ID; current];
        path.Cost = [path.Cost; current_cost];

        % Delete it from the open list 
        open_list.ID = open_list.ID(open_list.ID ~= current);
        open_list.Cost = open_list.Cost(open_list.ID ~= current);

        % Add it to the close list
        close_list.ID = [close_list.ID; current];
        close_list.Cost = [close_list.Cost; current_cost];
        close_list.G = [close_list.G; current_G];

        % Check the generated radiative pattern 
        goal_flag = false; 

        if (~goal_flag)
            % Explore the connectivity matrix 
            adjecentNodes = mesh.AdjencyMatrix(current,:);
            adjecentNodes = find(adjecentNodes == 1);
            adjecentPoints = mesh.Points(adjecentNodes.',:); 

            for i = 1:size(adjecentNodes,1)
                % Check if it has been closed already or if it hasn't been opened
                children.ID = adjecentNodes(i);
                close_node = any(close_list.ID == children.ID);
    
                if (~close_node)
                    % Directivity cost
                    children.H = sum(abs(mesh.Points(destination,:)-adjecentPoints(i,:)));

                    % Length of the wire cost
                    children.G = current_cost + sum(abs(adjecentPoints(i,:)-mesh.Points(current,:)));

                    % Total cost
                    children.F = children.G + children.H;
    
                    % Add it to the open list 
                    open_node = any(open_list.ID == children.ID);
                    cost_flag = all(children.G < close_list.G);
                    if (~open_node && cost_flag)
                        open_list.ID = [open_list.ID; children.ID];
                        open_list.Cost = [open_list.Cost; children.F];
                        open_list.G = [open_list.G; children.G];
                    end
                end
            end
    
            % Convergence 
            if (isempty(open_list.ID))
                GoOn = false;
            end
        else
            GoOn = false; 
        end
    end

    % Final output
    cost = path.Cost;
    path = mesh.Points(path.ID,:);
    source = mesh.Points(source,:); 
    destination = mesh.Points(destination,:);
end


% Function to optimize the problem 
function [wire, A, S, dP, psi] = conformalOptimization(mesh, gamma, num_nodes, params, method)
    % Sanity check on the number of nodes 
    L = params(1); 
    if (num_nodes*gamma > L)
        num_nodes = L/gamma;
    end

    % Branch the different solver methods 
    switch (method)
        case 'fmincon'
            % Linear constraints
            A = []; 
            b = []; 
            Aeq = []; 
            beq = [];
            
            % Upper and lower bounds
            lb = 1*ones(1,num_nodes);                       % Lower bound for the mesh indices
            ub = size(mesh.Points,1)*ones(1,num_nodes);     % Upper bound for the mesh indices

            % Optimization options 
            options = optimset('Display', 'on');

            % Initial guess
            sol0 = 1:num_nodes; 

            % Design the wire
            wire = fmincon(@(index)costfunc(mesh, params, index), sol0, A, b, Aeq, beq, lb, ub, @(index)nonlcon(mesh, params, index), options); 
            wire = mesh.Points(round(wire),:);

        otherwise
            error('No valid solver was selected')
    end

    % Output 
    dP = 1; 
    psi = zeros(1,100); 
    A = 1; 
    S = 1; 
end

% Cost function 
function [cost] = costfunc(mesh, params, index)
    % Constants
    L = params(1);          % Total length of the wire

    % Points in the mesh 
    r = mesh.Points(round(index),:);

    % Cost function
    cost = darclength(r);
    cost = cost(end)-L;
end

% Nonlinear constraints
function [c, ceq] = nonlcon(mesh, params, index)    
    % Points in the mesh 
    r = mesh.Points(round(index),:);

    % Constants
    tol = params(end);      % Continuity constraint tolerance 
    L = params(1);          % Total length of the wire
    
    % Continuity analysis 
    error = zeros(1,length(index));
    
    for i = 1:length(index)-1
        error(i) = norm(r(i+1,:)-r(i,:))-tol;
    end

    % Total length constraint 
    s = darclength(r);
    error(end) = s(end)-L;

    % Nonlinear constraints
    ceq = [];                   % Equality constraint
    c = error;                  % Inequality constraint
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
end

% Compute the integral of a given function using the compound trapezoidal rule
function [I] = compound_trapz(s, func)
    % Constants 
    N = length(s);              % Number of intervals to analise

    % Compound trapezoidal rule
    fi = zeros(1,length(s));    % Preallocation of the internal values

    for i = 1:length(s)
        fi = func(s);
    end

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


%% Directive radiation functions 

% Radiative pattern for a single radiative element 
function [phi] = radiative_element(A, ds, r)
    % Compute the direction angle between the given direction and the current element 
    alpha = acos(dot(ds, r)); 

    % Compute the radiative intensity 
    phi = 15*pi*A^2*sin(alpha)^2;
end

% Space factor of a rectangular array 
function [S] = sf_recarray(a, b, I, lambda, theta, phi)
    F = sin(pi*a/lambda*sin(theta)*cos(phi))^2;
    F = F*sin(pi*b/lambda*sin(theta)*cos(phi))^2;
    F = F/(pi*a/lambda*sin(theta)*cos(phi))^2;
    F = F/(pi*a/lambda*sin(theta)*cos(phi))^2;
    S = (15*pi*a^2*I^2/lambda^2)*F*(1-sin(theta)^2*cos(phi)^2);
end

% Space factor of a given array 
function [S] = space_factor(A, ds, phi)
    % Compute the space factor 
    F = 1;
    
    for i = 1:size(ds,2)
        theta = 0;                              % Compute the angle between the two elements
        L = norm(ds(:,1)-ds(:,i));              % Distance between the two elements
        eps = 2*pi*L*cos(phi);                  % Interference angle
        F = F + A(i)*exp(1i*(i*eps+theta));     % Space factor 
    end

    % Final space factor
    S = sqrt(F*conj(F));
end