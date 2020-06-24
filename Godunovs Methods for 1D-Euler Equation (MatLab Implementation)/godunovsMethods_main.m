%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Mathematical Methods in Fluid Dynamics Assignment 2 (Part 1) - Diana Doctor %%%%%

% Program Setup
clear all
testNum = 3;
fluxType = 'R';
N = 10; % number of partitions

% Get initial data
[w,gamma,h,CFL,x,T] = godunovsMethods_getInitialData(testNum,N);

% Solve for Values of w for each time
tk = 0;
k = 1; 
W(:,:,1) = w;
t = 0;

% General Numerical Scheme
while abs(T - tk) > 1e-12
    % Solve for numerical flux H based on type
    [H,lambda_max] = godunovsMethods_getNumericalFlux(fluxType,w,gamma,N);
    % Solve for tau
    tau_k = min(CFL*h/lambda_max, .1*T);
    
    % Set time step
    if tk + tau_k > T
        tau_k = T - tk;
    end    
    tk = tau_k + tk;
    
    % General Numerical Scheme - getting all values of w at time k
    for j = 2:N+1
        w(:,j) = w(:,j) - (tau_k/h)*(H(:,j) - H(:,j-1));
    end    
    w(:,1) = w(:,2);
    w(:,N+2) = w(:,N+1);
    
    t(k) = tk;
    k = k + 1; % increases time
    W(:,:,k) = w; % storing all values of w for each point x and time t    
end

% Plotting the density, pressure, and velocity with respect to x and t
S = size(W);
Rho_Data(:,:) = W(1,2:S(2),2:S(3));
Rho_U_Data(:,:) = W(2,2:S(2),2:S(3));
U_Data(:,:) = Rho_U_Data./Rho_Data;
E_Data(:,:) = W(3,2:S(2),2:S(3));
P_Data(:,:) = (gamma - 1)*(E_Data-.5*Rho_Data.*U_Data.^2);

[K,J] = meshgrid(t,x);
% Density Plot
figure(1)
surf(J,K,Rho_Data);
title('Fluid Density ({\rho}) Behavior at Point (j) and Time (k)');
xlabel('Point (j)');
ylabel('Time (k)');
zlabel('Density ({\rho})');

% Velocity Plot
figure(2)
surf(J,K,U_Data);
title('Fluid Velocity (u) Behavior at Point (j) and Time (k)');
xlabel('Point (j)');
ylabel('Time (k)');
zlabel('Velocity (u)');

% Pressure Plot
figure(3)
surf(J,K,P_Data);
title('Fluid Pressure (p) Behavior at Point (j) and Time (k)');
xlabel('Point (j)');
ylabel('Time (k)');
zlabel('Pressure (p)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%