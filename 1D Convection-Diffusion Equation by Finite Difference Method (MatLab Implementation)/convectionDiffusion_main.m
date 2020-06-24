%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mathematical Methods in Fluid Dynamics Assignment 1 (Part 1) - Diana Doctor %%%

% Setup and parameters
clear all
N = 10; % number of partitions
eps = .1; % set value for epsilon
combNum = 4; % set type of combination

x0 = 0; % initial value of x
xN = 1; % final value of x
u0 = 0; % u for initial value of x
uN = 0; % u for final value of x
h = (xN-x0)/N;
uApprox = [u0 zeros(1,N-1) uN]';

% Computes Exact and Approximate values of u
x = (x0:h:xN)';
uExact = x+(1-exp(x/eps))/(exp(1/eps)-1);

% Calculate Approximated values of u
uApprox = convectionDiffusion_getuApprox(N,eps,h,combNum,uApprox);

% Check for difference in the exact and approximation of u
uDiff = uExact - uApprox;

% Plot the graph of uApprox and uExact
figure(1)
plot(x,uApprox,'b--*','LineWidth',1.2)
hold on
plot(x,uExact,'g','LineWidth',1.5)
legend('Approximate Solution','Exact Solution')
title('One-dimensional Convection-Diffusion Equation')
ylabel('u')
xlabel('x')
hold off

% Estimate order of the combinations
a = 5; % some constants
b = 2;
n = 100;
hk = zeros(1,n);

figure(2)
% For all schemes combinations
for allComb=1:5
    % Compute error difference vector
    error = zeros(1,n)';
    for k=1:n
       % Setup parameters
       Nk = a + b*k;
       hk(k) = 1/Nk; % mesh size
       xk=(x0:hk(k):xN)';
       uApprox_k = [u0 zeros(1,Nk-1) uN]';

       % Gets exact and approximate solution for each k   
       uExact_k = xk + (1-exp(xk/eps))/(exp(1/eps)-1);
       uApprox_k = convectionDiffusion_getuApprox(Nk,eps,hk(k),allComb,uApprox_k);

       % Check for error difference in the exact and approximation of u
        uDiff_k = uApprox_k - uExact_k;
        error(k) = norm(uDiff_k,'inf'); % infinity norm
    end

    % Approximate order p
    for k=1:n-1
       q(k) = log(error(k+1)/error(k))/log(hk(k+1)/hk(k));
    end
    
    txt = strcat('Combination =',' ', num2str(allComb));
    plot(q,'-.','LineWidth',2,'DisplayName',txt)
    hold on
end

% Plot the graph of the convergence of the order to p
title('Estimation of the Numerical Method Order for All Combinations')
ylabel('Order')
xlabel('Number of Partitions')
hold off
legend show

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%