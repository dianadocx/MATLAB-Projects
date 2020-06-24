%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Mathematical Methods in Fluid Dynamics Assignment 2 (Part 2) - Diana Doctor %%%%%

% Get data from main method
function [w,gamma,h,CFL,x,T] = godunovsMethods_getInitialData(testNum,N)

% Set entry values for Initial Data Table
%              l   T    x0  rho_L    u_L        p_L      rho_R    u_R        p_R               
InitialData = [1  .2    .5  1        0          1        .125     0          .1;
               1  .2    .3  1        .75        1        .125     0          .1;
               1  .15   .5  1        -2         .4       1        2          .4;
               1  .012  .5  1        0          1000     1        0          .01;
               1  .035  .4  5.99924  19.5975    460.894  5.99942  -6.19633   46.0950;
               1  .012   .8  1        -19.59754  1000     1        -19.59749  .01];

% Set initial data for each parameters based on chosen test
row = testNum+1;

l = InitialData(row,1);
T = InitialData(row,2);
x0 = InitialData(row,3);
rho_L = InitialData(row,4);
u_L = InitialData(row,5);
p_L = InitialData(row,6);
rho_R = InitialData(row,7);
u_R = InitialData(row,8);
p_R = InitialData(row,9);

% Set common parameters
gamma = 1.4; % Poisson adiabatic constant
h = l/N; % length of each partitions
CFL = .9; % CFL constant
x = 0:h:l; % domain

% Set initial data for w_L and w_R
E_L = (p_L/(gamma-1))+(.5*rho_L*(u_L^2));
E_R = (p_R/(gamma-1))+(.5*rho_R*(u_R^2));

w_L = [rho_L, rho_L*u_L, E_L]';
w_R = [rho_R, rho_R*u_R, E_R]';

% Build w values for each point in each time starting at the initial values
w = zeros(3,N+2); 
for i=1:length(x)
    if x(i) <= x0 
        w(:,i) = w_L;
    else
        w(:,i) = w_R;
    end
end

w(:,N+2) = w_R; % For w_N+1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%