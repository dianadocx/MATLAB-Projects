%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mathematical Methods in Fluid Dynamics Assignment 1 (Part 2) - Diana Doctor %%%

%  Get data from main method
function uApprox = convectionDiffusion_getuApprox(N,eps,h,combNum,uApprox)

% Approximations for 1st derivative of u based on combination
switch combNum
    
    % 1 : Combination (1) + (6)
    case 1
        d0 = -eps-h;
        d1 = 2*eps+h;
        d2 = -eps;
        b = ones(1,N-1)'* h^2; 
        
        % Create coefficient matrix A for Au = b
        for i = 1:N-1
            for j = 1:N-1
                if i == j
                    A(i,j) = d1;
                elseif i == j+1
                    A(i,j) = d0;
                elseif i == j-1
                    A(i,j) = d2;
                end
            end
        end
        
    % 2 : Combination (2) + (6)
    case 2
        d0 = -eps;
        d1 = 2*eps-h;
        d2 = -eps+h;
        b = ones(1,N-1)'* h^2; 
        
        % Create coefficient matrix A for Au = b
        for i = 1:N-1
            for j = 1:N-1
                if i == j
                    A(i,j) = d1;
                elseif i == j+1
                    A(i,j) = d0;
                elseif i == j-1
                    A(i,j) = d2;
                end
            end
        end
    
    % 3 : Combination (3) + (6)
    case 3
        d0 = -2*eps-h;
        d1 = 4*eps;
        d2 = h-2*eps;
        b = ones(1,N-1)'* 2*h^2; 
        
        % Create coefficient matrix A for Au = b
        for i = 1:N-1
            for j = 1:N-1
                if i == j
                    A(i,j) = d1;
                elseif i == j+1
                    A(i,j) = d0;
                elseif i == j-1
                    A(i,j) = d2;
                end
            end
        end
        
    % 4 : Combination (4) + (6)
    case 4
        d0 = h;
        d1 = -2*eps-4*h;
        d2 = 4*eps+3*h;
        d3 = -2*eps;
        b = ones(1,N-1)'* 2*h^2;

        % Create coefficient matrix A for Au = b
        A(1,:) = [4*eps h-2*eps zeros(1,N-3)];
        
        for i = 2:N-1
            for j = 1:N-1
                if i == j+1
                    A(i,j) = d1;
                elseif i == j+2
                    A(i,j) = d0;
                elseif i == j
                    A(i,j) = d2;
                elseif i == j-1
                    A(i,j) = d3;
                end
            end
        end
       
    % 5 : Combination (5) + (7)
    case 5
        d0 = -4*h;
        d1 = eps+18*h;
        d2 = -16*eps-36*h;
        d3 = 30*eps+22*h;
        d4 = -16*eps;
        d5 = eps;
        b = ones(1,N-1)'* 12*h^2;

        % Create coefficient matrix A for Au = b
        A(1,:) = [20*eps-6*h    12*h-6*eps  -4*eps-2*h  eps  zeros(1,N-5)];
        A(2,:) = [-16*eps-12*h  30*eps+6*h  4*h-16*eps  eps  zeros(1,N-5)];
                
        for i = 3:N-2
            for j = 1:N-1
                if i == j+2
                    A(i,j) = d1;
                elseif i == j+3
                    A(i,j) = d0;
                elseif i == j+1
                    A(i,j) = d2;
                elseif i == j
                    A(i,j) = d3;
                elseif i == j-1
                    A(i,j) = d4;
                elseif i == j-2
                    A(i,j) = d5;
                end
            end
        end  
        A(N-1,:) = [zeros(1,N-5) eps-4*h -4*eps+18*h -6*eps-36*h 20*eps+22*h];        
end

% Calculate Approximated values of u
uApprox(2:N) = A\b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%