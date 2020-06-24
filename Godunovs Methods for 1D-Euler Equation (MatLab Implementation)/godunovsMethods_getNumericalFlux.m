%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Mathematical Methods in Fluid Dynamics Assignment 2 (Part 3) - Diana Doctor %%%%%

% Get data from main method
function [numFlux, lambda_max] = godunovsMethods_getNumericalFlux(fluxType,w,gamma,N)

% Initialize numerical flux
numFlux = zeros(3,N+1);
lambda_max = 0;

% For a chosen numerical flux type
switch fluxType
    
    %%%%%%%%%%%% Roe Numerical Flux %%%%%%%%%%%%
    case 'R'
        for j = 1:N+1 % For each point x
            % Set variables
            w_L = w(:,j); % w at j
            w_R = w(:,j+1); % w at j+1 
            
            % Based on definition of w
            rho_L = w_L(1);
            u_L = w_L(2)/rho_L;
            E_L = w_L(3);
            
            rho_R = w_R(1);
            u_R = w_R(2)/rho_R;
            E_R = w_R(3);

            % Based on definition at the beginning
            p_L = (gamma - 1)*(E_L - .5*rho_L*(u_L^2));
            a_L = sqrt(gamma*p_L/rho_L);
            H_L = a_L^2/(gamma - 1) + .5*(u_L^2);
            
            p_R = (gamma - 1)*(E_R - .5*rho_R*(u_R^2));
            a_R = sqrt(gamma*p_R/rho_R);
            H_R = a_R^2/(gamma - 1) + .5*(u_R^2);
            
            % Define hat parameters
            rho_hat = (.5*(sqrt(rho_L) + sqrt(rho_R)))^2;
            u_hat = (sqrt(rho_L)*u_L + sqrt(rho_R)*u_R)/(sqrt(rho_L) + sqrt(rho_R));
            H_hat = (sqrt(rho_L)*H_L + sqrt(rho_R)*H_R)/(sqrt(rho_L) + sqrt(rho_R));
            E_hat = (1/gamma)*rho_hat*H_hat + ((gamma - 1)/(2*gamma))*rho_hat*(u_hat^2);
            a_hat = sqrt((gamma - 1)*(H_hat - .5*(u_hat^2)));
            
            if (gamma - 1)*(H_hat - .5*(u_hat^2)) < 0
                disp('The computation will result to a complex number');
                return
            end
            
                   
            % Compute gamma coefficients
            delta = w_R - w_L;
            gamma_hat = zeros(1,3);
            
            gamma_hat(2) = ((gamma - 1)/(a_hat^2))*((H_hat - u_hat^2)*delta(1) + u_hat*delta(2) - delta(3));
            gamma_hat(1) = (.5/a_hat)*((u_hat + a_hat)*delta(1) - delta(2) - a_hat*gamma_hat(2));
            gamma_hat(3) = delta(1) - gamma_hat(2) - gamma_hat(1);

            % By eigenvalues definition
            lambda_L = [u_L-a_L, u_L, u_L+a_L];
            lambda_R = [u_R-a_R, u_R, u_R+a_R];
            lambda_hat = [u_hat-a_hat, u_hat, u_hat+a_hat];
            
            % Compute eigenvectors
            r1 = [1 u_hat-a_hat H_hat-(a_hat*u_hat)]';
            r2 = [1 u_hat .5*(u_hat^2)]';
            r3 = [1 u_hat+a_hat H_hat+(a_hat*u_hat)]';
            
            % Compute f for w_L and w_R
            f_L = [rho_L*u_L rho_L*(u_L^2)+p_L (E_L+p_L)*u_L]'; 
            f_R = [rho_R*u_R rho_R*(u_R^2)+p_R (E_R+p_R)*u_R]';
            
            % Calculate Roe numerical flux
            if lambda_hat(2) > 0
              % Set parameters from w_L_star
              w_L_star = w_L + gamma_hat(1)*r1;
              rho_L_star = w_L_star(1);
              u_L_star = w_L_star(2)/rho_L_star;
              E_L_star = w_L_star(3);
              
              p_L_star = (gamma - 1)*(E_L_star - .5*rho_L_star*(u_L_star^2));
              a_L_star = sqrt(gamma*p_L_star/rho_L_star);
              lambda_L_star = [u_L_star - a_L_star, u_L_star, u_L_star + a_L_star];
                            
              % Set value for lambda_tilda
              if lambda_L(1) < 0 && lambda_L_star(1) > 0
                  lambda_tilda = lambda_L(1)*((lambda_L_star(1) - lambda_hat(1))/(lambda_L_star(1) - lambda_L(1)));    
              else
                  lambda_tilda = lambda_hat(1);
              end
              
              % For the numerical flux
              numFlux(:,j) = f_L + gamma_hat(1)*min(lambda_tilda,0)*r1;
              
            else % lambda_hat(2) <= 0
              % Set parameters from w_R_star
              w_R_star = w_R - gamma_hat(3)*r3;
              rho_R_star = w_R_star(1);
              u_R_star = w_R_star(2)/rho_R_star;
              E_R_star = w_R_star(3);
              
              p_R_star = (gamma - 1)*(E_R_star - .5*rho_R_star*(u_R_star^2));
              a_R_star = sqrt(gamma*p_R_star/rho_R_star);
              lambda_R_star = [u_R_star - a_R_star, u_R_star, u_R_star + a_R_star];
                            
              % Set value for lambda_tilda
              if lambda_R_star(3) < 0 && lambda_R(3) > 0
                  lambda_tilda = lambda_R(3)*((lambda_hat(3) - lambda_R_star(3))/(lambda_R(3) - lambda_R_star(3)));   
              else
                  lambda_tilda = lambda_hat(3);
              end
              
              % For the numerical flux
              numFlux(:,j) = f_R - gamma_hat(3)*max(lambda_tilda,0)*r3;
            end         
            lambda_max = max([lambda_max, lambda_hat]);
        end
        
    %%%%%%%%%%%% Vijayasundaram Numerical Flux %%%%%%%%%%%%
    case 'V'
        for j = 1:N+1 % For each point x
            % Set variables
            w_L = w(:,j); % w at j
            w_R = w(:,j+1); % w at j+1
            
            % Based on definition of w_hat
            w_hat = (w_L + w_R)/2;
            rho_hat = w_hat(1);
            u_hat = w_hat(2)/rho_hat;
            E_hat = w_hat(3);
            
            % Based on definition at the beginning
            p_hat = (gamma - 1)*(E_hat - .5*rho_hat*(u_hat^2));
            a_hat = sqrt(gamma*p_hat/rho_hat);
            
            if gamma*p_hat/rho_hat < 0
                disp('The computation will result to a complex number');
                return
            end
            
            H_hat = a_hat^2/(gamma - 1) + .5*(u_hat^2);
            
            % Compute eigenvalues of matrix A(w)
            lambda1 = u_hat - a_hat;
            lambda2 = u_hat;
            lambda3 = u_hat + a_hat;
            
            % Build Diagonal Matrix D_plus and D_minus
            D_plus = diag([max(lambda1,0), max(lambda2,0), max(lambda3,0)]);
            D_minus = diag([min(lambda1,0), min(lambda2,0), min(lambda3,0)]);
            
            % Compute eigenvectors of matrix A(w)
            r1 = [1 u_hat-a_hat H_hat-(a_hat*u_hat)]';
            r2 = [1 u_hat .5*(u_hat^2)]';
            r3 = [1 u_hat+a_hat H_hat+(a_hat*u_hat)]';
            
            % Build the matrix of eigenvectors T
            T = [r1 r2 r3];
            
            % Build Jacobi matrix A_plus and A_minus
            A_plus = (T*D_plus)/T;
            A_minus = (T*D_minus)/T;
            
            % Solve for the numerical flux
            numFlux(:,j) = A_plus*w_L + A_minus*w_R;            
            lambda_max = max([lambda_max, lambda1, lambda2, lambda3]);
        end
        
    %%%%%%%%%%%% Steger-Warming Numerical Flux %%%%%%%%%%%%
    case 'SW' 
        for j = 1:N+1 % For each point x
            % Set variables
            w_L = w(:,j); % w at j
            w_R = w(:,j+1); % w at j+1 
            
            % Based on definition of w
            rho_L = w_L(1);
            u_L = w_L(2)/rho_L;
            E_L = w_L(3);
            
            rho_R = w_R(1);
            u_R = w_R(2)/rho_R;
            E_R = w_R(3);

            % Based on definition at the beginning
            p_L = (gamma - 1)*(E_L - .5*rho_L*(u_L^2));
            a_L = sqrt(gamma*p_L/rho_L);
            H_L = a_L^2/(gamma - 1) + .5*(u_L^2);
            
            p_R = (gamma - 1)*(E_R - .5*rho_R*(u_R^2));
            a_R = sqrt(gamma*p_R/rho_R);
            
            if gamma*p_R/rho_R < 0
                disp('The computation will result to a complex number');
                return
            end
            
            H_R = a_R^2/(gamma - 1) + .5*(u_R^2);
            
            % Compute eigenvalues of matrix A(w_L) and A(w_R)
            lambda_L = [u_L-a_L u_L u_L+a_L];
            lambda_R = [u_R-a_R u_R u_R+a_R];
            
            % Build Diagonal Matrix D_plus_L and D_minus_R
            D_plus_L = diag([max(lambda_L(1),0), max(lambda_L(2),0), max(lambda_L(3),0)]);
            D_minus_R = diag([min(lambda_R(1),0), min(lambda_R(2),0), min(lambda_R(3),0)]);
            
            % Compute eigenvectors of matrix A(w_L) and A(w_R)
            r1_L = [1 u_L-a_L H_L-(a_L*u_L)]';
            r2_L = [1 u_L .5*(u_L^2)]';
            r3_L = [1 u_L+a_L H_L+(a_L*u_L)]';
            
            r1_R = [1 u_R-a_R H_R-(a_R*u_R)]';
            r2_R = [1 u_R .5*(u_R^2)]';
            r3_R = [1 u_R+a_R H_R+(a_R*u_R)]';
            
            % Build the matrix of eigenvectors T_L and T_R
            T_L = [r1_L r2_L r3_L];
            T_R = [r1_R r2_R r3_R];
            
            % Build Jacobi matrix A_plus(w_L) and A_minus(w_R)
            A_plus_L = (T_L*D_plus_L)/T_L;
            A_minus_R = (T_R*D_minus_R)/T_R;
            
            % Solve for the numerical flux
            numFlux(:,j) = A_plus_L*w_L + A_minus_R*w_R;            
            lambda_max = max([lambda_max, lambda_L, lambda_R]);
        end
        
    %%%%%%%%%%%% Van Leer Numerical Flux %%%%%%%%%%%%
    case 'VL'
        for j = 1:N+1 % For each point x
            % Set variables
            w_L = w(:,j); % w at j
            w_R = w(:,j+1); % w at j+1
            
            % Based on definition of w_hat
            w_hat = (w_L + w_R)/2;
            rho_hat = w_hat(1);
            u_hat = w_hat(2)/rho_hat;
            E_hat = w_hat(3);
            
            % Based on definition at the beginning
            p_hat = (gamma - 1)*(E_hat - .5*rho_hat*(u_hat^2));
            a_hat = sqrt(gamma*p_hat/rho_hat);
            
            if gamma*p_hat/rho_hat < 0
                disp('The computation will result to a complex number');
                return
            end
            
            H_hat = a_hat^2/(gamma - 1) + .5*(u_hat^2);
            
            % Compute eigenvalues of matrix A(w)            
            lambda_hat = [u_hat-a_hat u_hat u_hat+a_hat];
            
            % Build Diagonal Matrix D_plus and D_minus
            phi = zeros(1,3);
            delta = .5;
            for i = 1:length(lambda_hat)
                if abs(lambda_hat(i)) < delta
                   phi(i) = (lambda_hat(i)^2 + delta^2)/(2*delta); 
                else
                   phi(i) = abs(lambda_hat(i));                    
                end
            end     
            D_abs = diag(phi);
            
            % Compute eigenvectors of matrix A_abs
            r1 = [1 u_hat-a_hat H_hat-(a_hat*u_hat)]';
            r2 = [1 u_hat .5*(u_hat^2)]';
            r3 = [1 u_hat+a_hat H_hat+(a_hat*u_hat)]';
            
            % Build the matrix of eigenvectors T
            T = [r1 r2 r3];
            
            % Build Jacobi matrix A_abs
            A_abs = (T*D_abs)/T;

            % Get values of f for w_L and w_R
            rho_L = w_L(1);
            u_L = w_L(2)/rho_L;
            E_L = w_L(3);
            p_L = (gamma - 1)*(E_L - (.5*rho_L*(u_L^2)));
            
            rho_R = w_R(1);
            u_R = w_R(2)/rho_R;
            E_R = w_R(3);
            p_R = (gamma - 1)*(E_R - (.5*rho_R*(u_R^2)));
            
            f_L = [rho_L*u_L rho_L*(u_L^2)+p_L (E_L+p_L)*u_L]'; 
            f_R = [rho_R*u_R rho_R*(u_R^2)+p_R (E_R+p_R)*u_R]';
            
            % Solve for the numerical flux
            numFlux(:,j) = .5*[(f_L + f_R) - (A_abs*(w_R - w_L))];
            lambda_max = max([lambda_max, lambda_hat]);
        end     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
