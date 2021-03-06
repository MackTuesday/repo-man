% Least angle regression (LAR) algorithm
% Author: Xiaohui Chen (xiaohuic@ece.ubc.ca)
% Version: 2012-Feb

function [beta, A, mu, C, c, gamma] = LARS0(X, Y, option, t)

% Least Angle Regression (LAR) algorithm.
% Ref: Efron et. al. (2004) Least angle regression. Annals of Statistics.
% option = 'lar' implements the vanilla LAR algorithm (default);
% option = 'lasso' solves the lasso path with a modified LAR algorithm.
% t -- a vector of increasing positive real numbers. If given, LARS stops and 
% returns the solution at t.
%
% Output:
% A -- a sequence of indices that indicate the order of variable inclusions;
% beta: history of estimated LARS coefficients;
% mu -- history of estimated mean vector;
% C -- history of maximal current absolute correlations;
% c -- history of current correlations;
% gamma: history of LARS step size.
% Note: history is traced by rows. If t is given, beta is just the
% estimated coefficient vector at the constraint ||beta||_1 = t.
%
% Remarks:
% 1. LARS is originally proposed to estimate a sparse coefficient vector in
% a noisy over-determined linear system. LARS outputs estimates for all
% shrinkage/constraint parameters (homotopy).
%
% 2. LARS is well suited for Basis Pursuit (BP) purpose in the real case. This lars function
% automatically terminates when the current correlations for inactive set are
% all zeros. The recovered coefficient vector is the last column of beta 
% with the *lasso* option. Hence, this function provides a fast and 
% efficient solution for the ell_1 minimization problem. 
% Ref: Donoho and Tsaig (2006). Fast solution of ell_1 norm minimization problems when the
% solution may be sparse.
    if (nargin < 4)
        t = Inf;
    endif
    
    if (nargin < 3)
        option = 'lar';
    endif

    if (strcmpi(option, 'lasso'))
        lasso = 1;
    else
        lasso = 0;
    end

    eps = 1e-10;    % Effective zero

    [n,p] = size(X);
    m = min(p,n-1); % Maximal number of variables in the final active set
    T = length(t);

    beta = zeros(1,p);
    mu = zeros(n,1);    % Mean vector
    gamma = []; % LARS step lengths
    A = [];
    Ac = 1:p;
    nVars = 0;
    signOK = 1;
    i = 0;
    mu_old = zeros(n,1);
    t_prev = 0;
    beta_t = zeros(T,p);
    ii = 1;
    tt = t;


    % LARS loop
    while (nVars < m)
        i = i+1
        c = X'*(Y-mu);  % Current correlation
        C = max(abs(c));    % Maximal current absolute correlation
        if (C < eps || isempty(t))
            break;
        endif    % Early stopping criteria
        if (1 == i)
            addVar = find(C==abs(c));
        endif
        if (signOK)
            A = [A,addVar]; % Add one variable to active set
            nVars = nVars+1;
        endif
        s_A = sign(c(A));
        Ac = setdiff(1:p,A);    % Inactive set
        nZeros = length(Ac);
        X_A = X(:,A);
        G_A = X_A'*X_A; % Gram matrix
        invG_A = inv(G_A);
        qqq = s_A'*invG_A*s_A;
        
        if (qqq <= 0)
            return;
        endif
        
        L_A = sqrt(1/qqq);
        w_A = L_A*invG_A*s_A;   % Coefficients of equiangular vector u_A
        u_A = X_A*w_A;  % Equiangular vector
        a = X'*u_A; % Angles between x_j and u_A
        beta_tmp = zeros(p,1);
        gammaTest = zeros(nZeros,2);
        if (nVars == m)
            gamma(i) = C/L_A;   % Move to the least squares projection
        else
            for j = 1:nZeros
                jj = Ac(j);
                gammaTest(j,:) = [(C-c(jj))/(L_A-a(jj)), (C+c(jj))/(L_A+a(jj))];
            endfor
            gammaTest
            [gamma(i) min_i min_j] = minplus(gammaTest);
            addVar = unique(Ac(min_i));
        endif
        beta_tmp(A) = beta(i,A)' + gamma(i)*w_A;    % Update coefficient estimates
        % Check the sign feasibility of lasso
        if (lasso)
            signOK = 1;
            gammaTest = -beta(i,A)'./w_A;
            [gamma2 min_i min_j] = minplus(gammaTest);
            if (gamma2 < gamma(i))   % The case when sign consistency gets violated
                gamma(i) = gamma2;
                beta_tmp(A) = beta(i,A)' + gamma(i)*w_A;    % Correct the coefficients
                beta_tmp(A(unique(min_i))) = 0;
                A(unique(min_i)) = [];  % Delete the zero-crossing variable (keep the ordering)
                nVars = nVars-1;
                signOK = 0;
            endif
        endif
        if (Inf ~= t(1))
            t_now = norm(beta_tmp(A),1);
            if (t_prev < t(1) && t_now >= t(1))
                % Compute coefficient estimates corresponding to a specific t
                beta_t(ii,A) = beta(i,A) + L_A*(t(1)-t_prev)*w_A';
                t(1) = [];
                ii = ii+1;
            endif
            t_prev = t_now;
        endif
        mu = mu_old + gamma(i)*u_A; % Update mean vector
        mu_old = mu;
        beta = [beta; beta_tmp'];
    endwhile

    if (1 < ii)
        noCons = (tt > norm(beta_tmp,1));
        if (0 < sum(noCons))
            beta_t(noCons,:) = repmat(beta_tmp',sum(noCons),1);
        endif
        beta = beta_t;
    endif

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the minimum and its index over the (strictly) positive part of X
% matrix
function [m, I, J] = minplus(X)

    % Remove complex elements and reset to Inf
    [I,J] = find(0~=imag(X));
    for i = 1:length(I)
        X(I(i),J(i)) = Inf;
    endfor

    X(X<=0) = Inf;
    m = min(min(X));
    [I,J] = find(X==m);

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least-Angle Regression
% Matteo Matteucci and Luigi Malago
% Politecnico di Milano
% 2012 May 25

function coeffs = LARS___XXX(in, predictors)
% predictors is a collection of column vectors


keyboard();
    in = in(:);
    num_predictors = size(predictors,2);
    
    % initialize number of iterations
    iteration_count = 0;
    
    % this vector is used to save the signs of the correlations
    s = ones(num_predictors, 1);
    % value of beta at the first iteration
    beta_a = zeros(num_predictors, 1);
    % \hat y at the first iteration
    mu_a = predictors * beta_a;
    
    % this matrix contains the \beta vector at each iteration
    B = [zeros(num_predictors,1)];
    % the t vector contains the l1-norms of \hat \beta at each iteration.
    % at the first iteration, the norm is zero, since \beta = 0
    t(1) = 0;
    
    % this vector contains the step size, at first iteration is zero
    gammas(1) = 0;
    
    % vector of current selected columns
    col = zeros(num_predictors, 1)
    
    % main loop, cycle until all variables are selected
    while (norm(col, 1) != norm(ones(num_predictors,1), 1))
        
        iteration_count = iteration_count + 1;
        printf("iteration %d\n", iteration_count);
        
        % evaluate correlations among variables and in - \hat in
        cor = predictors' * (in - mu_a);
        % initialize the vectors of selected columns
        col = zeros(num_predictors, 1);
        
        % find the most correlated variables
        [maxVal maxI] = max(abs(cor));
        for i = 1:size(cor,1)
            if (abs(abs(cor(i)) - maxVal) < 1e-10)
                % mark column as selected
                col(i) = 1;
                
                % save sign of correlation
                s(i) = sign(cor(i));
            else
                col(i) = 0;
            endif
        endfor
        printf("Selected variables at current iteration:");
        col'
keyboard();

        % build predictors_a matrix from predictors
        predictors_a = [];
        for i = 1:size(col,1)
            if (col(i) == 1)
                predictors_a = [predictors_a s(i)*predictors(:,i)];
            endif
        endfor
keyboard();
        
        % evaluate G_a
        G_a = predictors_a' * predictors_a;
        Ones_a = ones(norm(col,1), 1);
        M = Ones_a' * inv(G_a) * Ones_a;
        % the following code is used to evaluate the square matrix of M by
        % diagonalization: let M = V D inv(V), where D is diagonal, then
        % sqrt(M) = V sqrt(D) inv(V), where sqrt(D) has diagonal elements
        % that are the square roots of the diagonal elements of D.
        % first compute eigenvalues and eigenvectors of M to obtain
        % M = V D inv(V)
        [V,D] = eig(M);
        % remember that we evaluate sqrt(inv(M)) not just sqrt(M)
        A_a = inv(V * sqrt(D) * inv(V));
        w_a = A_a * inv(G_a) * Ones_a;
        u_a = predictors_a * w_a;
        
        a = predictors' * u_a;
        
        flag = 0;
        for i = 1:size(col,1)
            if (col(i) == 0)
                flag = flag + 1;
                % see formula 2.13 of B. Efron, T. Hastie, I. Johnstone,
                % and R. Tibshirani. Least angle regression. The Annals of
                % statistics, 32(2):407-499, 2004.  [1]
                el1 = (maxVal - cor(i)) / (A_a - a(i));
                el2 = (maxVal + cor(i)) / (A_a + a(i));
                if (el1 > 0 && el2 > 0)
                    temp = min(el1, el2);
                endif
                if (el1 <= 0 && el2 > 0)
                    temp = el2;
                endif
                if (el2 <- 0 && el1 > 0)
                    temp = el1;
                endif
                if (el1 < 0 && el2 < 0)
                    printf("ERROR unexpected condition\n", i);
                endif
                % gamma is the step size for the update of mu_a
                if (flag == 1)
                    gamma = temp;
                else
                    gamma = min(gamma, temp);
                endif
            endif
        endfor
keyboard();
        % check if this is the last step
        if (col == ones(num_predictors,1))
            gamma = maxVal / A_a;
        endif
        
        gamma
        
        % solve linear system: predictors beta_aTemp = mu_a, corresponds to
        % linsolve(predictors, mu_a) in Matlab
        beta_aTemp = predictors \ mu_a;
        beta_aTemp
        
        % save gamma at the current iteration
        gammas(iteration_count + 1) = gamma;
        mu_a = mu_a + gamma*u_a;
        
        %beta at current iteration
        beta_a = predictors \ mu_a
        
        % evaluate new value for t
        t(iteration_count+1) = norm(abs(beta_a), 1);
        % B contains all beta vectors at different iterations
        B = [B beta_a];
        
    endwhile

    B'
    gammas'
    printf("Number of iterations %d\n", iteration_count);
    
endfunction