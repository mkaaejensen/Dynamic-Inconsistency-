function [opts] = algorithm_call(a, b, n, MaxIt, innerstep, opts, RESTART, dissaving, dis, p, u, tol, f, offset, dsp, dyncon, sm, tola)
% Two-step best response CMA-ES algorithm for dynamic decision problems
%
% Inputs:
%   a, b      - lower and upper bounds of the grid
%   n         - grid resolution (number of subintervals)
%   MaxIt     - maximum number of iterations
%   innerstep - initial step size for CMA-ES
%   opts      - initial policy vector
%   RESTART   - policy iteration control (1=iterative)
%   dissaving - speed optimization flag
%   dis       - discount sequence (function handle)
%   p         - probability parameter
%   u         - utility function handle
%   tol       - tolerance parameter
%   f         - income function handle
%   offset    - offset function handle
%   dsp       - display flag
%   dyncon    - dynamic constraint for initialization
%   sm        - threshold for grid points to preserve from dyncon
%   tola      - convergence tolerance
%
% Output:
%   opts      - updated policy vector after algorithm convergence

    %% Initialization
    nVar = n + 1;                           % number of variables
    VarSize = [1, nVar];                    % variable dimensions
    G = linspace(a, b, nVar);               % grid points
    initial_innerstep = innerstep;          % store initial step size

    % Variable bounds
    VarMax = min(f(G), b);                  % upper bounds
    VarMin = a;                             % lower bounds

    % Interpolate (only relevant if restrating with different grid)
    opts = interp1(linspace(a, b, length(opts)), opts, G, 'pchip');
    opts = max(min(opts, VarMax), VarMin);

    % Handle inactive grid points below threshold
    active_idx = true(1, nVar);             % preset boolean nVar
    if sm > a
        active_idx = (G > sm);        
        opts(~active_idx) = dyncon(~active_idx);        
    end
    
    % Cost function with constraint preservation
    CostFunction = @(x) Lossfunction(a, b, n, mandate(x, dyncon, active_idx), ...
                                     dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);

    % Initial candidate generation via supervised policy iteration
    % By default RESTART = 0 when using tsa
    if RESTART == 1
        for m = 1:100
            [current_cost, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);
            opts2 = mandate(transpose(Br), dyncon, active_idx);
            
            if ~isequal(opts2, opts)
                new_cost = CostFunction(opts2);
                if new_cost < current_cost
                    opts = opts2; 
                else
                    break;
                end
                if new_cost == 0
                    opts = opts2;
                    disp('Input policy is a Markov policy with Ego Loss = 0. Finishing.');
                    return;
                end
            else
                break;
            end
        end
    end

    %% CMA-ES parameters
    lambda = 8;                                                 % population size
    mu = round(lambda/2);                                       % parent number
    w = log(mu + 0.5) - log(1:mu);                              % selection weights
    w = w / sum(w);                                             % norma 
    mu_eff = 1 / sum(w.^2);                                     % effective parent number
    sigma0 = innerstep * (max(VarMax) - VarMin);                % initial step size
    cs = (mu_eff + 2) / (nVar + mu_eff + 5);                    % cumulation parameter
    ds = 1 + cs + 2*max(sqrt((mu_eff - 1)/(nVar+1)) - 1, 0);    % damping parameter
    ENN = sqrt(nVar)*(1 - 1/(4*nVar) + 1/(21*nVar^2));          % expected normal length
    cc = (4 + mu_eff / nVar) / (4 + nVar + 2*mu_eff / nVar);    % covariance cumulation
    c1 = 2 / ((nVar + 1.3)^2 + mu_eff);                         % rank-one update
    alpha_mu = 2;                                               % rank-mu factor
    cmu = min(1 - c1, alpha_mu*(mu_eff - 2 + 1/mu_eff) / ((nVar + 2)^2 + alpha_mu*mu_eff/2));
                                                                % rank-mu update
    hth = (1.4 + 2/(nVar+1)) * ENN;                             % stall threshold

    % Initialize CMA-ES state
    ps = cell(MaxIt, 1);                    % evolution path
    pc = cell(MaxIt, 1);                    % covariance path
    C = cell(MaxIt, 1);                     % covariance matrix
    sigma = cell(MaxIt, 1);                 % step size
    ps{1} = zeros(VarSize);
    pc{1} = zeros(VarSize);
    sigma{1} = sigma0;
    bias_factor = 0.5;                      % initial covariance bias
    initialCov = eye(nVar) + bias_factor*(ones(nVar) - eye(nVar));
    C{1} = initialCov;

    % Initialize population
    M(1).Position = mandate(opts, dyncon, active_idx);
    M(1).Step     = zeros(VarSize);
    M(1).Cost     = CostFunction(M(1).Position);
    BestSol       = M(1);                   % best solution
    BestCost      = zeros(MaxIt, 1);        % cost history
    stall_counter = 0;                      % short stall counter
    major_stall_counter = 0;                % long stall counter

    % Initialize waitbar
    hWaitbar = waitbar(0, 'Initiating Algorithm', 'Name', 'Solving (cancel to skip)', 'CreateCancelBtn','delete(gcbf)');
    set(hWaitbar, 'Units','Pixels','Position',[50 500 380 100]);

    %% Main CMA-ES loop
    for g = 1:MaxIt
        % Generate population with two-step best response
        pop = repmat(struct('Position',[],'Step',[],'Cost',[]), lambda, 1);
        for i = 1:lambda
            step_i = mvnrnd(zeros(VarSize), C{g});
            cand   = M(g).Position + sigma{g} * step_i;
            cand   = max(min(cand, VarMax), VarMin);
            cand   = mandate(cand, dyncon, active_idx);

            % Apply two-step best response transformation
            x = cand;
            [~, y] = Lossfunction(a, b, n, x, dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);
            y = max(min(transpose(y), VarMax), VarMin);
            y = mandate(y, dyncon, active_idx);
            [~, z] = Lossfunction(a, b, n, y, dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);
            z = max(min(transpose(z), VarMax), VarMin);
            z = mandate(z, dyncon, active_idx);

            val = CostFunction(z);
            pop(i).Step     = step_i;
            pop(i).Position = z;
            pop(i).Cost     = val;
            if val < BestSol.Cost && isreal(val) && val > 0
                BestSol = pop(i);
            end
        end

        % Sort population by cost
        [~, idxSo] = sort([pop.Cost]);
        pop        = pop(idxSo);
        BestCost(g) = BestSol.Cost;
        opts       = BestSol.Position;
        opts       = mandate(opts, dyncon, active_idx);

        % Check convergence and display progress
        if g > 1
            if BestCost(g)/nVar < tola
                disp(['Ego Loss equals zero within tolerance (' num2str(BestCost(g)/nVar) '). Finishing']);
                break;
            end
            if ~ishandle(hWaitbar)
                break;
            else
                waitbar((g-1)/MaxIt, hWaitbar, [' Ego Loss = ' num2str(BestCost(g-1)/nVar) ' (Iteration ' num2str(g-1) ' of ' num2str(MaxIt) ' )']);
            end
        end
        if ~isreal(BestCost(g))
            break;
        end
        if g == MaxIt
            close(hWaitbar);
            break;
        end

        % Update CMA-ES mean with two-step best response
        M(g+1).Step = zeros(VarSize);
        for jj = 1:mu
            M(g+1).Step = M(g+1).Step + w(jj)*pop(jj).Step;
        end

        x = M(g).Position + sigma{g}*M(g+1).Step;
        x = max(min(x, VarMax), VarMin);
        x = mandate(x, dyncon, active_idx);
        [~, y] = Lossfunction(a, b, n, x, dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);
        y = max(min(transpose(y), VarMax), VarMin);
        y = mandate(y, dyncon, active_idx);
        [~, z] = Lossfunction(a, b, n, y, dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);
        z = max(min(transpose(z), VarMax), VarMin);
        z = mandate(z, dyncon, active_idx);
        testVal = CostFunction(z);
        
        if testVal > 0 && isreal(testVal)
            M(g+1).Position = z;
            M(g+1).Cost     = testVal;
            if testVal < BestSol.Cost
                BestSol = M(g+1); 
            end
        else
            M(g+1).Position = M(g).Position;
            M(g+1).Cost     = CostFunction(M(g).Position);
        end

        % Adapt step size and covariance matrix
        ps{g+1} = (1-cs)*ps{g} + sqrt(cs*(2-cs)*mu_eff)*(M(g+1).Step / chol(C{g})');
        sigma{g+1} = sigma{g} * exp(cs/ds*(norm(ps{g+1})/ENN - 1))^0.3;
        hs = norm(ps{g+1}) / sqrt(1-(1-cs)^(2*(g+1))) < hth;
        delta_c = (1 - hs)*cc*(2 - cc);
        pc{g+1} = (1-cc)*pc{g} + hs*sqrt(cc*(2-cc)*mu_eff)*M(g+1).Step;

        C{g+1} = (1 - c1 - cmu)*C{g} + c1*(pc{g+1}'*pc{g+1} + delta_c*C{g});
        for jj = 1:mu
            C{g+1} = C{g+1} + cmu*w(jj)*(pop(jj).Step'*pop(jj).Step);
        end
        [VV, EE] = eig(C{g+1});
        if any(diag(EE) < 0)
            EE = max(EE, 0);
            C{g+1} = VV*EE/VV;
        end

        % Stall detection and restart
        if g > 1
            if BestCost(g) >= BestCost(g-1)
                stall_counter = stall_counter + 1;
                major_stall_counter = major_stall_counter + 1;
            else
                stall_counter = 0;
                major_stall_counter = 0;
            end
        end
        
        % Handle short-term stalls
        if stall_counter >= 5
            innerstep = innerstep/10;
            sigma{g+1} = innerstep*(max(VarMax) - VarMin);
            C{g+1} = initialCov;
            ps{g+1} = zeros(VarSize);
            pc{g+1} = zeros(VarSize);
            M(g+1).Position = BestSol.Position;
            stall_counter = 0;
        end
        
        % Handle long-term stalls
        if major_stall_counter >= 25
            innerstep = initial_innerstep;
            sigma{g+1} = innerstep*(max(VarMax) - VarMin);
            C{g+1} = initialCov;
            ps{g+1} = zeros(VarSize);
            pc{g+1} = zeros(VarSize);
            M(g+1).Position = BestSol.Position;
            major_stall_counter = 0;
        end
    end

    % Cleanup
    if exist('hWaitbar','var') && ishandle(hWaitbar)
        close(hWaitbar); 
    end
    opts = BestSol.Position;
end

function result = mandate(candidate, original, active_idx)
    % Preserve inactive entries from original policy
    result = candidate;
    result(~active_idx) = original(~active_idx);
end