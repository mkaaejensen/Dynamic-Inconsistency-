function [opts] = algorithm_call_tsa(a, b, n, MaxIt, innerstep, opts, RESTART, dissaving, dis, p, u, tol, f, offset,dsp)
% algorithm_call implements the CMA-ES algorithm with policy iteration.
% Inputs:
%   a, b        - lower and upper bounds of the grid
%   n           - grid resolution (number of subintervals)
%   MaxIt       - maximum number of iterations
%   innerstep   - initial step size for CMA-ES
%   opts        - initial policy vector
%   RESTART     - determines type of policy iteration
%   speedy      - speed optimization flag
%   dis         - discount sequence (as a function handle)
%   p           - probability parameter
%   u           - utility function handle
%   tol         - tolerance parameter
%   f           - income function handle
%   offset      - offset function handle
%
% Output:
%   opts        - updated policy vector after algorithm convergence

    %% Optional Parameters
    WB = 1;           
    pelts = 0;        
    policyboost = 1;  
    pols = 1;         
    innerpolicyboost = 1;
    FLB = a;  % Feasible lower bound



    %% Initialization
    eta = ones(1, n+1);
    omega = 1;
    % Use Lossfunction directly as the cost function.
    CostFunction = @(x) Lossfunction(a, b, n, x, dissaving, dis, u, p, 0, tol, f, offset,dsp);
    
    nVar = n + 1;
    VarSize = [1, nVar];

    % Compute grid points G(j)
    G = linspace(a, b, nVar);
    
    % Compute variable-specific upper bounds (clamped by b)
    VarMax = min(f(G), b);
    VarMin = a;
    FLB = VarMin;

    % Possibly re-interpolate opts to match dimension
    opts = interp1(linspace(a, b, length(opts)), opts, G, 'pchip');
    opts = max(min(opts, VarMax), VarMin);

    % Single policy iteration if RESTART == 2
    if RESTART == 2
        [~, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 0, tol, f, offset,dsp);
        opts = max(min(transpose(Br), VarMax), VarMin);
    end

    % If RESTART == 1, do iterative policy improvement
    if RESTART == 1
        for m = 1:100
            [~, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 0, tol, f, offset,dsp);
            opts2 = max(min(transpose(Br), VarMax), VarMin);
            if CostFunction(opts2) < CostFunction(opts)
                opts = opts2; 
            else
                break;
            end
            if CostFunction(opts2) == 0
                opts = opts2;
                disp('Input policy is a Markov policy with Ego Loss = 0. Finishing.');
                return;
            end
        end
    end

    %% CMA-ES Specific Settings
    lambda = 8;%(4 + round(3*log(nVar))) * 2; 
    mu = round(lambda/2);
    w = log(mu + 0.5) - log(1:mu);
    w = w / sum(w);
    mu_eff = 1 / sum(w.^2);
    
    sigma0 = innerstep*(max(VarMax) - VarMin);
    cs = (mu_eff + 2) / (nVar + mu_eff + 5);
    ds = 1 + cs + 2*max(sqrt((mu_eff - 1)/(nVar+1)) - 1, 0);
    ENN = sqrt(nVar)*(1 - 1/(4*nVar) + 1/(21*nVar^2));
    cc = (4 + mu_eff / nVar) / (4 + nVar + 2*mu_eff / nVar);
    c1 = 2 / ((nVar + 1.3)^2 + mu_eff);
    alpha_mu = 2;
    cmu = min(1 - c1, alpha_mu*(mu_eff - 2 + 1/mu_eff) / ((nVar + 2)^2 + alpha_mu*mu_eff/2));
    hth = (1.4 + 2/(nVar+1))*ENN;

    ps = cell(MaxIt, 1);
    pc = cell(MaxIt, 1);
    C = cell(MaxIt, 1);
    sigma = cell(MaxIt, 1);
    ps{1} = zeros(VarSize);
    pc{1} = zeros(VarSize);
    sigma{1} = sigma0;
    % Bias the covariance matrix in the direction of increasing solutions.
    bias_factor = 0.5;
    initialCov = eye(nVar) + bias_factor*(ones(nVar) - eye(nVar));
    C{1} = initialCov;

    % For CMA-ES: prepare an empty individual
    empty_individual.Position = [];
    empty_individual.Step = [];
    empty_individual.Cost = [];
    
    M(1).Position = opts;
    M(1).Step = zeros(VarSize);
    M(1).Cost = CostFunction(M(1).Position);
    
    if M(1).Cost < 0
        disp('Ego Loss negative - initial solution stems from finer ORACLE. Stopping.');
        MaxIt = 0;
    end

    BestSol = M(1);
    BestCost = zeros(MaxIt, 1);
    count = 0;
    setflag = 0;
    stall_counter = 0;
    
    if WB == 1
        hWaitbar = waitbar(0, 'Initiating Algorithm', 'Name', 'Solving (cancel to skip)', ...
                           'CreateCancelBtn','delete(gcbf)');
        set(hWaitbar, 'Units','Pixels','Position',[50 500 380 100]);
    end

    for g = 1:MaxIt
        count = count + 1;
        % Generate population
        pop = repmat(empty_individual, lambda, 1);
        for i = 1:lambda
            step_i = mvnrnd(zeros(VarSize), C{g});
            cand = M(g).Position + sigma{g} * step_i;
            cand = max(min(cand, VarMax), VarMin);
            % --- Two-step lookahead begins ---
            % (i) x is the candidate from CMA-ES:
            x = cand;
            % (ii) Compute best response y to x:
            [~, y] = Lossfunction(a, b, n, x, dissaving, dis, u, p, 0, tol, f, offset,dsp);
            y = max(min(transpose(y), VarMax), VarMin);
            % (iii) Compute best response z to y:
            [~, z] = Lossfunction(a, b, n, y, dissaving, dis, u, p, 0, tol, f, offset,dsp);
            z = max(min(transpose(z), VarMax), VarMin);
            % (iv) Evaluate the loss at z:
            val = CostFunction(z);
            % --- Two-step lookahead ends ---
            pop(i).Step = step_i;  % (Step is as originally generated.)
            pop(i).Position = z;   % The effective candidate is z.
            pop(i).Cost = val;
            if val < BestSol.Cost && isreal(val) && (val > 0)
                BestSol = pop(i);
            end
        end
        
        % Sort population by cost
        Costs = [pop.Cost];
        [Costs, idxSo] = sort(Costs);
        pop = pop(idxSo);
        
        BestCost(g) = BestSol.Cost;
        opts = BestSol.Position;
        
        if g > 1
            if setflag == 1
                if pelts == 1
                    disp(['Iteration ' num2str(g-1) ' [' num2str(n) ' gp]: Ego Loss = ' num2str(BestCost(g-1)/nVar) ' [reweighting]']);
                end
                setflag = 0;
            else
                if pelts == 1
                    disp(['Iteration ' num2str(g-1) ' [' num2str(n) ' gp]: Ego Loss = ' num2str(BestCost(g-1)/nVar)']);
                end
            end
            
            if BestCost(g)/nVar<1e-17
                g = MaxIt;
            end
           
            if WB==1 && ~ishandle(hWaitbar)
                g = MaxIt;                              
                break;
            elseif WB==1
                waitbar((g-1)/MaxIt, hWaitbar, [' Ego Loss = ' num2str(BestCost(g-1)/nVar) ' (Iteration ' num2str(g-1) ' of ' num2str(MaxIt) ' )']);
            end
        end
        
        if BestCost(g) < 0
            disp('Divergence: Cost < 0. Check domain or step length.');
            break;
        end
        if ~isreal(BestCost(g))
            disp('Divergence: Cost is imaginary. Domain issue?');
            break;
        end
        
        if g == MaxIt
            disp(['Iteration ' num2str(g) ' [final] [' num2str(n) ' gp]: Ego Loss = ' num2str(BestCost(g)/nVar)']);
            if WB==1, close(hWaitbar); end
            break;
        end
        
        % Update mean using selected population steps
        M(g+1).Step = zeros(VarSize);
        for jj = 1:mu
            M(g+1).Step = M(g+1).Step + w(jj)*pop(jj).Step;
        end
        
        % --- Two-step lookahead for the weighted candidate ---
        x = M(g).Position + sigma{g}*M(g+1).Step;
        x = max(min(x, VarMax), VarMin);
        [~, y] = Lossfunction(a, b, n, x, dissaving, dis, u, p, 0, tol, f, offset,dsp);
        y = max(min(transpose(y), VarMax), VarMin);
        [~, z] = Lossfunction(a, b, n, y, dissaving, dis, u, p, 0, tol, f, offset,dsp);
        z = max(min(transpose(z), VarMax), VarMin);
        testVal = CostFunction(z);
        
        if ~isreal(testVal)
            disp('Imaginary test candidate encountered.');
        end
        if testVal > 0 && isreal(testVal)
            M(g+1).Position = z;
            M(g+1).Cost = testVal;
            if testVal < BestSol.Cost && testVal > 0 && isreal(testVal)
                BestSol = M(g+1);
            end
        else
            M(g+1).Position = M(g).Position;
            M(g+1).Cost = CostFunction(M(g).Position);
        end
        
        % Update step size and covariance matrix
        ps{g+1} = (1-cs)*ps{g} + sqrt(cs*(2-cs)*mu_eff)*(M(g+1).Step / chol(C{g})');
        sigma{g+1} = sigma{g} * exp(cs/ds*(norm(ps{g+1})/ENN - 1))^0.3;
        
        if norm(ps{g+1}) / sqrt(1-(1-cs)^(2*(g+1))) < hth
            hs = 1;
        else
            hs = 0;
        end
        delta_c = (1 - hs)*cc*(2 - cc);
        pc{g+1} = (1-cc)*pc{g} + hs*sqrt(cc*(2-cc)*mu_eff)*M(g+1).Step;
        
        C{g+1} = (1 - c1 - cmu)*C{g} + c1*(pc{g+1}'*pc{g+1} + delta_c*C{g});
        for jj = 1:mu
            C{g+1} = C{g+1} + cmu*w(jj)*(pop(jj).Step'*pop(jj).Step);
        end
        
        % Ensure C is positive definite
        [VV, EE] = eig(C{g+1});
        if any(diag(EE) < 0)
            EE = max(EE, 0);
            C{g+1} = VV*EE/VV;
        end
        
        %  
        if g > 1
            if BestCost(g) >= BestCost(g-1)
                stall_counter = stall_counter + 1;
            else
                stall_counter = 0;
            end
        end
        if stall_counter >= 5
            %disp('Stall detected, restarting with reduced step length');
            innerstep = innerstep/10;
            sigma{g+1} = innerstep*(max(VarMax) - VarMin);
            C{g+1} = initialCov;
            ps{g+1} = zeros(VarSize);
            pc{g+1} = zeros(VarSize);
            M(g+1).Position = BestSol.Position;
            stall_counter = 0;
            %disp(['Adjusting innerstep to ' num2str(innerstep) ' and restarting CMA-ES'])
        end
    end
    
    EgoLoss = BestCost(g)/nVar;
    
end