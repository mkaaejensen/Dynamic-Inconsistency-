function [opts] = algorithm_call_sep(a, b, n, MaxIt, innerstep, opts, RESTART, dissaving, dis, P, u, tol, f_low, f_high, dsp, power_param, tolalg, G, IncomeG_low, IncomeG_high)
% Heavy-tailed separable CMA-ES algorithm with policy iteration
% Uses Student-t distribution for heavy-tailed sampling and diagonal covariance matrix
% Optimizes monotonic savings policies on asset-income grids

% Algorithm control parameters
wtc = 0;                    % Console output flag
WB = 1;                     % Waitbar display flag
policyinj = 1;              % Policy injection flag
gradinj = 1;                % Gradient injection flag

% Heavy-tailed distribution parameters
nu = 3;                     % Degrees of freedom for Student-t distribution

%% Pre-compute transformation bounds
upper_low = min(IncomeG_low, b) - a;
upper_high = min(IncomeG_high, b) - a;

%% Store initial loss
[initial_loss, ~] = Lossfunction(a, b, n, opts, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);
%% Policy iteration restart
%if RESTART == 1
 %   for m = 1:100
  %      [old_cost, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);
  %      opts2 = Br(:);
  %      [new_cost, ~] = Lossfunction(a, b, n, opts2, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);
  %      if new_cost < old_cost
 %           opts = opts2;
 %       else
 %           break;
 %       end
 %   end
%end

if RESTART == 2    % Single policy iteration step
    [~, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);
    opts = transpose(Br(:));
end
opts = opts(:);

%% Define cost function
CostFunction = @(z) Lossfunction(a, b, n, ...
    ztox(z, upper_low, upper_high, a, n), ...
    dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);

nVar = 2 * (n + 1);

%% Initialize CMA-ES parameters
z_init = xtoz(opts, upper_low, upper_high, a, n);

% CMA-ES configuration
lambda = 3;
mu = 2;
w = [0.6 0.4];
mu_eff = 1 / sum(w.^2);

% Strategy parameters
c_c = 4 / (nVar + 4);
c_1 = 0;                    % No rank-one update for diagonal scheme
alpha_mu = 2;
c_mu = min(1, alpha_mu * (mu_eff - 2 + 1/mu_eff) / ((nVar + 2)^2 + alpha_mu*mu_eff/2));

sigma0 = innerstep;
cs = (mu_eff + 2) / (nVar + mu_eff + 5);
ds = 1 + cs + 2*max(sqrt((mu_eff - 1)/(nVar+1)) - 1, 0);
ENN = sqrt(nVar)*(1 - 1/(4*nVar) + 1/(21*nVar^2));
cc = c_c;
c1 = c_1;
cmu = c_mu;
hth = (1.4 + 2/(nVar+1))*ENN;

% Separable CMA-ES state variables
ps_curr = zeros(1, nVar);
D_curr = ones(nVar, 1);     % Diagonal elements of covariance matrix
sigma_curr = sigma0;

% State tracking
M_pos = z_init(:);
M_step = zeros(nVar, 1);
M_cost = CostFunction(M_pos);

BestSol_pos = M_pos;
BestSol_cost = M_cost;
BestCost_prev = Inf;
stall_counter = 0;

% Pre-allocate arrays
pop_costs = zeros(lambda, 1);
pop_steps = zeros(lambda, nVar);
pop_positions = zeros(lambda, nVar);
S = zeros(nVar, lambda);

if WB == 1
    hWaitbar = waitbar(0, 'Initiating', ...
        'Name', 'Solving (cancel to skip)', 'CreateCancelBtn','delete(gcbf)');
    set(hWaitbar, 'Units','Pixels','Position',[50 500 380 100]);
end

%% Main CMA-ES loop
brcount = 0;
grad_alpha = 1.0;           % Gradient step size parameter

for g = 1:MaxIt
    % Generate candidates using Student-t distribution
    Z = trnd(nu, lambda, nVar);
    stepsAll = bsxfun(@times, Z, sqrt(D_curr)');
    candMat = repmat(M_pos', lambda, 1) + sigma_curr * stepsAll;
    
    brcount = brcount + 1;
    
    % Handle special candidates
    max_special = 2;
    special_costs = zeros(max_special, 1);
    special_steps = zeros(max_special, nVar);
    special_positions = zeros(max_special, nVar);
    num_special = 0;
    
    % Policy injection
    if brcount == 2 && policyinj == 1
        best_policy = ztox(BestSol_pos, upper_low, upper_high, a, n);
        
        [~, Br] = Lossfunction(a, b, n, best_policy, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);
        
        policy_cand = xtoz(Br, upper_low, upper_high, a, n);
        policy_step = zeros(nVar, 1);
        policy_cost = CostFunction(policy_cand);
        
        % Report improvement
        if policy_cost < BestSol_cost && isreal(policy_cost) && policy_cost > 0
            improvement = BestSol_cost - policy_cost;
            if wtc == 1
                fprintf('Policy injection successful: improvement = %.4e\n', improvement);
            end
        end
        
        num_special = num_special + 1;
        special_costs(num_special) = policy_cost;
        special_steps(num_special, :) = policy_step';
        special_positions(num_special, :) = policy_cand';
    end
    
    % Gradient injection
    if brcount == 3 && gradinj == 1
        grad_estimate = estimateGradient(M_pos, M_cost, pop_positions, pop_costs, sigma_curr, lambda);
        
        if norm(grad_estimate) > 1e-12
            grad_step = -grad_estimate / norm(grad_estimate);
            grad_cand = M_pos + grad_alpha * sigma_curr * grad_step;
            grad_step_store = grad_alpha * grad_step;
            grad_cost = CostFunction(grad_cand);
            
            % Report improvement
            if grad_cost < BestSol_cost && isreal(grad_cost) && grad_cost > 0
                improvement = BestSol_cost - grad_cost;
                if wtc == 1
                    fprintf('Gradient injection successful: improvement = %.4e (grad_norm = %.3e)\n', ...
                            improvement, norm(grad_estimate));
                end
            end
        else
            % Use standard candidate as fallback
            grad_step_store = stepsAll(lambda - num_special, :)';
            grad_cand = candMat(lambda - num_special, :)';
            grad_cost = CostFunction(grad_cand);
        end
        
        num_special = num_special + 1;
        special_costs(num_special) = grad_cost;
        special_steps(num_special, :) = grad_step_store';
        special_positions(num_special, :) = grad_cand';
        
        brcount = 0;
    end
    
    % Determine number of standard candidates
    lambda_standard = lambda - num_special;
    
    % Pre-allocate arrays for standard candidates
    standard_costs = zeros(lambda_standard, 1);
    standard_steps = zeros(lambda_standard, nVar);
    standard_positions = zeros(lambda_standard, nVar);
    
    % Generate standard candidates in parallel
    for i = 1:lambda_standard
        step_i = stepsAll(i, :)';
        cand_i = candMat(i, :)';
        cost_i = CostFunction(cand_i);
        
        standard_costs(i) = cost_i;
        standard_steps(i, :) = step_i';
        standard_positions(i, :) = cand_i';
    end
    
    % Combine standard and special candidates
    if num_special > 0
        pop_costs = [standard_costs; special_costs(1:num_special)];
        pop_steps = [standard_steps; special_steps(1:num_special, :)];
        pop_positions = [standard_positions; special_positions(1:num_special, :)];
    else
        pop_costs = standard_costs;
        pop_steps = standard_steps;
        pop_positions = standard_positions;
    end
    
    % Update best solution
    for i = 1:lambda
        if pop_costs(i) < BestSol_cost && isreal(pop_costs(i)) && (pop_costs(i) > 0)
            BestSol_pos = pop_positions(i, :)';
            BestSol_cost = pop_costs(i);
        end
    end
    
    % Sort population by cost
    [pop_costs, idxSo] = sort(pop_costs);
    pop_steps = pop_steps(idxSo, :);
    pop_positions = pop_positions(idxSo, :);
    
    BestCost_curr = BestSol_cost;
    opts = ztox(BestSol_pos, upper_low, upper_high, a, n);
    
    % Store sorted steps
    S = pop_steps';
    
    if (BestCost_curr/nVar) < tolalg
        disp(['Zero loss up to tolerance (' sprintf('%.4e', BestCost_curr/nVar) '). Finishing']);
        if WB==1, close(hWaitbar); end
        opts = roundToA(opts, a);
        clear cached_W_out cached_initialized VF_low_interp VF_high_interp tFine_power numGrid columnIndices_template vf_low_interp_vfi vf_high_interp_vfi last_grid_size;
        break;
    end
    
    if g == MaxIt
        disp(['Initial Loss = ' sprintf('%.4e', initial_loss/nVar) ', Final Loss = ' sprintf('%.4e', BestCost_curr/nVar)]);
        if WB==1, close(hWaitbar); end
        opts = roundToA(opts, a);
        clear cached_W_out cached_initialized VF_low_interp VF_high_interp tFine_power numGrid columnIndices_template vf_low_interp_vfi vf_high_interp_vfi last_grid_size;
        break;
    end
    
    % Compute weighted average step
    M_step_next = (w * pop_steps(1:mu, :))';
    
    M_pos_next = M_pos + sigma_curr * M_step_next;
    M_cost_next = CostFunction(M_pos_next);
    
    if M_cost_next > 0 && isreal(M_cost_next)
        M_pos = M_pos_next;
        M_cost = M_cost_next;
        M_step = M_step_next;
        if M_cost_next < BestSol_cost
            BestSol_pos = M_pos_next;
            BestSol_cost = M_cost_next;
        end
    else
        M_step = zeros(nVar, 1);
    end
    
    % Update current best cost
    BestCost_curr = BestSol_cost;
    
    % Termination checks
    if g > 1
        if WB==1 && ~ishandle(hWaitbar)
            disp('Canceled by user');
            disp(['Initial Loss = ' sprintf('%.4e', initial_loss/nVar) ', Final Loss = ' sprintf('%.4e', BestCost_curr/nVar)]);
            g = MaxIt;
            opts = roundToA(opts, a);
            clear cached_W_out cached_initialized VF_low_interp VF_high_interp tFine_power numGrid columnIndices_template vf_low_interp_vfi vf_high_interp_vfi last_grid_size;
            break;
        elseif WB==1
            waitbar((g-1)/MaxIt, hWaitbar, ...
                ['Ego Loss = ' sprintf('%.4e', BestCost_prev/nVar) ...
                 ' (Iteration ' num2str(g-1) ' of ' num2str(MaxIt) ')']);
        end
    end
    
    % Evolution path and step size update
    M_step_scaled = M_step ./ sqrt(D_curr);
    ps_next = (1-cs)*ps_curr + sqrt(cs*(2-cs)*mu_eff) * M_step_scaled';
    ps_next = ps_next(:)';
    sigma_next = sigma_curr * exp(cs/ds*(norm(ps_next)/ENN - 1));
    
    % Diagonal covariance update
    weighted_squares = bsxfun(@times, S(:,1:mu).^2, w);
    D_next = (1 - c_mu) * D_curr + c_mu * sum(weighted_squares, 2);
    
    % Ensure D stays positive and bounded
    D_next = max(1e-10, min(1e10, D_next));
    
    % Stall detection and restart
    if g > 1
        if BestCost_curr >= BestCost_prev
            stall_counter = stall_counter + 1;
        else
            stall_counter = 0;
        end
    end
    
    if stall_counter >= 50
        innerstep = 0.8*innerstep;
        sigma_next = innerstep;
        D_next = ones(nVar, 1);
        ps_next = zeros(1, nVar);
        M_pos = BestSol_pos;
        stall_counter = 0;
        if wtc==1
            disp(['Reducing innerstep to ' num2str(innerstep) ' and restarting']);
        end
    end
    
    % Update state variables
    ps_curr = ps_next;
    D_curr = D_next;
    sigma_curr = sigma_next;
    BestCost_prev = BestCost_curr;
end

if WB==1 && ishandle(hWaitbar), close(hWaitbar); end

opts = ztox(BestSol_pos, upper_low, upper_high, a, n);
opts = roundToA(opts, a);

end

 
%% Parameter Space Transformation and Gradient Estimation Functions
 

function x = ztox(z, upper_low, upper_high, a, n)
    % Transform parameters to monotonic savings policies
    c_low = cumsum(z(1:n+1));
    c_high = cumsum(z(n+2:2*(n+1)));
    
    f_low = 1 ./ (1 + exp(-c_low));
    f_high = 1 ./ (1 + exp(-c_high));
    
    x = [a + upper_low .* f_low;
         a + upper_high .* f_high];
end

function z = xtoz(x, upper_low, upper_high, a, n)
    % Inverse transform from policy to parameter space
    eps = 1e-15;
    
    % Process both states
    x_low = x(1:n+1);
    x_high = x(n+2:2*(n+1));
    
    % Low state
    f_low = (x_low - a) ./ upper_low;
    f_low = max(eps, min(1-eps, f_low));
    c_low = log(f_low ./ (1 - f_low));
    
    % High state
    f_high = (x_high - a) ./ upper_high;
    f_high = max(eps, min(1-eps, f_high));
    c_high = log(f_high ./ (1 - f_high));
    
    % Recover increments
    z = [[c_low(1); diff(c_low)];
         [c_high(1); diff(c_high)]];
end

function vals = roundToA(vals, a)
    % Correct floating point imprecision at borrowing limit
    closeToA = abs(vals - a) < 1e-8;
    vals(closeToA) = a;
end

function grad_estimate = estimateGradient(center_pos, center_cost, positions, costs, sigma, lambda)
    % Estimate pseudo-gradient from function evaluations
    % Uses finite difference approximation with improving samples
    
    nVar = length(center_pos);
    grad_estimate = zeros(nVar, 1);
    
    for i = 1:lambda
        if costs(i) < center_cost && isreal(costs(i)) && costs(i) > 0
            % Use only improving directions
            step_direction = (positions(i, :)' - center_pos) / sigma;
            improvement = center_cost - costs(i);
            
            % Accumulate weighted gradient estimate
            grad_estimate = grad_estimate - improvement * step_direction;
        end
    end
    
    % Normalize by number of improving samples
    n_improving = sum(costs < center_cost & isreal(costs') & costs > 0);
    if n_improving > 0
        grad_estimate = grad_estimate / n_improving;
    end
end