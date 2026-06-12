function [loss, br] = Lossfunction(a, b, n, x, dissaving, dis, u, P, rep, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high)
% Value function iteration with grid search optimization for two-state Markov process
% Implements caching and interpolation for computational efficiency
% Now accepts VFI grid G and pre-computed income values
x = x(:);

% Store extended policies when rep=3 (extended grid mode)
persistent extended_opts_low extended_opts_high extended_grid extended_bounds

if rep == 3
    extended_opts_low = x(1:(n+1));
    extended_opts_high = x((n+2):end);
    extended_grid = G;
    extended_bounds = [a, b];
end

% Persistent storage for computational efficiency
persistent cached_W_out cached_initialized
persistent VF_low_interp VF_high_interp
persistent DS_grid numGrid columnIndices_template
persistent cached_grid_hash cached_n cached_dsp

% Initialize computational cache on first call or when grid parameters change
current_grid_hash = G(1) + G(end) + length(G);  % Simple grid fingerprint
cache_invalid = isempty(cached_initialized) || numGrid ~= (n+1)*dsp || ...
                cached_grid_hash ~= current_grid_hash || cached_n ~= n || cached_dsp ~= dsp;

if cache_invalid
    cached_W_out = [];
    cached_initialized = true;
    cached_grid_hash = current_grid_hash;
    cached_n = n;
    cached_dsp = dsp;
    
    % Pre-compute grid search components
    numGrid = (n+1)*dsp;
    DS_grid = linspace(0, 1, numGrid);
    
    % Pre-compute indexing template for tie-breaking in optimization
    columnIndices_template = repmat(1:numGrid, n+1, 1);
    
    % Clear interpolants to force rebuild
    VF_low_interp = [];
    VF_high_interp = [];
end

% Extract discount parameters from input vector
delta = dis(2)/dis(1);
beta = dis(1)/delta;

% Parse state-dependent policy variables
x_low = x(1:n+1);
x_high = x(n+2:end);
s_mat = [x_low, x_high];

% Initialize value function when cache is empty or grid size changed
if isempty(cached_W_out) || size(cached_W_out, 1) ~= length(G)
    % Myopic consumption-based value function initialization
    V_low_init = u(f_low(0)) * ones(size(G));
    V_high_init = u(f_high(0)) * ones(size(G));
    
    init_W_out = [V_low_init, V_high_init];
else
    init_W_out = cached_W_out;
end

% Create extended policy interpolants if available (for boundary conditions)
% Only use for main grid calls, not during extended optimization (rep=3)
extended_low_interp = [];
extended_high_interp = [];
if rep ~= 3 && ~isempty(extended_opts_low) && ~isempty(extended_grid)
    extended_low_interp = griddedInterpolant(extended_grid, extended_opts_low, 'linear', 'linear');
    extended_high_interp = griddedInterpolant(extended_grid, extended_opts_high, 'linear', 'linear');
end

% Solve for optimal value function via iteration
[VFt, ~] = valf2(a, b, n, u, P, delta, tol, s_mat, f_low, f_high, G, IncomeG_low, IncomeG_high, init_W_out, extended_low_interp, extended_high_interp);

% Update cache
cached_W_out = VFt;

% Construct or update interpolation objects for value functions
if isempty(VF_low_interp) || length(VF_low_interp.GridVectors{1}) ~= length(G)
    VF_low_interp = griddedInterpolant(G, VFt(:,1), 'linear', 'linear');
    VF_high_interp = griddedInterpolant(G, VFt(:,2), 'linear', 'linear');
else
    VF_low_interp.Values = VFt(:,1);
    VF_high_interp.Values = VFt(:,2);
end
 
% Minimum consumption bound
min_consumption = 0.1*(f_low(a) - a); % setting bound at 10% of worst possible consumption (allows some floating point imprecision)

% Construct initial savings grids via broadcasting
% Matrix dimensions: (n+1) x numGrid
savingGrid_high = a + DS_grid .* (IncomeG_high - a);
savingGrid_low = a + DS_grid .* (IncomeG_low - a);

% Apply hard constraints for high income state
% Cannot save more than income less minimum consumption
max_saving_high = IncomeG_high - min_consumption;
savingGrid_constr_high = min(savingGrid_high, max_saving_high);

% Apply constraints for low income state
if dissaving == 1
    % When dissaving=1: cannot save more than current assets OR income minus min consumption
    max_saving_low = min(G, IncomeG_low - min_consumption);
    savingGrid_constr_low = min(savingGrid_low, max_saving_low);
else
    % When dissaving=0: only income constraint applies
    max_saving_low = IncomeG_low - min_consumption;
    savingGrid_constr_low = min(savingGrid_low, max_saving_low);
end

% Compute consumption from budget constraints
% Now guaranteed to be >= min_consumption due to hard constraints above
cons_low = IncomeG_low - savingGrid_constr_low;
cons_high = IncomeG_high - savingGrid_constr_high;

% Evaluate objective function for low-income state
cont_val_low = beta*delta*(P(1,1)*VF_low_interp(savingGrid_constr_low(:)) + P(1,2)*VF_high_interp(savingGrid_constr_low(:)));
obj_vals_low = reshape(u(cons_low(:)) + cont_val_low, n+1, numGrid);

% Evaluate objective function for high-income state
cont_val_high = beta*delta*(P(2,1)*VF_low_interp(savingGrid_constr_high(:)) + P(2,2)*VF_high_interp(savingGrid_constr_high(:)));
obj_vals_high = reshape(u(cons_high(:)) + cont_val_high, n+1, numGrid);

% Find optimal decisions via vectorized maximum search
[payoff_low, maxIdx_low] = max(obj_vals_low, [], 2);
[payoff_high, maxIdx_high] = max(obj_vals_high, [], 2);

% Implement tie-breaking rule using numerical tolerance
mask_low = abs(obj_vals_low - payoff_low) < 1e-12;
tieIndices_low = max(mask_low .* columnIndices_template, [], 2);

mask_high = abs(obj_vals_high - payoff_high) < 1e-12;
tieIndices_high = max(mask_high .* columnIndices_template, [], 2);

% Validate computed indices (making absolutely certain that there are no
% numerical issues with tie-braking rule)
validIndices_low = tieIndices_low > 0 & tieIndices_low <= numGrid;
validIndices_high = tieIndices_high > 0 & tieIndices_high <= numGrid;
tieIndices_low(~validIndices_low) = maxIdx_low(~validIndices_low);
tieIndices_high(~validIndices_high) = maxIdx_high(~validIndices_high);
payoff_low(~validIndices_low) = 1e10; % Assign penalty values for computational errors
payoff_high(~validIndices_high) = 1e10; % Assign penalty values for computational errors

% Extract optimal policies from indices 
policy_low = savingGrid_constr_low(sub2ind([n+1, numGrid], (1:n+1)', tieIndices_low));
policy_high = savingGrid_constr_high(sub2ind([n+1, numGrid], (1:n+1)', tieIndices_high));

% Compute policy loss functions for both states
% Low-income state loss calculation
if any(IncomeG_low - x_low <= 0)
   warning('Negative consumption detected in low-income state: min = %.2e', min(IncomeG_low - x_low));
end
first_mover_payoff_low = u(IncomeG_low - x_low) + beta*delta*(P(1,1)*VF_low_interp(x_low) + P(1,2)*VF_high_interp(x_low));
state_loss_low = payoff_low - first_mover_payoff_low;

% High-income state loss calculation
if any(IncomeG_high - x_high <= 0)
   warning('Negative consumption detected in high-income state: min = %.2e', min(IncomeG_high - x_high));
end
first_mover_payoff_high = u(IncomeG_high - x_high) + beta*delta*(P(2,1)*VF_low_interp(x_high) + P(2,2)*VF_high_interp(x_high));
state_loss_high = payoff_high - first_mover_payoff_high;

% Handle numerical edge cases in vectorized manner
nonRealMask_low = ~isreal(policy_low);
state_loss_low(nonRealMask_low) = 1e8;

negativeMask_low = state_loss_low < 0;
policy_low(negativeMask_low) = x_low(negativeMask_low);
state_loss_low(negativeMask_low) = 0;

nonRealMask_high = ~isreal(policy_high);
state_loss_high(nonRealMask_high) = 1e8;

negativeMask_high = state_loss_high < 0;
policy_high(negativeMask_high) = x_high(negativeMask_high);
state_loss_high(negativeMask_high) = 0;

% Aggregate loss function and construct best response policy
loss = sum(state_loss_low) + sum(state_loss_high);
br = [policy_low; policy_high];

% Return value function interpolant when requested
if rep == 2
    br = @(z) [VF_low_interp(z), VF_high_interp(z)];
end
end

function [W_out, grid] = valf2(a, b, n, u, P, delta, tol, s, f_low, f_high, G, IncomeG_low, IncomeG_high, init_W, extended_low_interp, extended_high_interp)
% Value function iteration 

% Iteration precision parameters
% tol : value-function iteration tolerance (set in main script)
iterlimit = 200; % maximum number of interations (can be said low because of caching)

% Persistent interpolation objects
persistent vf_low_interp_vfi vf_high_interp_vfi last_grid_hash

% Create or update interpolation objects based on grid changes
current_grid_hash = G(1) + G(end) + length(G);
if isempty(vf_low_interp_vfi) || isempty(last_grid_hash) || last_grid_hash ~= current_grid_hash
    vf_low_interp_vfi = griddedInterpolant(G, init_W(:,1), 'linear', 'linear');
    vf_high_interp_vfi = griddedInterpolant(G, init_W(:,2), 'linear', 'linear');
    last_grid_hash = current_grid_hash;
end

% High-precision mode with extended iteration limit
if tol == 1e-20
    iterlimit = 10000;    
end

% Initialize value function iteration
W_old = init_W;

iter = 0;
dist = Inf;

% Value function iteration loop
while dist > tol && iter < iterlimit
    s_low_policy = s(:,1);
    s_high_policy = s(:,2);
    
    cons_low = max(IncomeG_low - s_low_policy, 1e-14);
    cons_high = max(IncomeG_high - s_high_policy, 1e-14);
    
    % Update interpolant values  
    vf_low_interp_vfi.Values = W_old(:,1);
    vf_high_interp_vfi.Values = W_old(:,2);
    
    % Get continuation values - use extended policies for out-of-bounds savings
    v_low_low = vf_low_interp_vfi(s_low_policy);
    v_low_high = vf_high_interp_vfi(s_low_policy);
    v_high_low = vf_low_interp_vfi(s_high_policy);
    v_high_high = vf_high_interp_vfi(s_high_policy);
    
    % Use extended policies for sophisticated extrapolation when savings > b
    if ~isempty(extended_low_interp) && ~isempty(extended_high_interp)
        % Check for out-of-bounds savings
        out_low = s_low_policy > b;
        out_high = s_high_policy > b;
        
        if any(out_low)
            % For savings > b in low state, estimate value using extended policy
            s_out = s_low_policy(out_low);
            ext_savings_low = extended_low_interp(s_out);
            ext_savings_high = extended_high_interp(s_out);
            
            % Consumption if agent has assets s_out and follows extended policy
            cons_low_ext = f_low(s_out) - ext_savings_low;
            cons_high_ext = f_high(s_out) - ext_savings_high;
            
            % Myopic value estimate (sophisticated boundary condition)
            v_low_low(out_low) = u(max(cons_low_ext, 1e-14));
            v_low_high(out_low) = u(max(cons_high_ext, 1e-14));
        end
        
        if any(out_high)
            % For savings > b in high state, estimate value using extended policy
            s_out = s_high_policy(out_high);
            ext_savings_low = extended_low_interp(s_out);
            ext_savings_high = extended_high_interp(s_out);
            
            % Consumption if agent has assets s_out and follows extended policy
            cons_low_ext = f_low(s_out) - ext_savings_low;
            cons_high_ext = f_high(s_out) - ext_savings_high;
            
            % Myopic value estimate (sophisticated boundary condition)
            v_high_low(out_high) = u(max(cons_low_ext, 1e-14));
            v_high_high(out_high) = u(max(cons_high_ext, 1e-14));
        end
    end
    
    W_new_low = u(cons_low) + delta*(P(1,1)*v_low_low + P(1,2)*v_low_high);
    W_new_high = u(cons_high) + delta*(P(2,1)*v_high_low + P(2,2)*v_high_high);
    
    W_new = [W_new_low, W_new_high];
    
    dist = max(max(abs(W_new - W_old)));
    W_old = W_new;
    iter = iter + 1;
end

W_out = W_old;
grid = G;
end