function [SS, x, idr] = algorithm_ssr(a, b, n, u, dis, cct, opts, dissaving, f, offset, p, tol, con, branch)
%algorithm_renewal(a, b, n, u, dis, cct, opts, speedy, f, offset, p, tol, con, branch)
% Extended version that handles both deterministic and stochastic cases
%
% Parameters:
% a, b - lower and upper bounds for asset grid
% n - number of grid points
% u - utility function
% dis - discount factors array
% cct - =1 prints iteration number to console
% opts - initial policy function guess
% dissaving =1 searches only up to 45 degree line (speed-boost but may reduce precision)
% f - income function
% offset - maps low income into high income
% p - probability of low shock (0 <= p <= 1)
% tol - tolerance for value function iteration (default: 1e-6)
% con - continuation policy option: 
%     0 = use constant policy t for future periods (default)
%     1 = use opts for future periods
% branch - method selection parameter:
%     0 = always use value function iteration regardless of p (default)
%     1 = use summation for deterministic case (p=0 or p=1), VFI otherwise

 

% Determine delta and beta from dis
delta = dis(2)/dis(1);
beta = dis(1)/delta;

% Initialize variables
x = opts; % initial policy function
t = 1; count = 0; ELB = a; nss = 0; SS = -1; % starting up


% Pre-allocate for performance
idr = zeros(n+1, 1);
rlb = zeros(n+1, 1);
G_values = a + (0:n) * ((b-a)/n);
income_values = f(G_values);

% Upper bounds depending on dissaving
if dissaving > 0
    ub_values = G_values; % search only to 45 degree line
else
    ub_values = min(income_values, b); % full search
end

% Determine which method to use for value function computation
is_deterministic = (p == 0 || p == 1);
use_summation = (branch == 1 && is_deterministic);

% Main loop over asset grid
for j = 0:1:n    
    if count == 100 && cct == 1 && n < 1000
        display(['Iteration ' num2str(j) ' out of ' num2str(n)]);
        fprintf(1, '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b', n);
        count = 0;
    elseif count == 100 && cct == 1 && n >= 1000 && n < 10000 && j < 1000
        display(['Iteration ' num2str(j) ' out of ' num2str(n)]);
        fprintf(1, '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b', n);
        count = 0;
    elseif count == 100 && cct == 1 && n >= 1000 && n < 10000 && j >= 1000
        display(['Iteration ' num2str(j) ' out of ' num2str(n)]);
        fprintf(1, '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b', n);
        count = 0;
    elseif count == 100 && cct == 1 && n >= 10000 && j < 1000
        display(['Iteration ' num2str(j) ' out of ' num2str(n)]);
        fprintf(1, '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b', n);
        count = 0;
    elseif count == 100 && cct == 1 && n >= 10000 && j < 10000 && j >= 1000
        display(['Iteration ' num2str(j) ' out of ' num2str(n)]);
        fprintf(1, '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b', n);
        count = 0;
    elseif count == 100 && cct == 1 && n >= 10000 && j >= 10000
        display(['Iteration ' num2str(j) ' out of ' num2str(n)]);
        fprintf(1, '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b', n);
        count = 0;
    end
    
    if cct == 1
        count = count + 1;
    end
    
    % fetch grid point
    G = G_values(j+1);
    
    % Set up cost function 
    if use_summation
        % Use summation approach for deterministic case when branch=1
        Lossfunction = @(t) lossfunction_summation_based(a, b, n, j, x, t, 1, dissaving, ELB, dis, u, f, offset, opts, con);
    else
        % Use value function iteration approach (for stochastic or when branch=0)
        Lossfunction = @(t) Lossfunction_vfi(a, b, n, j, x, t, G, u, beta, delta, p, tol, f, offset, opts, con, income_values(j+1));
    end
    
    % Set lower and upper bounds for search
    if j == 0
        lb = a;
    else
        lb = max(x(j), ELB); % Markov policies are weakly increasing
    end
    
    % fetch upper bound
    ub = ub_values(j+1);
    
    % minimize loss (note that we are ignoring /(n+1) in definition of Lossfunction as this is irrelevant
    % for optimization)
    [amin, cfm] = fminbnd(Lossfunction, lb, ub);
    
    % Store response based on which method is being used
    if use_summation
        idr(j+1) = lossfunction_summation_based(a, b, n, j, x, amin, 2, dissaving, ELB, dis, u, f, offset, opts, con);
    else
        % Calculate best response using VFI approach
        fun = @(z) stochastic_best_response(a, b, n, j, x, amin, G, z, u, beta, delta, p, tol, f, offset, opts, con, income_values(j+1));
        [amax, ~] = fminbnd(fun, lb, ub, optimset('Display', 'off'));
        idr(j+1) = amax;
    end
    
    % Extra precision check against bounds
    if cfm >= Lossfunction(lb)
        amin = lb;
        cfm = Lossfunction(lb);
    end
    
    if dissaving == 2 && cfm >= Lossfunction(ub)
        amin = ub;    % precision check against upper bound and move lower bound if ub is optimal (captures steady states and improves precision)
    end
    
    % Set optimal saving
    x(j+1) = amin;
     
end
x=transpose(x);
idr=transpose(idr);
end


function [lossv] = lossfunction_summation_based(a, b, n, j, s, t, out, dissaving, ELB, dis, u, f, offset, opts, con, income_values)
    % out=1: Returns Loss, taking s(j+1)=s(j+2)=..s(n+1)=t
    % out=2: Second-mover policy
    % out=3: Returns Loss given continuation s
    
    cutoff = length(dis);  
    
    G = a + j * ((b - a) / n); % current wealth at grid point j
    fG = f(G); % simplifying notation slightly
    
    
    if out == 1
        % Create a copy of s to modify
        s_temp = s;
        
        % Always set the current grid point to t
        s_temp(j+1) = min(t, b);
        
        % Set future policy based on con parameter
        if con == 0
            s_temp(j+2:n+1) = min(t, b);            
        else
            s_temp(j+2:n+1) = opts(j+2:n+1);            
        end
        
        % Replace s with s_temp for the remainder of the function
        s = s_temp;
    end
    
    lb = ELB;
    if dissaving > 0
        ub = G;
    else
        ub = min(fG, b);
    end
    
    options = optimset('Display', 'off');
    
    fun = @(z) -u(fG - z) - V_deterministic(a, b, n, z, s, ELB, dis, u, f);
    [amax, vmin] = fminbnd(fun, lb, ub, options);
    
    lossv = -vmin - u(fG - s(j+1)) - V_deterministic(a, b, n, s(j+1), s, ELB, dis, u, f, offset);
    
    if out == 2
        lossv = amax;
    end
end

% Deterministic value function
function [Vt] = V_deterministic(a, b, n, y, s, ~, dis, u, f, ~)
    cutoff = length(dis);  
    delta = dis(2) / dis(1);
    beta = dis(1) / delta;
    run = funcen(s, y, a, b, n);
    val = dis(1) * u(f(y) - run);
    
    if cutoff > 1
        for j = 2:1:cutoff
            prev = run;
            run = funcen(s, prev, a, b, n);
            val = val + dis(j) .* u(f(prev) - run);
        end
    end
    
    cons = (delta^(cutoff + 1)) / (1 - delta);
    val = val + beta * cons * u(f(prev) - run);
    Vt = val;
end

% Stochastic loss function via vfi
function [loss] = Lossfunction_vfi(a, b, n, j, s, t, G, u, beta, delta, p, tol, f, offset, opts, con, fG)
    % This function mirrors the deterministic approach by creating a hypothetical
    % policy where all future choices follow t or opts, then computing the value function

    % Create temporary policy function
    s_temp = s;
    
    % Always set the current grid point to t
    s_temp(j+1) = min(t, b);
    
    % Set future policy based on con parameter - use loops for safety
    if con == 0
        % Default: All future choices follow t
        for r = j+2:1:n+1
            s_temp(r) = min(t, b);
        end
    else
        % Alternative: All future choices follow opts
        for r = j+2:1:n+1
            s_temp(r) = opts(r);
        end
    end
    
    % Get grid and initialize value function
    gridFine = linspace(a, b, n+1);
    W_init = u(max(0.1, f(a) - a)) * ones(size(gridFine));
    
    % Compute value function based on this policy
    W = compute_value_function(a, b, n, u, p, delta, tol, s_temp, f, offset, W_init, gridFine);
    
    % Create interpolation function
    VF = @(z) interp1(gridFine, W, real(z), 'linear', 'extrap');
    
    % safeguard (if lossfunctions are changed!)
    if nargin > 16 && ~isempty(fG)
        income_G = fG;
    else
        income_G = f(G);
    end
    
    % Find best response at current grid point/asset level
    fun = @(z) -u(income_G - z) - beta * delta * (p * VF(z) + (1 - p) * VF(offset(z)));
    [amax, vmin] = fminbnd(fun, a, min(income_G, b), optimset('Display', 'off'));
    
    % second-mover payoff
    utility_t = u(income_G - t) + beta * delta * (p * VF(t) + (1 - p) * VF(offset(t)));
    utility_opt = u(income_G - amax) + beta * delta * (p * VF(amax) + (1 - p) * VF(offset(amax)));
    
    % loss
    loss = utility_opt - utility_t;
    
    % If loss is negative, set to zero (can happen if second-mover optimum
    % is not actually an optimum)
    if loss < 0
        loss = 0;
    end
end

% For computing best response/second-mover policy in stochastic model
function [loss] = stochastic_best_response(a, b, n, j, s, t, G, z, u, beta, delta, p, tol, f, offset, opts, con, fG)
    s_temp = s;    
    s_temp(j+1) = min(t, b);
    
    % Set future policy based on con parameter - using loop for safety
    if con == 0
        % Default: All future choices follow t
        for r = j+2:1:n+1
            s_temp(r) = min(t, b);
        end
    else
        % Alternative: All future choices follow opts
        for r = j+2:1:n+1
            s_temp(r) = opts(r);
        end
    end
    
    % Get grid and initialize value function
    gridFine = linspace(a, b, n+1);
    W_init = u(max(0.1, f(a) - a)) * ones(size(gridFine));
    
    % Compute value function based on this policy
    W = compute_value_function(a, b, n, u, p, delta, tol, s_temp, f, offset, W_init, gridFine);
    
    % Create interpolation function
    VF = @(zz) interp1(gridFine, W, real(zz), 'linear', 'extrap');
    
    income_G = fG;
    
    
    % Return loss for minimization
    loss = -u(income_G - z) - beta * delta * (p * VF(z) + (1 - p) * VF(offset(z)));
end

% Compute value function (main valf from stochastic models)
function [W_out] = compute_value_function(a, b, n, u, p, delta, tol, s, f, offset, W_initial, gridFine)
    iterlimit = 100;
    
    if tol < 1e-10
        iterlimit = 1000;
    end
    
    % Use persistent variables to store previous value function
    persistent W_old_persistent grid_persistent
    if isempty(W_old_persistent) || length(grid_persistent) ~= length(gridFine) || any(grid_persistent ~= gridFine)
        W_old = W_initial;
    else
        W_old = W_old_persistent;
    end
    
    % Value function iteration
    iter = 0;
    dist = Inf;
    originalGrid = linspace(a, b, n+1);
    
    % Pre-calculate function values once
    fGridFine = f(gridFine);
    
    while dist > tol && iter < iterlimit
        % policy function
        gVals = interp1(originalGrid, s, gridFine, 'linear', 'extrap');
        gVals = max(a, gVals);
        
        % Extrapolate points above b
        above_b = gridFine > b;
        if any(above_b)
            num_segments = min(10, length(s) - 1);
            segment_slopes = diff(s(end-num_segments:end)) ./ diff(originalGrid(end-num_segments:end));
            avg_slope = mean(segment_slopes);
            extrap_vals = s(end) + avg_slope * (gridFine(above_b) - originalGrid(end));
            gVals(above_b) = extrap_vals;
        end
        
        % Calculate consumption and values
        cons = fGridFine - gVals;
        vLow = interp1(gridFine, W_old, gVals, 'linear', 'extrap');
        
        % Calculate offset values
        offsetGVals = offset(gVals);
        vHigh = interp1(gridFine, W_old, offsetGVals, 'linear', 'extrap');
        
        % New value function with dampening - vectorized operation
        W_new = 0.5 * W_old + 0.5 * (u(cons) + delta * (p * vLow + (1 - p) * vHigh));
        
        % Check convergence
        dist = max(abs(W_new - W_old));
        W_old = W_new;
        
        iter = iter + 1;
    end
    
    W_out = W_old;
    
    % Update persistent storage
    W_old_persistent = W_out;
    grid_persistent = gridFine;
end

% Interpolation function (handrolled for speed)
function [funcetn] = funcen(ss, xx, a, b, n)
    if xx > b
        funcetn = max(ss);
    else
        ngp = max(min(ceil((xx - a) / ((b - a) / n)), n), 0);
        if ngp == 0
            funcetn = ss(1);
        else
            alp = ngp + ((a - xx) .* n ./ (b - a));
            funcetn = alp .* ss(ngp) + (1 - alp) .* ss(ngp + 1);
        end
    end
end