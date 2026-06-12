function [opts] = algorithm_call_iid(a, b, n, MaxIt, innerstep, opts, RESTART, dissaving, dis, p, u, tol, f, offset, dsp, dyncon, sm)
% CMA-ES optimization on non-mandated ("active") policy components where assets > sm
% Uses Hansen's Cholesky update for covariance matrix adaptation
% - xfromz: maps z ∈ R^{active_dims} to weakly non-decreasing sequence in [a,b]
% - zfromx: inverse mapping for initialization and best response injection
tol_a = 1e-8;  % tolerance for proximity to borrowing limit
wtc=0;  % writes additional information to the console
WB = 1; % enable waitbar
opts=opts(:);
dyncon=dyncon(:);

% Compute grid points G(j)
G = linspace(a, b, n+1)';

% determines non-mandated ("active") and mandated indices
if sm > a
    active_idx = (G > sm);  
    opts(~active_idx) = dyncon(~active_idx);  
else
    active_idx = true(n+1, 1);
end

idxA = find(active_idx);          % numeric indices of active (free) points
nA   = numel(idxA);               % number of free points
nVar = nA;

% Define cost function
CostFunction = @(x) Lossfunction(a, b, n, x, dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);

% Initial candidate generation
if RESTART == 1
    % Initial evaluation 
    opts(~active_idx) = dyncon(~active_idx);
    [cost, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);
    
    for m = 1:100
        Br = Br(:);
        Br(~active_idx) = dyncon(~active_idx);
        [newcost, newBr] = Lossfunction(a, b, n, Br, dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);
        if newcost < cost
            opts = Br;
            cost = newcost;
            Br = newBr;
        else
            break;
        end        
    end    
end

% Create z_init from opts, but limit to active indices for the transform
z_init = zfromx(opts(idxA), a, b, nA-1, tol_a);

% Store initial loss for reporting
initial_loss = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);

lambda = 8; 
mu = round(lambda/2);
w = log(mu + 0.5) - log(1:mu);
w = w / sum(w);
mu_eff = 1 / sum(w.^2);

sigma = innerstep; 
cs = (mu_eff + 2) / (nVar + mu_eff + 5);
ds = 1 + cs + 2*max(sqrt((mu_eff - 1)/(nVar+1)) - 1, 0);
ENN = sqrt(nVar)*(1 - 1/(4*nVar) + 1/(21*nVar^2));
cc = (4 + mu_eff / nVar) / (4 + nVar + 2*mu_eff / nVar);
c1 = 2 / ((nVar + 1.3)^2 + mu_eff);
alpha_mu = 2;
cmu = min(1 - c1, alpha_mu*(mu_eff - 2 + 1/mu_eff) / ((nVar + 2)^2 + alpha_mu*mu_eff/2));
hth = (1.4 + 2/(nVar+1))*ENN;

% Evolution path variables and Cholesky factor
ps = zeros(1, nVar);
pc = zeros(1, nVar);
A = eye(nVar);

% Current best solution
current_pos = z_init;
current_step = zeros(1, nVar);
x = dyncon;
x(idxA) = xfromz(z_init, a, b, nA-1);
current_cost = CostFunction(x);

if current_cost < 0
    disp('Negative initial Ego Loss. Stopping.');
    MaxIt = 0;
end

best_pos = current_pos;
best_cost = current_cost;
BestCost = zeros(MaxIt, 1);
stall_counter = 0;

if WB == 1
    hWaitbar = waitbar(0, 'Initiating', ...
        'Name', 'Solving (cancel to skip)', 'CreateCancelBtn','delete(gcbf)');
    set(hWaitbar, 'Units','Pixels','Position',[50 500 380 100]);
end

brcount = 0;
cached_parent_cost = 0;

% Preallocate (for parfor)
candidate_costs = zeros(lambda, 1);
candidate_steps = zeros(lambda, nVar);
candidate_positions = zeros(lambda, nVar);

for g = 1:MaxIt
    % Sample using Cholesky factor
    Z = randn(lambda, nVar);
    stepsAll = Z * A';
    
    brcount = brcount + 1;
    if brcount == 2
        specialCandidateIndex = lambda;
    else
        specialCandidateIndex = -1;
    end

    % Parallel candidate evaluation with preallocated arrays
    parfor i = 1:lambda
        if (i == specialCandidateIndex)
            % Best response injection
            x = dyncon;
            x(idxA) = xfromz(best_pos, a, b, nA-1);
            [~, Br] = Lossfunction(a, b, n, x, dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);
            Br = Br(:); 
            Br(~active_idx) = dyncon(~active_idx);
            candidate_positions(i,:) = zfromx(Br(idxA), a, b, nA-1, tol_a);
            candidate_steps(i,:) = zeros(1, nVar);
        else
            candidate_steps(i,:) = stepsAll(i,:);
            candidate_positions(i,:) = current_pos + sigma * candidate_steps(i,:);
        end
        x = dyncon;
        x(idxA) = xfromz(candidate_positions(i,:), a, b, nA-1);
        candidate_costs(i) = CostFunction(x);
    end
    
    if specialCandidateIndex == lambda
        brcount = 0;
    end

    % Find best candidate
    [min_cost, best_idx] = min(candidate_costs);
    
    % Check only the best candidate against current best
    if min_cost < best_cost && isreal(min_cost) && (min_cost > 0)
        best_pos = candidate_positions(best_idx,:);
        best_cost = min_cost;
    end
    
    BestCost(g) = best_cost;

    % Extract final solution
    x = dyncon;
    x(idxA) = xfromz(best_pos, a, b, nA-1);
    opts = x;

    if g > 1
        if BestCost(g)/nVar<1e-14
                disp(['Ego Loss equals zero within tolerance (' num2str(BestCost(g)/nVar) '). Finishing']);
                if WB==1, close(hWaitbar); end
                break;
        end
        
        if WB==1 && ~ishandle(hWaitbar)
            disp('Canceled by user');
            g = MaxIt;
            break;
        elseif WB==1
            waitbar((g-1)/MaxIt, hWaitbar, ...
                ['Ego Loss = ' num2str(BestCost(g-1)/nVar) ...
                 ' (Iteration ' num2str(g-1) ' of ' num2str(MaxIt) ')']);
        end
    end

    if g == MaxIt
        disp(['Iteration ' num2str(g) ' [final] Ego Loss = ' num2str(BestCost(g)/nVar)]);
        if WB==1, close(hWaitbar); end
        break;
    end

    % Weighted step calculation
    weighted_step = zeros(1, nVar);
    for jj = 1:mu
        % Get jj-th best candidate by sorting only indices we need
        if jj == 1
            step_to_add = candidate_steps(best_idx,:);
        else
            [~, sorted_indices] = sort(candidate_costs);
            step_to_add = candidate_steps(sorted_indices(jj),:);
        end
        weighted_step = weighted_step + w(jj) * step_to_add;
    end
    
    new_pos = current_pos + sigma * weighted_step;
    x = dyncon;
    x(idxA) = xfromz(new_pos, a, b, nA-1);
    testVal = CostFunction(x);
    
    if ~isreal(testVal)
        disp('Imaginary candidate encountered.');
    end
    
    % Validity check 
    if testVal > 0 && isreal(testVal)
        current_pos = new_pos;
        current_step = weighted_step;
        current_cost = testVal;
        if testVal < best_cost
            best_pos = current_pos;
            best_cost = testVal;
        end
    else
        current_step = weighted_step;  % Keep step for path updates
        if g == 1
            cached_parent_cost = current_cost;
        end
        current_cost = cached_parent_cost;
    end

    % Evolution path updates
    ps = (1-cs)*ps + sqrt(cs*(2-cs)*mu_eff)* (current_step / A');
    sigma = sigma * exp(cs/ds*(norm(ps)/ENN - 1));

    if norm(ps) / sqrt(1-(1-cs)^(2*g)) < hth
        hs = 1;
    else
        hs = 0;
    end
    delta_c = (1 - hs)*cc*(2 - cc);
    pc = (1-cc)*pc + hs*sqrt(cc*(2-cc)*mu_eff)*current_step;

    % Hansen's Cholesky update
    A = sqrt(1 - c1 - cmu + delta_c*c1) * A;
    
    % Rank-one update for evolution path
    if norm(pc) > 0
        A = cholupdate(A, sqrt(c1) * pc', '+');
    end
    
    % Rank-μ updates for selected offspring
    [~, sorted_indices] = sort(candidate_costs);
    for jj = 1:mu
        step_j = candidate_steps(sorted_indices(jj),:);
        if w(jj) > 0 && norm(step_j) > 0
            A = cholupdate(A, sqrt(cmu * w(jj)) * step_j', '+');
        end
    end

    % Stall detection and restart
    if g > 1
        if BestCost(g) >= BestCost(g-1)
            stall_counter = stall_counter + 1;
        else
            stall_counter = 0;
        end
    end
    if stall_counter >= 20
        innerstep = innerstep/2;
        sigma = innerstep;
        % Restart: reinitialize
        A = eye(nVar);
        ps = zeros(1, nVar);
        pc = zeros(1, nVar);
        current_pos = best_pos;
        current_cost = best_cost;
        stall_counter = 0;
        if wtc==1
            disp(['Reducing innerstep to ' num2str(innerstep) ' and restarting']);
        end
    end
end

% Final solution extraction
x = dyncon;
x(idxA) = xfromz(best_pos, a, b, nA-1);
opts = roundToA(x, a, tol_a);

% Report initial and final loss
final_loss = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 0, tol, f, offset, dsp, sm);
fprintf('Initial loss %.6e; Final loss %.6e\n', initial_loss/(n+1), final_loss/(n+1));
end

% Maps z ∈ R^{active_dims} to weakly non-decreasing sequence in [a,b]
function x = xfromz(z, a, b, n)
    x1 = a + (b - a) / (1 + exp(-z(1)));
    alphas = exp(z(2:end));
    S = sum(alphas);
    remaining = (b - x1) * S / (S + 1);
    x_inc = (remaining / S) * alphas;
    x = [x1, x1 + cumsum(x_inc)];
end

% Maps weakly non-decreasing sequence in [a,b] to z ∈ R^{active_dims}
function z = zfromx(x, a, b, n, tol_a)
    z = zeros(1, n+1);
    if x(1) <= a + tol_a
        z(1) = -20;  
    else
        frac = (x(1) - a) / (b - a);
        z(1) = log(frac / (1 - frac));
    end
    
    D = max(diff(x), 0);
    leftover = b - x(1);
    total_inc = sum(D);
    S = total_inc / (leftover - total_inc);
    alphas = (S + 1) * D / leftover;
    z(2:end) = log(max(alphas, 1e-10));
end

function vals = roundToA(vals, a, tol_a)
    closeToA = abs(vals - a) < tol_a;
    vals(closeToA) = a;
end