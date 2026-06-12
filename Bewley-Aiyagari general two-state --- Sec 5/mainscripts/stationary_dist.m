function [Mean, asim, lsim, Tsim] = stationary_dist(a, b, n, ~, ~, lmin, lmax, P, s, ~, ~, power_param)
% STATIONARY_DIST Compute theoretical stationary distribution for Aiyagari-type models
%
% Computes exact stationary distribution using Markov chain iteration on a joint
% (asset, labor productivity) state space. Handles general 2-state labor productivity 
% processes via transition matrix P. Initializes all mass at borrowing constraint.
%
% Key features:
% - Fine grid refinement to handle discontinuous policy functions
% - Handles agents at borrowing constraint exactly
% - Grid-consistent moment calculations
% - Policy boundary check: errors if >5% of policy points exit asset grid meaningfully
% - Includes Monte Carlo simulation of asset paths on original grid
%
% Inputs: asset bounds [a,b], grid size n, power grid parameter, 
%         productivity bounds [lmin,lmax], transition matrix P, policy functions s

% Build asset grid
t_full = linspace(0, 1, n+1);
agrid = a + (b - a) * (t_full.^power_param)';

% Extract policies for low/high labor productivity states
if size(s, 1) == 2*(n+1)
    s_low = s(1:n+1);
    s_high = s(n+2:end);
else
    s_low = s;
    s_high = s;
end

% Check for problematic policy violations (>5% of points)
grid_width = b - a;
tol_meaningful = 0.001 * grid_width;
violations_low = sum(s_low < a - tol_meaningful | s_low > b + tol_meaningful);
violations_high = sum(s_high < a - tol_meaningful | s_high > b + tol_meaningful);

if violations_low > 0.05 * length(s_low) || violations_high > 0.05 * length(s_high)
    fprintf('Warning: %d policy violations in s_low, %d in s_high\n', violations_low, violations_high);
    fprintf('Max s_low violation: %.6f, Max s_high violation: %.6f\n', ...
            max([max(s_low) - b, a - min(s_low), 0]), ...
            max([max(s_high) - b, a - min(s_high), 0]));
    error('More than 5%% of policy points exit grid meaningfully; extend grid or tighten solver');
end

% Store policies in matrix form
sav = [s_low, s_high];

%% Theoretical calculation on finer grid
% Use denser grid to reduce discretization bias in transition matrix
% Critical for discontinuous policies: when s(a) jumps over [a_j, a_k], 
% coarse transition matrix has zero-mass rows for j+1,...,k-1.
% Fine grid reduces |s(a_i+1) - s(a_i)|, improving interpolation accuracy
% in transition probability allocation.
fgridm = 2; % Fine grid multiplier
na_fine = fgridm * (n + 1);
t_fine = linspace(0, 1, na_fine);
agrid_fine = a + (b - a) * (t_fine.^power_param)';

% Interpolate policies to fine grid
savinterp_fine = cell(1,2);
for il = 1:2
    savinterp_fine{il} = griddedInterpolant(agrid, sav(:, il), 'linear', 'nearest');
end
sav_fine = [savinterp_fine{1}(agrid_fine), savinterp_fine{2}(agrid_fine)];

% reapply bounds in case of floating-point precision errors
a_next_low = max(a, min(b, sav_fine(:, 1)));
a_next_high = max(a, min(b, sav_fine(:, 2)));

%% Build joint transition matrix
% Get labor productivity stationary distribution
[V_eig, D_eig] = eig(P');
[~, idx] = min(abs(diag(D_eig)-1));
ldist = V_eig(:, idx);
ldist = ldist / sum(ldist);

% Initialize joint transition matrix
total_states = na_fine * 2;
Q_joint = sparse(total_states, total_states);
state_idx = @(a_idx, l_state) a_idx + (l_state-1)*na_fine;

% Fill transition probabilities for low labor productivity state
a_next_low = sav_fine(:, 1);
for i = 1:na_fine
    from_state = state_idx(i, 1);
    
    if abs(a_next_low(i) - a) < 1e-12
        % Stay at borrowing constraint
        Q_joint(from_state, state_idx(1, 1)) = Q_joint(from_state, state_idx(1, 1)) + P(1,1);
        Q_joint(from_state, state_idx(1, 2)) = Q_joint(from_state, state_idx(1, 2)) + P(1,2);
    else
        % Normal transition via interpolation
        [j_low, weights_low] = interpolation_indices(a_next_low(i), agrid_fine);
        
        % To low labor productivity
        Q_joint(from_state, state_idx(j_low, 1)) = Q_joint(from_state, state_idx(j_low, 1)) + P(1,1)*(1-weights_low);
        Q_joint(from_state, state_idx(min(j_low+1, na_fine), 1)) = Q_joint(from_state, state_idx(min(j_low+1, na_fine), 1)) + P(1,1)*weights_low;
        
        % To high labor productivity
        Q_joint(from_state, state_idx(j_low, 2)) = Q_joint(from_state, state_idx(j_low, 2)) + P(1,2)*(1-weights_low);
        Q_joint(from_state, state_idx(min(j_low+1, na_fine), 2)) = Q_joint(from_state, state_idx(min(j_low+1, na_fine), 2)) + P(1,2)*weights_low;
    end
end

% Fill transition probabilities for high labor productivity state
a_next_high = sav_fine(:, 2);
for i = 1:na_fine
    from_state = state_idx(i, 2);
    
    if abs(a_next_high(i) - a) < 1e-12
        % Stay at borrowing constraint
        Q_joint(from_state, state_idx(1, 1)) = Q_joint(from_state, state_idx(1, 1)) + P(2,1);
        Q_joint(from_state, state_idx(1, 2)) = Q_joint(from_state, state_idx(1, 2)) + P(2,2);
    else
        % Normal transition via interpolation
        [j_high, weights_high] = interpolation_indices(a_next_high(i), agrid_fine);
        
        % To low labor productivity
        Q_joint(from_state, state_idx(j_high, 1)) = Q_joint(from_state, state_idx(j_high, 1)) + P(2,1)*(1-weights_high);
        Q_joint(from_state, state_idx(min(j_high+1, na_fine), 1)) = Q_joint(from_state, state_idx(min(j_high+1, na_fine), 1)) + P(2,1)*weights_high;
        
        % To high labor productivity
        Q_joint(from_state, state_idx(j_high, 2)) = Q_joint(from_state, state_idx(j_high, 2)) + P(2,2)*(1-weights_high);
        Q_joint(from_state, state_idx(min(j_high+1, na_fine), 2)) = Q_joint(from_state, state_idx(min(j_high+1, na_fine), 2)) + P(2,2)*weights_high;
    end
end

% Normalize rows to ensure mass conservation
Q_joint = spdiags(1./sum(Q_joint,2), 0, total_states, total_states) * Q_joint;

%% Solve for stationary distribution
% Initialize with all mass at borrowing constraint
pi_joint_old = zeros(1, total_states);
pi_joint_old(state_idx(1, 1)) = ldist(1);
pi_joint_old(state_idx(1, 2)) = ldist(2);
pi_joint_old = pi_joint_old / sum(pi_joint_old);

% Power iteration
tol = 1e-12;
max_iter = 1e7;
for iter = 1:max_iter
    pi_joint_new = pi_joint_old * Q_joint;
    if max(abs(pi_joint_new - pi_joint_old)) < tol
        break;
    end
    pi_joint_old = pi_joint_new;
end
if iter == max_iter
    warning('Stationary distribution did not converge');
    disp('If there are no signs of artefacts in policy functions, reduce tol and/or raise max_iter in stationary_dist')
end
pi_joint = pi_joint_new / sum(pi_joint_new);

% Extract distributions by labor productivity state
pi_low = pi_joint(1:na_fine);
pi_high = pi_joint(na_fine+1:2*na_fine);

%% Calculate moments directly on fine grid
% Direct computation on fine grid avoids binning artifacts when policies
% have discontinuities. Binning pi_fine to coarse grid can assign zero mass
% to entire coarse bins if they fall within policy jump regions.
mean_assets_fine = (pi_low + pi_high)* agrid_fine;
Mean = mean_assets_fine; % labor productivities normalized to 1 in stationary distribution

%% Monte Carlo simulation  
burn_in = 1500;
Nsim = 50000;           % Must match Nsim in figures_random.m for consistent results  
Tsim_full = 2000;
Tsim = Tsim_full - burn_in;
rng(2025);              % random seed (for reproducability)

% Labor productivity grid
lgrid = [lmin, lmax]';
nl = length(lgrid);

 
% Generate all random numbers at once
u_all = rand(Nsim, Tsim_full);

% Initialize productivity indices
lsim_idx_full = zeros(Nsim, Tsim_full);
lsim_idx_full(:, 1) = 1 + (u_all(:, 1) > ldist(1));

% Vectorized Markov transitions
for it = 2:Tsim_full
    % Low state transitions: stay in 1 with prob P(1,1), move to 2 with prob P(1,2)
    low_mask = (lsim_idx_full(:, it-1) == 1);
    lsim_idx_full(low_mask, it) = 1 + (u_all(low_mask, it) > P(1,1));  % > P(1,1) means go to state 2

    % High state transitions: move to 1 with prob P(2,1), stay in 2 with prob P(2,2)  
    high_mask = (lsim_idx_full(:, it-1) == 2);
    lsim_idx_full(high_mask, it) = 1 + (u_all(high_mask, it) > P(2,1));  % > P(2,1) means stay in state 2
end

% Convert indices to actual productivities
lsim_full = lgrid(lsim_idx_full);

% Pre-allocate and vectorize
asim_full = zeros(Nsim, Tsim_full);
asim_full(:, 1) = a;

% Create interpolants 
savinterp = cell(1,nl);
savinterp{1} = griddedInterpolant(agrid, s_low, 'linear', 'nearest');
savinterp{2} = griddedInterpolant(agrid, s_high, 'linear', 'nearest');

% Asset evolution 
for it = 1:Tsim_full-1
    current_assets = asim_full(:, it);
    
    % Separate by productivity state 
    low_mask = (lsim_idx_full(:, it) == 1);
    high_mask = (lsim_idx_full(:, it) == 2);
    
    % Vectorized interpolation by state
    asim_full(low_mask, it+1) = savinterp{1}(current_assets(low_mask));
    asim_full(high_mask, it+1) = savinterp{2}(current_assets(high_mask));
end

% Return post-burn-in simulation
asim = asim_full(:, burn_in+1:end);
lsim = lsim_full(:, burn_in+1:end);

Meansim = mean(asim(:)) / mean(lsim(:));
% Calculate percentage difference
pct_diff = 100 * abs(Mean - Meansim) / Meansim;
if pct_diff>10
    disp(['Deviation between theoretical and simulation mean: ' num2str(pct_diff, '%.1f') '%']);
end
Mean=Meansim;
 
 
% Gini
w = (pi_low + pi_high); 
w = w(:)/sum(w(:));
x = agrid_fine(:);
L = cumsum(w .* x) / sum(w .* x);
Gini = 1 - sum(w .* ([0; L(1:end-1)] + L));
 % Palma ratio
cw = cumsum(w);
S = w .* x;
tot = sum(S);
i40 = find(cw >= 0.4, 1);
cw_prev = 0; 
if i40 > 1, cw_prev = cw(i40-1); end
share40 = (sum(S(1:i40-1)) + (0.4 - cw_prev) * x(i40)) / tot;
i90 = find(cw >= 0.9, 1);
cw_prev90 = 0; 
if i90 > 1, cw_prev90 = cw(i90-1); end
shareTop10 = (sum(S(i90+1:end)) + (cw(i90) - 0.9) * x(i90)) / tot;
palma = shareTop10 / share40;

%for debugging purposes or to rely instead on theoretical distribution
%disp(['Gini (theory) ' num2str(Gini) ' Palma ratio  ' num2str(palma) ' Fraction constrained ' num2str(w(1))]);
 

%Mean=Meansim;
end

%% Helper functions
function [j, weights] = interpolation_indices(a_next, grid)
    n = length(grid);
    j = ones(size(a_next));
    for i = 2:n
        j(a_next >= grid(i-1) & a_next <= grid(i)) = i-1;
    end
    j(a_next > grid(end)) = n-1;
    next_j = min(j+1, n);
    denom = grid(next_j) - grid(j);
    weights = (a_next - grid(j)) ./ denom;
    weights(denom == 0) = 0;
    weights = min(max(weights, 0), 1);
end

function [yi] = lininterp1_vec(x, y, xi)
    yi = zeros(size(xi));
    for i = 1:length(xi)
        placeLow = find(xi(i) < x, 1) - 1;
        if placeLow == 0
            placeLow = 1;
        end
        if isempty(placeLow)
            placeLow = length(x) - 1;
        end
        placeHigh = placeLow + 1;
        xLow = x(placeLow);
        xHigh = x(placeHigh);
        yLow = y(placeLow);
        yHigh = y(placeHigh);
        yi(i) = yLow + (xi(i) - xLow) * (yHigh - yLow) / (xHigh - xLow);
    end
end