function [opts] = algorithm_call_default(a, b, n, MaxIt, innerstep, opts, RESTART, dissaving, dis, P, u, tol, f_low, f_high, dsp, power_param, tolalg, G, IncomeG_low, IncomeG_high)
% Standard CMA-ES algorithm 
% Uses Cholesky factor approach (Hansen, 2006) for full covariance matrix
% Optimizes monotonic savings policies on asset-income grids

% Algorithm control parameters 
wtc = 0;					% =1 enable console output
WB = 1;						% =1 enable waitbar          
policyinj = 1;				% =1 enable policy injections 
gradinj = 1;				% =1 enable gradient injections


%% Store initial loss
[initial_loss, ~] = Lossfunction(a, b, n, opts, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);

%% Policy iteration restart
if RESTART == 1
	for m = 1:100
		[old_cost, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);
		opts2 = Br(:);        
		[new_cost, ~] = Lossfunction(a, b, n, opts2, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);
		if new_cost < old_cost
			opts = opts2;
		else
			break;
		end
	end
end

if RESTART == 2		% Single policy iteration step for comparative statics
	[~, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);        
	opts = transpose(Br(:));       
end
opts = opts(:);

%% Cost function for CMA-ES optimization
CostFunction = @(z) Lossfunction(a, b, n, ztox(z, IncomeG_low, IncomeG_high, a, b, n), dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);



nVar = 2 * (n + 1);

%% Initialize CMA-ES parameters
z_init = xtoz(opts, IncomeG_low, IncomeG_high, a, b, n);

% CMA-ES configuration
lambda = 4;
mu = round(lambda/2);
w = log(mu + 0.5) - log(1:mu);
w = w / sum(w);
mu_eff = 1 / sum(w.^2);

% Standard CMA-ES parameters
c_c = 4 / (nVar + 4);
c_1 = 2 / ((nVar + 1.3)^2 + mu_eff);
alpha_mu = 2;
c_mu = min(1 - c_1, alpha_mu * (mu_eff - 2 + 1/mu_eff) / ((nVar + 2)^2 + alpha_mu*mu_eff/2));

sigma0 = innerstep; 
cs = (mu_eff + 2) / (nVar + mu_eff + 5);
ds = 1 + cs + 2*max(sqrt((mu_eff - 1)/(nVar+1)) - 1, 0);
ENN = sqrt(nVar)*(1 - 1/(4*nVar) + 1/(21*nVar^2));
cc = c_c;
c1 = c_1;
cmu = c_mu;
hth = (1.4 + 2/(nVar+1))*ENN;

% CMA-ES state variables using Cholesky factor approach
ps_curr = zeros(1, nVar);
pc_curr = zeros(1, nVar);
A_curr = eye(nVar);				% Cholesky factor: C = A * A^T
sigma_curr = sigma0;

% State tracking
M_pos = z_init(:);
M_step = zeros(nVar, 1);
M_cost = CostFunction(M_pos);

BestSol_pos = M_pos;
BestSol_cost = M_cost;
BestCost_prev = Inf;
stall_counter = 0;

% Pre-allocate for candidates (reused each iteration)
pop_costs = zeros(lambda, 1);
pop_steps = zeros(lambda, nVar);
pop_positions = zeros(lambda, nVar);
S = zeros(nVar, lambda);  % Pre-allocate for all samples

% Additional pre-allocations for performance (Optimization #3)
Z_storage = zeros(lambda, nVar);
weighted_steps_storage = zeros(mu, nVar);

if WB == 1
	hWaitbar = waitbar(0, 'Initiating', ...
		'Name', 'Solving (cancel to skip)', 'CreateCancelBtn','delete(gcbf)');
	set(hWaitbar, 'Units','Pixels','Position',[50 500 380 100]);
end

%% Main CMA-ES loop
brcount = 0;
grad_alpha = 1.0;			% Gradient step size parameter
for g = 1:MaxIt
	% Generate candidates using Cholesky factor (standard Gaussian) - vectorized
	Z = randn(lambda, nVar);
	stepsAll = Z * A_curr.';
	candMat = M_pos' + sigma_curr * stepsAll;
	
	brcount = brcount + 1;
	
	% Handle special candidates outside parallel loop
	% Pre-allocate for up to 2 special candidates (policy + gradient)
	max_special = 2;
	special_costs = zeros(max_special, 1);
	special_steps = zeros(max_special, nVar);
	special_positions = zeros(max_special, nVar);
	num_special = 0;
	
	% Policy injection
	if brcount == 2 && policyinj == 1
		best_policy = ztox(BestSol_pos, IncomeG_low, IncomeG_high, a, b, n);
		
		[~, Br] = Lossfunction(a, b, n, best_policy, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);
		
		policy_cand = xtoz(Br, IncomeG_low, IncomeG_high, a, b, n);
		policy_step = zeros(nVar, 1);
		policy_cost = CostFunction(policy_cand);
		
		% Report success
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
			
			% Report success
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
		
		brcount = 0;  % Reset after gradient injection
	end
	
	% Fixed population size issue (Issue #4)
	% Always evaluate lambda candidates total
	lambda_standard = lambda - num_special;
	
	% Generate standard candidates in parallel - vectorized
	standard_indices = 1:lambda_standard;
	standard_steps = stepsAll(standard_indices, :);
	standard_positions = candMat(standard_indices, :);
	
	% Evaluate costs in parallel
	standard_costs = zeros(lambda_standard, 1);
	for i = 1:lambda_standard
		standard_costs(i) = CostFunction(standard_positions(i, :)');
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

	% Update best solution - vectorized
	valid_mask = pop_costs < BestSol_cost & isreal(pop_costs) & pop_costs > 0;
	if any(valid_mask)
		[~, best_idx] = min(pop_costs(valid_mask));
		valid_indices = find(valid_mask);
		update_idx = valid_indices(best_idx);
		BestSol_pos = pop_positions(update_idx, :)';
		BestSol_cost = pop_costs(update_idx);
	end

	% Sort population by cost
	[pop_costs, idxSo] = sort(pop_costs);
	pop_steps = pop_steps(idxSo, :);
	pop_positions = pop_positions(idxSo, :);
	
	BestCost_curr = BestSol_cost;
	opts = ztox(BestSol_pos, IncomeG_low, IncomeG_high, a, b, n);

	% Fill pre-allocated S array with all sorted steps
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

	% Compute weighted average step for mean update - use pre-allocated storage
	weighted_steps_storage = pop_steps(1:mu, :) .* w';
	M_step_next = sum(weighted_steps_storage, 1)';
	
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

	% Update current best cost after all potential updates
	BestCost_curr = BestSol_cost;

	% Other termination checks
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

	% CMA-ES parameter updates: evolution paths and step size
	% Solve A_curr' * x = M_step for triangular system
	A_inv_M_step = A_curr' \ M_step;
	ps_next = (1-cs)*ps_curr + sqrt(cs*(2-cs)*mu_eff) * A_inv_M_step';
	ps_next = ps_next(:)';
	sigma_next = sigma_curr * exp(cs/ds*(norm(ps_next)/ENN - 1));

	if norm(ps_next) / sqrt(1-(1-cs)^(2*(g+1))) < hth
		hs = 1;
	else
		hs = 0;
	end
	delta_c = (1 - hs)*cc*(2 - cc);
	pc_next = (1-cc)*pc_curr + hs*sqrt(cc*(2-cc)*mu_eff)*M_step(:)';
	pc_next = pc_next(:)';

	% CMA-ES covariance matrix update using Cholesky rank-one updates
	
	% Start with scaled current factor
	scale_factor = sqrt(1 - c1 - cmu);
	A_next = scale_factor * A_curr;
	
	% Rank-one update: c1 * pc_next * pc_next'
	if c1 > 0
    A_next = choleskyRankOneUpdate(A_next, sqrt(c1) * pc_next(:), 1);
end

% Rank-mu update - use optimized BLAS operations (Optimization #4)
if cmu > 0
    weighted_steps_mu = S(:, 1:mu) .* sqrt(cmu * w);
    for i = 1:mu
        if abs(w(i)) > 1e-12
            A_next = choleskyRankOneUpdate(A_next, weighted_steps_mu(:,i), 1);
        end
    end
end

	% Stall detection and restart
	if g > 1
		if BestCost_curr >= BestCost_prev
			stall_counter = stall_counter + 1;
		else
			stall_counter = 0;
		end
	end

	if stall_counter >= 50
		innerstep = 0.5*innerstep;
		sigma_next = innerstep;
		% Mix with identity matrix and recompute Cholesky factor
		% Compute C_curr only when needed (Optimization #5)
		C_curr = A_curr * A_curr';
		C_restart = 0.5 * C_curr + 0.5 * eye(nVar);
		A_next = chol(C_restart, 'lower');
		ps_next = zeros(1, nVar);
		pc_next = zeros(1, nVar);
		M_pos = BestSol_pos;
		stall_counter = 0;
		if wtc==1
			disp(['Reducing innerstep to ' num2str(innerstep) ' and restarting']);
		end
	end
	
	% Update state variables
	ps_curr = ps_next;
	pc_curr = pc_next;
	A_curr = A_next;
	sigma_curr = sigma_next;
	BestCost_prev = BestCost_curr;
end

if WB==1 && ishandle(hWaitbar), close(hWaitbar); end

opts = ztox(BestSol_pos, IncomeG_low, IncomeG_high, a, b, n);
opts = roundToA(opts, a);

end


%function loss = CostFunctionWrapper(z, a, b, n, G, IncomeG_low, IncomeG_high, dissaving, dis, u, P, tol, f_low, f_high, dsp)
	% Wrapper for Lossfunction to compute cost from parameter vector z
%	opts = transform(z, IncomeG_low, IncomeG_high, a, b, n);
%	[loss, ~] = Lossfunction(a, b, n, opts, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);
%end

function x = ztox(z, IncomeG_low, IncomeG_high, a, b, n)
   % Transform parameters to almost monotone savings policies (almost is so
   % that flat segments are interior)
   c_low = cumsum(z(1:n+1));                    % Cumulative scores for low income state
   c_high = cumsum(z(n+2:2*(n+1)));            % Cumulative scores for high income state
   
   f_low = 1 ./ (1 + exp(-c_low));             % Logistic transform to [0,1]
   f_high = 1 ./ (1 + exp(-c_high));           % Logistic transform to [0,1]
   
   x = [a + (min(IncomeG_low, b) - a) .* f_low;    % Scale to feasible savings range
        a + (min(IncomeG_high, b) - a) .* f_high];  % Scale to feasible savings range
end

function z = xtoz(x, IncomeG_low, IncomeG_high, a, b, n)
   % Inverse transform (inverse of transform)
   eps = 1e-10;
   
   % Process low income state
   upper_low = min(IncomeG_low, b);
   denom_low = upper_low - a;
   f_low = (x(1:n+1) - a) ./ denom_low;
   f_low = max(eps, min(1-eps, f_low)); % Clamp to avoid logit singularities
   c_low = log(f_low ./ (1 - f_low));   % Logit transform
   
   % Process high income state
   upper_high = min(IncomeG_high, b);
   denom_high = upper_high - a;
   f_high = (x(n+2:2*(n+1)) - a) ./ denom_high;
   f_high = max(eps, min(1-eps, f_high)); % Clamp to avoid logit singularities
   c_high = log(f_high ./ (1 - f_high));   % Logit transform
   
   % Recover increments
   z = [[c_low(1); diff(c_low)]; 
        [c_high(1); diff(c_high)]];
end

function vals = roundToA(vals, a)
	% Correcting floating point imprecision at borrowing limit
	closeToA = abs(vals - a) < 1e-8;
	vals(closeToA) = a;
end

function grad_estimate = estimateGradient(center_pos, center_cost, positions, costs, sigma, lambda)
	% Estimate pseudo-gradient from CMA-ES function evaluations
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
	
	% Normalize by number of improving samples for stability
	n_improving = sum(costs < center_cost & isreal(costs') & costs > 0);
	if n_improving > 0
		grad_estimate = grad_estimate / n_improving;
	end
end

function A_new = choleskyRankOneUpdate(A, v, sign_update)
	% Rank-one Cholesky update: A_new such that A_new*A_new' = A*A' + sign_update*v*v'
	% Uses Givens rotations for numerical stability
	% A: lower triangular Cholesky factor (n×n)
	% v: update vector (n×1)  
	% sign_update: +1 for update, -1 for downdate
	
	n = size(A, 1);
	A_new = A;
	v_work = v(:);
	
	if sign_update > 0
		% Rank-one update
		for k = 1:n
			a_kk = A_new(k,k);
			v_k = v_work(k);
			
			if abs(v_k) < 1e-15
				continue;
			end
			
			r = sqrt(a_kk^2 + v_k^2);
			c = a_kk / r;
			s = v_k / r;
			
			A_new(k,k) = r;
			
			if k < n
				A_k_rest = A_new(k+1:n, k);
				v_rest = v_work(k+1:n);
				
				A_new(k+1:n, k) = c * A_k_rest + s * v_rest;
				v_work(k+1:n) = c * v_rest - s * A_k_rest;
			end
		end
		
	else
		% Rank-one downdate
		for k = 1:n
			a_kk = A_new(k,k);
			v_k = v_work(k);
			
			if abs(v_k) < 1e-15
				continue;
			end
			
			r_squared = a_kk^2 - v_k^2;
			if r_squared <= 1e-15
				r = 1e-8;
			else
				r = sqrt(r_squared);
			end
			
			c = a_kk / sqrt(a_kk^2 + r^2 - v_k^2);
			s = -v_k / sqrt(a_kk^2 + r^2 - v_k^2);
			
			A_new(k,k) = r;
			
			if k < n
				A_k_rest = A_new(k+1:n, k);
				v_rest = v_work(k+1:n);
				
				A_new(k+1:n, k) = c * A_k_rest + s * v_rest;
				v_work(k+1:n) = c * v_rest - s * A_k_rest;
			end
		end
	end
end