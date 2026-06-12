function [loss, br] = Lossfunction(a, b, n, x, dissaving, dis, u, p, rep, tol, f, offset, dsp, sm)
% Compute loss and best response for dynamic decision problem
%
% Inputs:
%   a, b       - lower and upper bounds of asset space
%   n          - grid resolution (number of subintervals)
%   x          - current policy vector
%   dissaving  - asset constraint flag (1=asset-constrained, 0=income-constrained)
%   dis        - discount sequence [delta*beta, delta]
%   u          - utility function handle
%   p          - probability parameter
%   rep        - return type (1=best response, 2=value function)
%   tol        - convergence tolerance for value iteration
%   f          - income function handle
%   offset     - asset transformation function handle
%   dsp        - grid density multiplier for optimization
%   sm         - threshold for zeroing losses below this level
%
% Outputs:
%   loss       - aggregate loss
%   br         - best response policy or value function

    directsearch = 1;                       % use grid search optimization

    %% Basic parameters
    delta = dis(2)/dis(1);                  % discount factor
    beta = dis(1)/delta;                    % preference parameter

    %% Value function computation
    [VFt, grid] = valf(a, b, n, u, p, delta, tol, x, f, offset);
    VF = @(z) interp1(grid, VFt, real(z), 'linear', 'extrap');

    %% Grid setup
    G = linspace(a, b, n+1)';               % asset grid for states
    Income = f(G);                          % pre-compute income values

    %% Best response computation
    if directsearch == 1
        % Grid search optimization
        numGrid = (n+1)*dsp;                % optimization grid size
        t = linspace(0, 1, numGrid);        % grid parameter
        
        if dissaving == 1
            savingGrid = a + t .* (G - a);  % asset-constrained choice set
        else
            savingGrid = a + t .* (Income - a); % income-constrained choice set
        end
        
        consVec = max(Income - savingGrid, 1e-10); % consumption choices
        
        % Compute objective function values
        VF_low = VF(savingGrid);            % value at saving choice
        offsetSavingGrid = offset(savingGrid); % transformed saving choice
        VF_high = VF(offsetSavingGrid);     % value at transformed choice
        
        objVals = u(consVec) + beta*delta*(p*VF_low + (1-p)*VF_high);
        objVals(isnan(objVals)) = -Inf;     % handle numerical issues
        
        % Find optimal choices with tie-breaking
        [maxValPerRow, ~] = max(objVals, [], 2);
        mask = abs(objVals - maxValPerRow) < 1e-14;
        [numRows, numCols] = size(objVals);
        columnIndices = repmat(1:numCols, numRows, 1);
        tieIndices = max(mask .* columnIndices, [], 2);
        
        % Handle invalid optimization results
        validIndices = tieIndices > 0 & tieIndices <= numCols;
        
        if ~all(validIndices)
            tieIndices(~validIndices) = 1;  % default to first choice
            
            for i = find(~validIndices)'
                validVals = objVals(i, ~isnan(objVals(i,:)));
                if ~isempty(validVals)
                    maxValPerRow(i) = max(validVals);
                else
                    maxValPerRow(i) = -Inf;
                end
            end
        end
        
        linearIndices = sub2ind(size(savingGrid), (1:numRows)', tieIndices);
        ORresp = savingGrid(linearIndices);  % optimal responses
        bestObj = maxValPerRow;              % optimal objective values
        
        % Handle boundary cases at a=0
        if a == 0
            infRows = isinf(bestObj) & bestObj < 0;
            if any(infRows)
                ORresp(infRows) = a;
                bestObj(infRows) = -1e10;   % finite penalty value
            end
        end
        
    else
        % Optimization using fminbnd
        ORresp = zeros(size(G));            % optimal responses
        bestObj = zeros(size(G));           % optimal objective values
        
        for i = 1:length(G)
            if dissaving == 1
                lb = a;                     % lower bound
                ub = G(i);                  % upper bound (asset-constrained)
            else
                lb = a;                     % lower bound
                ub = Income(i);             % upper bound (income-constrained)
            end
            
            income_i = Income(i);           % current income
            
            % Skip problematic cases
            if income_i <= 1e-8
                ORresp(i) = a;
                bestObj(i) = -1e10;
                continue;
            end
            
            % Optimize objective function
            obj = @(s) obj_function(s, income_i, u, beta, delta, p, VF, offset);
            [s_opt, fval] = fminbnd(obj, lb, ub);
            ORresp(i) = s_opt;
            bestObj(i) = -fval;
        end
    end

    %% Loss computation
    t_values = x(:);                        % current policy choices
    offset_t_values = offset(t_values);     % transformed choices
    VAL_t = beta*delta*(p*VF(t_values) + (1-p)*VF(offset_t_values));
    u_cons = u(max(Income - t_values, 1e-10)); % utility from current policy
    LOSSt = bestObj - u_cons - VAL_t;       % losses

    % Handle numerical issues
    nonRealMask = ~isreal(ORresp);          % non-real responses
    LOSSt(nonRealMask) = 1e8;
    negativeMask = LOSSt < 0;               % negative losses
    LOSSt(negativeMask) = 0;
    ORresp(negativeMask) = t_values(negativeMask);

    % Handle boundary case at a=0
    if a == 0
        if G(1) == 0 && (f(0) <= 1e-8 || isinf(LOSSt(1)))
            LOSSt(1) = 0;                   % zero loss contribution
            ORresp(1) = t_values(1);        % preserve original policy
        end
    end

    % Zero losses below threshold sm
    if nargin >= 14 && ~isempty(sm) && sm > a
        below_sm_mask = G <= sm;
        LOSSt(below_sm_mask) = 0;           % zero losses below threshold
    end

    %% Output assignment
    loss = sum(LOSSt);                      % aggregate loss
    br = ORresp;                            % best response policy
    if rep == 2
        br = VF;                            % return value function instead
    end
end

function val = obj_function(s, income, u, beta, delta, p, VF, offset)
    % Objective function for optimization
    cons = max(income - s, 1e-10);          % consumption level
    offset_s = offset(s);                   % transformed saving choice
    val = -(u(cons) + beta*delta*(p*VF(s) + (1-p)*VF(offset_s)));
end

function [W_out, gridFine] = valf(a, b, n, u, p, delta, tol, s, f, offset)
    % Value function iteration with adaptive grid refinement
    
    iterlimit = 100;                        % iteration limit
    gridFine = linspace(a, b, n+1);         % initial grid
    maxPoints = 2 * n;                      % maximum grid points
    
    if tol == 1e-20
        iterlimit = 10000;                  % high precision iteration limit
        gridFine = linspace(a, b, 5*(n+1)); % high precision grid
    end
    
    % Persistent storage for value function
    persistent W_old_persistent grid_persistent
    if isempty(W_old_persistent) || length(grid_persistent) ~= length(gridFine) || any(grid_persistent ~= gridFine)
        W_old = u(max(0.1, f(a)-a)) * ones(size(gridFine)); % initial guess
    else
        W_old = W_old_persistent;           % use stored value function
    end

    %% Value function iteration
    iter = 0;                               % iteration counter
    dist = Inf;                             % convergence distance
    originalGrid = linspace(a, b, n+1);    % policy grid
    fGridFine = f(gridFine);                % pre-compute income values
    
    while dist > tol && iter < iterlimit
        % Interpolate policy function
        gVals = interp1(originalGrid, s, gridFine, 'linear', 'extrap');
        gVals = max(a, gVals);              % enforce lower bound
        above_b = gridFine > b;             % extrapolation region
        
        if any(above_b)
            num_segments = min(10, length(s)-1); % segments for slope estimation
            segment_slopes = diff(s(end-num_segments:end)) ./ diff(originalGrid(end-num_segments:end));
            avg_slope = mean(segment_slopes); % average slope
            extrap_vals = s(end) + avg_slope*(gridFine(above_b) - originalGrid(end));
            gVals(above_b) = extrap_vals;   % extrapolated policy values
        end
        
        % Compute value function update
        cons = max(fGridFine - gVals, 1e-10); % consumption levels
        vLow = interp1(gridFine, W_old, gVals, 'linear', 'extrap');
        offsetGVals = offset(gVals);        % transformed policy values
        vHigh = interp1(gridFine, W_old, offsetGVals, 'linear', 'extrap');
        
        W_new = 0.5 * W_old + 0.5*(u(cons) + delta*(p*vLow + (1-p)*vHigh));
        dist = max(abs(W_new - W_old));     % convergence metric
        W_old = W_new;                      % update value function
        
        % Adaptive grid refinement
        if mod(iter, 5) == 0 && length(gridFine) < maxPoints
            oldGrid = gridFine;             % store current grid
            gradients = abs(diff(W_old));   % value function gradients
            avg_gradient = mean(gradients(end-4:end)); % recent average gradient
            refineThreshold = 0.2 * avg_gradient; % refinement threshold
            refineIdx = gradients > refineThreshold; % points to refine
            newPoints = (oldGrid(refineIdx) + oldGrid([false, refineIdx]))/2;
            
            if ~isempty(newPoints)
                gridFine = unique([oldGrid, newPoints(1:min(3, end))]); % add new points
                W_old = interp1(oldGrid, W_old, gridFine, 'linear', 'extrap');
                fGridFine = f(gridFine);    % update income values
            end
        end
        
        iter = iter + 1;                    % increment counter
    end

    W_out = W_old;                          % output value function
    
    % Update persistent storage
    W_old_persistent = W_out;
    grid_persistent = gridFine;
end