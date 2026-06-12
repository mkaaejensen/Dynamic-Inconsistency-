function [loss, br] = Lossfunction(a, b, n, x, dissaving, dis, u, p, rep, tol, f, offset,dsp,sm)

directsearch = 1; % default is a direct search over discretized saving space. Set = 0, to instead use fminbnd. This is only useful for debugging as direct search if both faster and more precise.

% Basic parameters
delta = dis(2)/dis(1);
beta = dis(1)/delta;

% Compute value function given second-mover policy x
[VFt, grid] = valf(a, b, n, u, p, delta, tol, x, f, offset);
VF = @(z) interp1(grid, VFt, real(z), 'linear', 'extrap');

% Common variables
G = linspace(a, b, n+1)';          % asset grid for states
Income = f(G);

if directsearch == 1
    % --- Grid search method ---
    numGrid = (n+1)*dsp;
    t = linspace(0, 1, numGrid);
    if dissaving == 1
        savingGrid = a + t .* (G - a);
    else
        savingGrid = a + t .* (Income - a);
    end
    consVec = max(Income - savingGrid, 1e-10);
    
    % Pre-compute offset 
    VF_low = VF(savingGrid);
    offsetSavingGrid = offset(savingGrid);
    VF_high = VF(offsetSavingGrid);
    
    objVals = u(consVec) + beta*delta*(p*VF_low + (1-p)*VF_high);
    objVals(isnan(objVals)) = -Inf;
    
    [maxValPerRow, ~] = max(objVals, [], 2);
    mask = abs(objVals - maxValPerRow) < 1e-14;
    [numRows, numCols] = size(objVals);
    columnIndices = repmat(1:numCols, numRows, 1);
    tieIndices = max(mask .* columnIndices, [], 2);
    
    % Ensure tieIndices are within valid range
    validIndices = tieIndices > 0 & tieIndices <= numCols;
    
    % Handle borderline case where no valid maximum was found (e.g., all -Inf)
    if ~all(validIndices)
        % For rows with invalid indices, set to first column as default
        tieIndices(~validIndices) = 1;
        
        
        for i = find(~validIndices)'
            % Find the first non-NaN value, or use -Inf if all are NaN
            validVals = objVals(i, ~isnan(objVals(i,:)));
            if ~isempty(validVals)
                maxValPerRow(i) = max(validVals);
            else
                maxValPerRow(i) = -Inf;
            end
        end
    end
    
    linearIndices = sub2ind(size(savingGrid), (1:numRows)', tieIndices);
    ORresp = savingGrid(linearIndices);
    bestObj = maxValPerRow;
    
    % For cases where a=0 and all utilities are -Inf, substitute with appropriate values
    if a == 0
        infRows = isinf(bestObj) & bestObj < 0;
        if any(infRows)
            % Use minimum allowed saving for these cases
            ORresp(infRows) = a;
            % Set a finite but very low value for the objective
            bestObj(infRows) = -1e10;
        end
    end
else
    % --- fminbnd method ---
    ORresp = zeros(size(G));
    bestObj = zeros(size(G));
    for i = 1:length(G)
        if dissaving == 1
            lb = a;
            ub = G(i);
        else
            lb = a;
            ub = Income(i);
        end
                
        income_i = Income(i);
        
        % Skip optimization if income is too low (f(0)=0 case)
        if income_i <= 1e-8
            ORresp(i) = a;
            bestObj(i) = -1e10;  % Use a very negative but finite value
            continue;
        end
        
        % Optimize the objective function
        obj = @(s) obj_function(s, income_i, u, beta, delta, p, VF, offset);
        
        [s_opt, fval] = fminbnd(obj, lb, ub);
        ORresp(i) = s_opt;
        bestObj(i) = -fval;
    end
end

% Compute loss for each state based on first-mover policy x
t_values = x(:);

% Compute offset values
offset_t_values = offset(t_values);
VAL_t = beta*delta*( p*VF(t_values) + (1-p)*VF(offset_t_values) );
u_cons = u(max(Income - t_values, 1e-10));  % Add max to avoid -Inf
ORACLEt = bestObj - u_cons - VAL_t;

% Adjust for non-real or negative losses
nonRealMask = ~isreal(ORresp);
ORACLEt(nonRealMask) = 1e8;
negativeMask = ORACLEt < 0;
ORACLEt(negativeMask) = 0;
ORresp(negativeMask) = t_values(negativeMask);

% Handle special case for a=0 (especially if f(0)=0)
if a == 0
    % For the first grid point, if it's problematic, ignore its contribution
    if G(1) == 0 && (f(0) <= 1e-8 || isinf(ORACLEt(1)))
        ORACLEt(1) = 0;  % Don't add to loss
        ORresp(1) = t_values(1);  %# Keep the original policy
    end
end
% Zero out losses for states below sm (only for low state)
if sm > a
    below_sm_mask = G <= sm;
    ORACLEt(below_sm_mask) = 0;  % Set loss to zero for low state below sm
    % Note: high state losses are NOT zeroed out (if targeted adjustment is
    % needed here)
end
loss = sum(ORACLEt);
br = ORresp;
if rep==2
    br = VF;
end
end

% Auxiliary function for fminbnd
function val = obj_function(s, income, u, beta, delta, p, VF, offset)
    cons = max(income - s, 1e-10);
    offset_s = offset(s);
    val = -( u(cons) + beta*delta*(p*VF(s) + (1-p)*VF(offset_s)) );
end

function [W_out, gridFine] = valf(a, b, n, u, p, delta, tol, s, f, offset)
 iterlimit = 400;
 gridFine = linspace(a, b, 2*(n+1));
 maxPoints = 2 * n;
 
if tol==1e-20
 iterlimit = 10000;
 gridFine = linspace(a, b, 5*(n+1));
end
 
% Use persistent variables to store previous value function
persistent W_old_persistent grid_persistent
if isempty(W_old_persistent) || length(grid_persistent) ~= length(gridFine) || any(grid_persistent ~= gridFine)
 W_old = u(max(0.1,f(a)-a)) * ones(size(gridFine));
else
 W_old = W_old_persistent;
end

% Value function iteration
iter = 0;
dist = Inf;
originalGrid = linspace(a, b, n+1);
 
% Fine income grid
fGridFine = f(gridFine);
 
while dist > tol && iter < iterlimit
 % Policy function interpolation
 gVals = interp1(originalGrid, s, gridFine, 'linear', 'extrap');
 gVals = max(a, gVals);
 above_b = gridFine > b;
 
 if any(above_b)
   num_segments = min(10, length(s)-1);
   segment_slopes = diff(s(end-num_segments:end)) ./ diff(originalGrid(end-num_segments:end));
   avg_slope = mean(segment_slopes);
   extrap_vals = s(end) + avg_slope*(gridFine(above_b) - originalGrid(end));
   gVals(above_b) = extrap_vals;
 end
 
 % consumption and values
 cons = max(fGridFine - gVals, 1e-10);  % Added max to avoid negative consumption
 vLow = interp1(gridFine, W_old, gVals, 'linear', 'extrap');
 vHigh = interp1(gridFine, W_old, offset(gVals), 'linear', 'extrap');
 
 W_new = 0.5 * W_old + 0.5*(u(cons) + delta*(p*vLow + (1-p)*vHigh));
 dist = max(abs(W_new - W_old));
 W_old = W_new;
 
 if mod(iter, 5) == 0 && length(gridFine) < maxPoints
   oldGrid = gridFine;
   gradients = abs(diff(W_old));
   avg_gradient = mean(gradients(end-4:end));
   refineThreshold = 0.2 * avg_gradient;
   refineIdx = gradients > refineThreshold;
   newPoints = (oldGrid(refineIdx) + oldGrid([false, refineIdx]))/2;
   
   if ~isempty(newPoints)
     gridFine = unique([oldGrid, newPoints(1:min(3, end))]);
     W_old = interp1(oldGrid, W_old, gridFine, 'linear', 'extrap');
     
     % Update 
     fGridFine = f(gridFine);
   end
 end
 
 iter = iter + 1;
end

W_out = W_old;
 
% Update persistent storage
W_old_persistent = W_out;
grid_persistent = gridFine;
end