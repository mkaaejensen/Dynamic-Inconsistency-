function [opts] = algorithm_call_standard(a, b, n, MaxIt, innerstep, opts, RESTART, dissaving, dis, p, u, tol, f, offset,dsp)
% CMA-ES algorithm with policy iteration, using z = [x1, D1, ..., Dn].
% x is monotonic in [a,b] via ztox
%
% Rounding near a happens just before we save or return opts, including
% any time we also set dyncon = opts.
wtc=0;  % writes additional information to the consolue
WB = 1; % enable waitbar          
 

CostFunction = @(z) Lossfunction(a, b, n, ztox(z,a,b), ...
                                 dissaving, dis, u, p, 0, tol, f, offset,dsp);
nVar = n + 1;
VarSize = [1, nVar];
% Compute grid points G(j)
    G = linspace(a, b, n+1);
    
    % Compute variable-specific upper bounds (clamped by b)
    VarMax = min(f(G), b);
    VarMin = a;
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


x_old = interp1(linspace(a, b, length(opts)), opts, G, 'pchip');
x_old = max(min(x_old, b), a);

z_init = xtoz(x_old, a, b);


lambda = 8; 
mu = round(lambda/2);
w = log(mu + 0.5) - log(1:mu);
w = w / sum(w);
mu_eff = 1 / sum(w.^2);

sigma0 = innerstep; 
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

initialCov = eye(nVar);
C{1} = initialCov;
sigma{1} = sigma0;

empty_individual.Position = [];
empty_individual.Step = [];
empty_individual.Cost = [];

M(1).Position = z_init;
M(1).Step = zeros(VarSize);
M(1).Cost = CostFunction(M(1).Position);

if M(1).Cost < 0
    disp('Negative initial Ego Loss. Stopping.');
    MaxIt = 0;
end

BestSol = M(1);
BestCost = zeros(MaxIt, 1);
stall_counter = 0;
setflag = 0;

if WB == 1
    hWaitbar = waitbar(0, 'Initiating', ...
        'Name', 'Solving (cancel to skip)', 'CreateCancelBtn','delete(gcbf)');
    set(hWaitbar, 'Units','Pixels','Position',[50 500 380 100]);
end

brcount = 0;
for g = 1:MaxIt
    stepsAll = mvnrnd(zeros(1,nVar), C{g}, lambda);
    
    brcount = brcount + 1;
    if brcount == 2
        specialCandidateIndex = -1;%lambda;
    else
        specialCandidateIndex = -1;
    end

    % Preallocate cell array for candidate results
    pop_temp = cell(lambda, 1);
    
    % Parallel candidate evaluation loop
    parfor i = 1:lambda
        step_i = stepsAll(i,:);
        cand = M(g).Position + sigma{g} * step_i;
        cost_val = CostFunction(cand);
        pop_temp{i} = struct('Step', step_i, 'Position', cand, 'Cost', cost_val);
    end
    
    pop = [pop_temp{:}];
    
    if specialCandidateIndex == lambda
        brcount = 0;
    end

    for i = 1:lambda
        if pop(i).Cost < BestSol.Cost && isreal(pop(i).Cost) && (pop(i).Cost > 0)
            BestSol = pop(i);
        end
    end

    Costs = [pop.Cost];
    [~, idxSo] = sort(Costs);
    pop = pop(idxSo);
    BestCost(g) = BestSol.Cost;

    z_best = BestSol.Position;
    opts = ztox(z_best,a,b);

    if g > 1
        %disp(['Iteration ' num2str(g-1) ' Ego Loss = ' num2str(BestCost(g-1)/nVar)]); 
        if BestCost(g)/nVar<1e-14
                disp(['Ego Loss equals zero within tolerance (' num2str(BestCost(g)/nVar) '). Finishing']);
                if WB==1, close(hWaitbar); end
                opts = ztox(z_best,a,b);
                opts = roundToA(opts, a);
                break;
        end
        
        if WB==1 && ~ishandle(hWaitbar)
            disp('Canceled by user');
            disp(['Iteration ' num2str(g) ' [final] Ego Loss = ' num2str(BestCost(g)/nVar)]);
            g = MaxIt;
            opts = roundToA(opts, a);
            break;
        elseif WB==1
            waitbar((g-1)/MaxIt, hWaitbar, ...
                ['Ego Loss = ' num2str(BestCost(g-1)/nVar) ...
                 ' (Iteration ' num2str(g-1) ' of ' num2str(MaxIt) ')']);
        end
    end

    if BestCost(g) < 0
        disp('Divergence: Cost < 0.');
        break;
    end
    if ~isreal(BestCost(g))
        disp('Divergence: Cost is imaginary.');
        break;
    end

    if g == MaxIt
        disp(['Iteration ' num2str(g) ' [final] Ego Loss = ' num2str(BestCost(g)/nVar)]);
        if WB==1, close(hWaitbar); end
        opts = ztox(z_best,a,b);
        opts = roundToA(opts, a);
        break;
    end

    M(g+1).Step = zeros(VarSize);
    for jj = 1:mu
        M(g+1).Step = M(g+1).Step + w(jj)*pop(jj).Step;
    end
    new_z = M(g).Position + sigma{g}*M(g+1).Step;

    testVal = CostFunction(new_z);
    if ~isreal(testVal)
        disp('Imaginary candidate encountered.');
    end
    if testVal > 0 && isreal(testVal)
        M(g+1).Position = new_z;
        M(g+1).Cost = testVal;
        if testVal < BestSol.Cost && testVal > 0 && isreal(testVal)
            BestSol = M(g+1);
        end
    else
        M(g+1).Position = M(g).Position;
        M(g+1).Cost = CostFunction(M(g).Position);
    end

    ps{g+1} = (1-cs)*ps{g} + sqrt(cs*(2-cs)*mu_eff)* (M(g+1).Step / chol(C{g})');
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

    [VV, EE] = eig(C{g+1});
    if any(diag(EE) < 0)
        EE = max(EE, 0);
        C{g+1} = VV*EE/VV;
    end

    if g > 1
        if BestCost(g) >= BestCost(g-1)
            stall_counter = stall_counter + 1;
        else
            stall_counter = 0;
        end
    end
    if stall_counter >= 20
        innerstep = innerstep/2;
        sigma{g+1} = innerstep;
        C{g+1} = eye(nVar);
        ps{g+1} = zeros(VarSize);
        pc{g+1} = zeros(VarSize);
        M(g+1).Position = BestSol.Position;
        stall_counter = 0;
        if wtc==1
        disp(['Reducing innerstep to ' num2str(innerstep) ' and restarting']);
        end
    end
end

z_best = BestSol.Position;
opts = ztox(z_best,a,b);
opts = roundToA(opts, a);

end

function x = ztox(z, a, b)
nVar = length(z);
n = nVar - 1;
x1 = a + (b - a)/(1 + exp(-z(1)));
if n < 1
    x = x1;
    return;
end
alphas = exp(z(2:end));
S = sum(alphas);
if S < 1e-14
    x_inc = zeros(1,n);
else
    Dsum = (b - x1)*(S/(S+1));
    x_inc = (Dsum/S)*alphas;
end
x = zeros(1,n+1);
x(1) = x1;
x(2:end) = x1 + cumsum(x_inc);
end

function z = xtoz(x, a, b)
    x = max(min(x, b), a);
    nVar = length(x);
    n = nVar - 1;
    z = zeros(1, nVar);
    
    % Vectorized first computation
    frac = (x(1) - a) / (b - a + 1e-14);
    frac = max(min(frac, 1 - 1e-14), 1e-14);
    z(1) = log(frac / (1 - frac));
    
    if n < 1
        return;
    end
    
    % Vectorized difference and filtering
    D = max(diff(x), 0);
    Dtot = sum(D);
    leftover = b - x(1);
    
    if leftover < 1e-14
        z(2:end) = -10;
        return;
    end
    
    % Vectorized S computation
    if leftover <= Dtot + 1e-12
        S = 1e3;
    else
        S = max(1e-12, Dtot / (leftover - Dtot));
    end
    
    % Fully vectorized alpha and log computation
    alphas = max((S + 1) * D / leftover, 1e-14);
    z(2:end) = log(alphas);
end

function vals = roundToA(vals, a)
closeToA = abs(vals - a) < 1e-2;
vals(closeToA) = a;
end