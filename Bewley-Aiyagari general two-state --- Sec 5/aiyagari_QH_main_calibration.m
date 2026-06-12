clear;
%clear all;
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
addpath(fullfile(pwd, 'figurescripts'));
addpath(fullfile(pwd, 'mainscripts'));

%% Initiation, storing and loading
loadoutput  = 1;                    % 0: initiates from horizontal policy; 1: load opts from output.mat with handle specified next
    Handle = 'aiyagari_QH_calb_new_ver';  
                                    % Appended to save and load filenames (saves
                                    % in "savedmatfiles" with a separate file called dyncon is saved when beta=1)
figonly     = 1;                    % = 1 draws figures from saved solution(s) [no algorithm call]

%% Preferences
delta        = 0.963149;      % Calibrated via wrapper
calmode      = 1;             % Compare DC and Markov policies  
    deltaDC  = 0.955;         % discount factor for DC case
    EGPfcal  = 1;             % = 0: loads dyncon.m computed with Loss algorithm for DC comparison, = 1: uses EGP to compute DC comparison
beta         = 0.8;           % 
risk_aver    = 2;             % Rate of risk aversion
                              % Aiyagari (1994) looks at the [1,5] range (e.g. https://mark-ponder.com/tutorials/heterogeneous-agents-aiyagari/)

%% Production  
A        = 1.00;              % Technology index
capshare = 0.36;              % Capital share
depr     = 0.08;              % Rate of depreciation    
k        =5.257;              % Capital-labor ratio

%% Bounds 
a           = 0;              % Lower bound (borrowing limit)
b           = 25;             % Upper bound for assets (set high when beta<1 because extension is not approximately linear in this case and so extrapolation is less good approximation than when beta = 1)
    bspecfig    = 8;          % Separate upper bound for Figure 1
power_param = 1.1;            % Power-spaced grid parameter:
                              % power_param > 1: more points near a
                              % power_param = 1: uniform grid

%% Labor productivity process
CV_target = 1-0.7551;         % from Aiyagari's calibration via Rouwenhorst method (Aiyagari's AR(1) process to a two-state chain
                              % AR(1): log(z_{t+1}) = ρ log(z_t) + σ (1-ρ^2)^(0.5) ε_t with ρ = 0.6, σ = 0.2) 
                              % Rouwenhorst method implemented in aiyagari_1994_full.m in folder /full
                              % Aiyagari model simulation 
rho = 0.6;
p = (1 + rho)/2;
P = [p, 1-p; 1-p, p];

%P = [0.8,0.2;0.1,0.9];       % optional overwrite for asymmetric specifications
                              % (lmin can be overwritten below)
% Stationary distribution
[V_eig, D_eig] = eig(P');
[~, idx] = min(abs(diag(D_eig)-1));
stat = V_eig(:, idx);
stat = stat / sum(stat);
p_low = stat(1);
p_high = stat(2);

% Solve for lmin to match CV
obj = @(lmin) ...
    sqrt(p_low*(lmin - 1)^2 + p_high*((1 - p_low*lmin)/p_high - 1)^2) - CV_target;
lmin_sol = fzero(obj, [0.01, 0.99]);

% Labor productivities (can be overwritten for direct specifications)
lmin = lmin_sol;                    % When P is symmetric, lmin = 1 - CV_target
lmax = (1 - p_low*lmin) / p_high;   % Normalization (implies average labor supply 1)

if rcond(P)<0.1
    disp(['This algorithm is not recommended for near-singular Markov processes as it flattens the search landscape and explodes the noise-to-signal ratio for CMA-ES.']);
    disp(['In cases like this first run the i.i.d. implementation, then import opts and run this script (disabling this warning). ']);
    return
end

% Prices (from production function)

F  = @(a) A*a^(capshare);           % Production function
DF = @(a) A*capshare*a^(capshare-1); 
                                    % MPK
r = DF(k)-depr;                     % Interest rate (R-1)
W = F(k) - DF(k)*k;                 % Wage rate
  
% Income at each labor state 
f_low = @(ass) (1+r)*ass + W*lmin;  % Income in low state
f_high = @(ass) (1+r)*ass + W*lmax; % Income in high state


                              


%% Output/figure settings
shortv      = 0;                    % = 1 prints mean and breaks after first figure (showing Markov policies only)
showval     = 1;                    % = 1 displays value function among figures                              
extrads     = 1;                    % full suite of inequality statistics (Teil, Atkinson, Palma, etc)
extratailm  = 0;                    % full tail analysis

%% Choice of algorithm
algorithm   = 'default';      % Options: 'default' or 'sep'  


%% Precision parameters
dissaving   = 1;              % = 1 restricts second-mover low Markov policy to weak dissaving policies
gridprec    = 20;             % Grid intervals per unit of assets
dsp         = 2;              % Second-mover direct search precision
innerstep   = 0.1 / sqrt(gridprec*(b-a)); 
                              % Initial step length for CMA-ES
tol = 1e-8;                   % Tolerance (value-function iteration)
tolalg= 3e-3;                 % Tolerance (loss function)
MaxIt = 6000;                  % Maximum iterations

if beta == 1                  % DC precision overwrites
    tolalg = 5e-7;
    dsp =10;
end

%% Internal settings (changes from now on might lead to errors)                              
n = round(gridprec*(b-a),0);   % Number of grid intervals (n+1 is number of grid points)

if figonly == 1
    RESTART = 0;
    cutoff = 300;
    MaxIt = 0;
    tol = 1e-16;              % Tolerance for value-function iteration (allows very precise figures)
end
    
  

if figonly == 1               % Concistency overwrite
    shortv = 0;
    loadoutput = 1;
end

showkaplan = 0;               % Shows standard Kaplan aiyagari simulation
shownaive  = 0;               % Shows naive policy in main plot (only low productivity)
bspecfigr   = b;              % Separate b for other figures (usually unnecessary to change)
bspecfigcv  = b; 
bspecfigcss = b;
 
% Utility function and discount sequence
if risk_aver == 1
    u = @(c) log(c);
else    
    u = @(c) (c.^(1-risk_aver)-1)./(1-risk_aver);
end

cutoff = 150;                 % Cutoff for finite-sum approximations 
dis = @(j) beta.*(delta.^j);
dis = dis(1:cutoff);

%% Compute grid and income values
t = linspace(0, 1, n+1);
G = a + (b-a)*(t.^power_param);        % Asset grid
G = G(:);

% Income at each asset level on the grid
IncomeG_low = f_low(G);                % Income in low state
IncomeG_high = f_high(G);              % Income in high state

if exist('r', 'var')
    disp(['Capital-output ratio  ' num2str(k/F(k))]);
    disp(['Interest rate: ' num2str(100*r) '%, Rstar ' num2str(100*((1-delta)/(beta*delta))) '%, Wage: ' num2str(W) ]);
    if r>(1-delta)/(beta*delta)
        disp('Warning: Interest rate is higher than Rstar which may result in unlimited accumulation');
    end
    disp(['Borrowing limit: ' num2str(a) ', Consumption at borrowing limit ' num2str(r*a+W*lmin) ', Natural borrowing limit: ' num2str(-((W*lmin)/(r)))]);
end

actual_CV = sqrt(p_low*(lmin - 1)^2 + (1 - p_low)*(lmax - 1)^2);
% Check symmetry of transition matrix
is_symmetric = all(abs(P(1,:) - P(2,end:-1:1)) < 1e-10); % robust check

if is_symmetric && exist('rho', 'var')
    disp(['Symmetric transition matrix with persistence (rho): ' num2str(rho) ', implied Markov transition matrix P = [' ...
        num2str(P(1,1)) ', ' num2str(P(1,2)) '; ' num2str(P(2,1)) ', ' num2str(P(2,2)) ']']);
else
    % Compute implied persistence
    rho_implied = P(1,1) + P(2,2) - 1;
    disp(['Asymmetric transition matrix, implied persistence (rho): ' num2str(rho_implied) ', transition matrix P = [' ...
        num2str(P(1,1)) ', ' num2str(P(1,2)) '; ' num2str(P(2,1)) ', ' num2str(P(2,2)) ']']);
end

disp(['Assumed cross-sectional coefficient of variation (CV): ' num2str(CV_target) ', implied productivity shocks (lmin/lmax): ' ...
    num2str(lmin) ' / ' num2str(lmax)]);
disp(['Stationary distribution over productivity states: [' num2str(stat(1)) ', ' num2str(stat(2)) ']']);




 

%% Retrieve dynamically consistent solution via EGP
if exist('r', 'var')
    disp('Computing saving functions in dynamically consistent (DC) case with EGP algorithm')
    if calmode == 1
        [dyncon_kaplan, grid_kaplan] = aiyagari(a, b, n, r, W, lmin, lmax, risk_aver, P, deltaDC, power_param);
    else 
        [dyncon_kaplan, grid_kaplan] = aiyagari(a, b, n, r, W, lmin, lmax, risk_aver, P, delta, power_param);
    end
end

%% Save and load handles and folder specifications
algorithm_handle = str2func(['algorithm_call_', algorithm]);
if isempty(Handle)
    optssave = 'output.mat';
    dynconsave = 'dyncon.mat';
else
    optssave = sprintf('output_%s.mat', Handle);
    dynconsave = sprintf('dyncon_%s.mat', Handle);
end
 
subfolder1 = 'savedmatfiles';
subfolder2 = sprintf('output_%s', Handle);

% Check if the subfolder structure exists, if not, create it
if ~exist(fullfile(subfolder1, subfolder2), 'dir')
    mkdir(fullfile(subfolder1, subfolder2));
end

% Define the full file path
optssave = fullfile(subfolder1, subfolder2, optssave);
dynconsave = fullfile(subfolder1, subfolder2, dynconsave);
 
%% Initiation 

if beta == 1 
    dyncon = a*ones(2*(n+1), 1);   
else
    if calmode == 1        
        if EGPfcal == 1 % Use the dyncon_kaplan solution directly
            dyncon = dyncon_kaplan;            
         elseif isfile(dynconsave)
            prev = load(dynconsave, 'dyncon', 'saved_power_param');
            nold = (length(prev.dyncon)/2) - 1;       % 2*(nold+1) = length(prev.dyncon)
            if prev.saved_power_param ~= power_param || nold ~= n
                told = linspace(0, 1, nold+1);
                Gold = a + (b-a) * (told.^(prev.saved_power_param));   
                tnew = linspace(0, 1, n+1);
                Gnew = a + (b-a) * (tnew.^power_param);
                dynconold_low  = prev.dyncon(1 : (nold+1));
                dynconold_high = prev.dyncon((nold+2) : (2*(nold+1)));
                dyncon_low  = interp1(Gold, dynconold_low,  Gnew, 'linear');
                dyncon_high = interp1(Gold, dynconold_high, Gnew, 'linear');
                dyncon = [dyncon_low dyncon_high];
            else
                dyncon = prev.dyncon;
end

        end
    end
end    

RESTART = 0;     % not in use (included for backward compatibility)
opts = a * ones(2*(n+1), 1); % initiation of policy iteration following next
%% Extended policy generation

if figonly ~= 1 %&& beta ~= 1 
    if loadoutput == 0
        disp('Generating initial guess and extension above b via supervised policy iteration');
    else
        disp('Generating sophisticated extension above upper bound for algorithm');
    end
    % Extended parameters    
    b_extended = 5*b;
    gridprec_extended = gridprec;
    dsp_extended = 1*dsp;
    opts = extended_policy(a, b, n, gridprec, power_param, b_extended, gridprec_extended, dsp_extended, dissaving, dis, u, P, tol, f_low, f_high, opts);
    
    if b_extended == b && gridprec_extended == gridprec && dsp_extended == dsp % Explanation: This allows us to, by not extending or refining, revert to naive (linear) extrapolation in vfi (but we still use the generated candidate). Useful for diagnostics.
        clear Lossfunction;
    end
    if loadoutput == 0
        
         
        [cost_interp, ~] = Lossfunction(a, b, n, opts, dissaving, dis, u, P, 0, tol, ...
                                        f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);
        loss_interp=cost_interp/(2*(n+1));
        disp(['Initial candidate loss: ', sprintf('%.4e', loss_interp)]);
    end
end

if loadoutput == 1 % overwriting [a,b] (but leaving extension in place for Lossfunction) with loaded policy
    if isfile(optssave)
        prev = load(optssave, 'opts', 'saved_power_param');
        nold=(length(prev.opts)/2)-1;       % 2*(nold+1) = length(prev.opts)
        if prev.saved_power_param ~= power_param || nold ~= n
            told = linspace(0, 1, nold+1);

            Gold = a + (b-a) * (told.^(prev.saved_power_param));   
            tnew = linspace(0, 1, n+1);
            Gnew = a + (b-a) * (tnew.^power_param);
            optsold_low  = prev.opts(1 : (nold+1));
            optsold_high = prev.opts((nold+2) : (2*(nold+1)));
            opts_low=interp1(Gold, optsold_low, Gnew, 'linear');
            opts_high=interp1(Gold, optsold_high, Gnew, 'linear');
            opts=[opts_low opts_high];
        else
            opts=prev.opts;
        end
        [cost_interp, ~] = Lossfunction(a, b, n, opts, dissaving, dis, u, P, 0, tol, ...
                                        f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);
        loss_interp=cost_interp/(2*(n+1));
        disp(['Initial candidate loss: ', sprintf('%.4e', loss_interp)]);
    
      
    else
        disp(['No ' optssave ' file found. Initiating from initial guess']);
        %opts = a * ones(2*(n+1), 1); % no overwrite or initial condition
        %necessay as extended policy can be used
    end
end
  
 

%% Executing algorithm
tic;
if MaxIt > 0
    if beta == 1
        disp(['Initiating ' algorithm ' CMA-ES and computing Markov policy with beta = ' num2str(beta) ' (dynamic consistency) and delta = ' num2str(delta,'%.7f')]);
    else
        disp(['Initiating ' algorithm ' CMA-ES and computing Markov policy with beta = ' num2str(beta) ' (dynamic inconsistency) and delta = ' num2str(delta,'%.7f')]);
    end
    
    % Pass grid and income values to algorithm_call
    opts = algorithm_handle(a, b, n, MaxIt, innerstep, opts, RESTART, dissaving, dis, P, u, tol, f_low, f_high, dsp, power_param, tolalg, G, IncomeG_low, IncomeG_high);
    

    time = round(toc);
    disp(['...Total Loss Algorithm run time: ' num2str(time) ' seconds.']);
    % clearing persistent variables from Lossfunction
    clear cached_W_out cached_initialized VF_low_interp VF_high_interp tFine_power numGrid columnIndices_template vf_low_interp_vfi vf_high_interp_vfi last_grid_size;

end

%% Saving results
if beta == 1 && figonly == 0
    dyncon = opts;
    saved_power_param = power_param;  % Save the power_param used for this grid
    save(dynconsave, 'dyncon','saved_power_param');
    dyncon = dyncon_kaplan;
elseif figonly == 0
    saved_power_param = power_param;  % Save the power_param used for this grid
    save(optssave, 'opts','saved_power_param');
end

% second-mover strategy   
[~, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);

%disp('Computing naive solution')
[~, naive] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, P, 0, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);

% computing value functions
[~, VF]   = Lossfunction(a, b, n, opts, dissaving, dis, u, P, 2, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high); 
[~, VFDC] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, P, 2, tol, f_low, f_high, dsp, G, IncomeG_low, IncomeG_high);


disp('Computing stationary distribution and simulating asset path');
if exist('r', 'var')
    [MA, asim, lsim, Tsim] = stationary_dist(a, b, n, r, W, lmin, lmax, P, opts, f_low, f_high, power_param);
    [MAK, asimK, lsimK, TsimK] = stationary_dist(a, b, n, r, W, lmin, lmax, P, dyncon, f_low, f_high, power_param);
    %[MAKap, asimKap, lsimKap, TsimKap] = stationary_dist(a, b, n, r, W, lmin, lmax, P, dyncon_kaplan, f_low, f_high, power_param);
end
 
if shownaive==1
        [MAN, asimN, ysimN, TsimN] = stationary_dist(a, b, n, r, W, lmin, lmax, P, naive, f_low, f_high, power_param);
    else
        [MAN, asimN, ysimN, TsimN]= deal([]);
end
   

%% Plotting results
% Save workspace in structure for easy export
vars = who;
S = cell2struct(cellfun(@(v) evalin('caller', v), vars, 'UniformOutput', false), vars, 1);

figures_random(S);

clear S;

function opts = extended_policy(a, b, n, gridprec, power_param, b_extended, gridprec_extended, dsp_extended, dissaving, dis, u, P, tol, f_low, f_high, opts_initial)
 diagnostics = 1; % set to 1 to see values and iteration path
% Calculate extended range grid size
 n_extended = round(gridprec_extended*(b-a), 0);
% Create extended range grid
 t_extended = linspace(0, 1, n_extended+1);
 G_extended = a + (b_extended-a)*(t_extended.^power_param);
 G_extended = G_extended(:);
% Extended range income values
 Income_extended_low = f_low(G_extended);
 Income_extended_high = f_high(G_extended);
% Initialize with flat policy at borrowing limit on extended range
 opts_extended = a * ones(2*(n_extended+1), 1);
% Storage for loss tracking and policy evolution
 old_costs = [];
 new_costs = [];
 policy_history_low = [];
% Supervised policy iteration
for m = 1:100
 [old_cost, Br] = Lossfunction(a, b_extended, n_extended, opts_extended, dissaving, dis, u, P, 3, tol, f_low, f_high, dsp_extended, G_extended, Income_extended_low, Income_extended_high);
 opts_extended2 = Br(:);
 [new_cost, ~] = Lossfunction(a, b_extended, n_extended, opts_extended2, dissaving, dis, u, P, 3, tol, f_low, f_high, dsp_extended, G_extended, Income_extended_low, Income_extended_high);
% Store loss values
 old_costs(m) = old_cost;
 new_costs(m) = new_cost;
% Store low productivity policy from current iteration
 policy_history_low(:,m) = opts_extended(1:(n_extended+1));
if new_cost < old_cost
 opts_extended = opts_extended2;
if new_cost < 1e-7
break;
end
else
break;
end
end
% Make sure Lossfunction's cache hols the optimal policy (and not the
% one when the loss increased)
 [~, ~] = Lossfunction(a, b_extended, n_extended, opts_extended, dissaving, dis, u, P, 3, tol, f_low, f_high, dsp_extended, G_extended, Income_extended_low, Income_extended_high);
% Extract policies for low and high productivity states
 opts_extended_low = opts_extended(1:(n_extended+1));
 opts_extended_high = opts_extended((n_extended+2):(2*(n_extended+1)));
% Create original grid for interpolation
 t = linspace(0, 1, n+1);
 G = a + (b-a)*(t.^power_param);
 G = G(:);
% Interpolate back to original grid
 opts_low = interp1(G_extended, opts_extended_low, G, 'linear', 'extrap');
 opts_high = interp1(G_extended, opts_extended_high, G, 'linear', 'extrap');
% Ensure column vectors and combine
 opts_low = opts_low(:);
 opts_high = opts_high(:);
 opts = [opts_low; opts_high];
% Display diagnostics
if diagnostics == 1
 costpath = [old_costs(1), new_costs];
 Loss = costpath./(2*(n_extended+1));
 fmt = arrayfun(@(x) sprintf('%.4e', x), Loss, 'UniformOutput', false);
 fmt{end} = ['(', fmt{end}, ')'];
 disp(['Loss path: ', strjoin(fmt, ' ')]);
% Plot evolution of low productivity saving functions
 figure(99); clf;
 set(gcf, 'Name', 'Supervised Policy Iteration', 'NumberTitle', 'off');
 colors = lines(size(policy_history_low,2));
 hold on;
% Add 45 degree line
 plot(G_extended, G_extended, 'k--', 'LineWidth', 2, 'DisplayName', '45° Line');
% Add vertical line at b
 plot([b b], ylim, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', 'Original Grid Boundary (b)');
for i = 1:size(policy_history_low,2)
    if i == size(policy_history_low,2)
        % Final iteration - solid line, thicker, green
        line_style = '-';
        line_width = 2.5;
        line_color = [0 0.5 0]; % Green
    else
        % All other iterations - dashed line
        line_style = '--';
        line_width = 1.5;
        line_color = colors(i,:);
    end
    
    if i == 1
        label_name = 'Initialization';
    elseif i == size(policy_history_low,2)
        label_name = sprintf('Iteration %d (final candidate)', i-1);
    else
        label_name = sprintf('Iteration %d', i-1);
    end
    
 plot(G_extended, policy_history_low(:,i), line_style, 'Color', line_color, 'LineWidth', line_width, ...
'DisplayName', label_name);
end
 xlabel('Assets');
 ylabel('Savings (Low Productivity)');
 title('Supervised Policy Iteration Path (Low Productivity Saving Function)');
 legend('show', 'Location', 'best');
 grid on;
 hold off;
 drawnow
% Check if final candidate exceeds 45 degree line
 final_policy = policy_history_low(:,end);
 if any(final_policy(2:end) >= G_extended(2:end))
     disp('Final candidate policy function exceeds 45° line at some point(s). Typically this happens if dsp_extended is too low');
 end
end
end