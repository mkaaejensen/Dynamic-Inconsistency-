clear;
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
addpath(fullfile(pwd, 'figurescripts'));
addpath(fullfile(pwd, 'mainscripts'));

%% Model Parameters
% Preferences
delta     = 0.95;    % Discount factor
beta      = 0.93;    % Present-bias parameter
risk_aver = 1;       % Rate of risk aversion

%% Initiation
loadoutput  = 1;     % 0: horizontal policy; 1: load from saved file
Handle = 'stoch_opt_growth';

%% Algorithm settings
algorithm = 'default';
algorithm_handle = str2func(['algorithm_call_', algorithm]);

%% Save and load folder structure
subfolder1 = 'savedmatfiles';
subfolder2 = sprintf('output_%s', Handle);

if ~exist(fullfile(subfolder1, subfolder2), 'dir')
    mkdir(fullfile(subfolder1, subfolder2));
end

if isempty(Handle)
    optssave = fullfile(subfolder1, 'output.mat');
    dynconsave = fullfile(subfolder1, 'dyncon.mat');
else
    optssave = fullfile(subfolder1, subfolder2, sprintf('output_%s.mat', Handle));
    dynconsave = fullfile(subfolder1, subfolder2, sprintf('dyncon_%s.mat', Handle));
end

%% Numerical Parameters
% Grid parameters
a = 0;                                    % Lower bound
b = 20;                                   % Upper bound
bspecfig = 12;                            % Upper bound for figures
gridprec = 30;                            % Grid precision
dsp =4;                                   % Direct search precision
n = round(gridprec*(b-a),0);              % Number of grid points

% Algorithm parameters
tol = 1e-6;                              % Convergence tolerance
innerstep = 0.1 / sqrt(gridprec*(b-a));  % Initial step size
tola = 1e-6;                             % CMA-ES tolerance
MaxIt = 5000;                            % Maximum iterations

%% Labor Endowment Process
lmin = 0.9;                    % Low productivity state
p = 0.5;                       % Probability of low state
lmax = (1-p*lmin)/(1-p);      % High productivity state (normalized mean = 1)

% Discount sequence parameters
cutoff = 149;                  % Cutoff for finite-sum approximations

%% Utility function and discount sequence
if risk_aver == 1
    u = @(c) log(c);
else    
    u = @(c) (c.^(1-risk_aver)-1)./(1-risk_aver);
end

dis = @(j) beta.*(delta.^j);
dis = dis(1:cutoff);

%% Production Technology
% Production parameters
A = 5;                                          % Technology level
capshare = 0.36;                                % Capital share

% State-dependent production functions
f_low = @(a) A .* lmin .* a.^capshare;         % Output with low productivity
f_high = @(a) A .* lmax .* a.^capshare;        % Output with high productivity

% Productivity offset function
% For any a, f_high(a) = f_low(offset(a)) - maintains equivalence between states
offset = @(a) a .* (lmax/lmin).^(1/capshare);  % Derived for power function specification
f = f_low;                                      % Base production function
disp(['Quasi-hyperbolic Brock-Mirman with beta = ' num2str(beta) ' and delta = ' num2str(delta)]);
disp(['Probability of low labor endowment: p = ' num2str(p) ' , low / high endowment ' num2str(lmin) ' / ' num2str(lmax)]);

%% Load or create dynamically consistent solution
if beta == 1 
    dyncon = a*ones(1, n+1);
else
    if isfile(dynconsave)
        load(dynconsave);
    else
        disp('Please run with beta=1 to generate dynamically consistent solution for comparison plots');
        return;
    end
    dyncon = interp1(linspace(a, b, length(dyncon)), dyncon, linspace(a, b, n+1), 'linear');
end

%% Initialize or load quasi-hyperbolic policy
if loadoutput == 1
    if exist(optssave, 'file') == 2
        load(optssave); 
        disp(['Initiating from opts saved in ' optssave])
    else
        disp(['No ' optssave ' file found. Please change loadoutput = 1 to loadoutput = 0']);
        return;
    end
else
    xx = linspace(a, b, n+1);
    opts = 0.2.*f(xx);
end

%% Run algorithm
if MaxIt > 0
    tic;
    disp(['Initiating ' algorithm ' CMA-ES and computing Markov policy'])
    opts = algorithm_handle(a, b, n, MaxIt, innerstep, opts, 1, 0, dis, p, u, tol, f, offset, dsp,tola);
    
    time = round(toc);
    disp(['...Total Ego Loss algorithm run time: ' num2str(time) ' seconds.']);
    
    % Save results
    if beta == 1
        dyncon = opts;
        save(dynconsave,'dyncon');
    else
        save(optssave,'opts');
    end
end

%% Compute naive policy
[~, naive] = Lossfunction(a, b, n, dyncon, 0, dis, u, p, 0, tol, f, offset, 4);
naive = transpose(naive);

%% Computing stationary distributions and statistics
disp(' ');
disp('=== COMPARISON SUMMARY ===');

[meanass_qh, meanout_qh, meancons_qh, meansavrate_qh] = stationary_dist(a, b, n, p, opts, f, offset, lmin, lmax); 
[meanass_dc, meanout_dc, meancons_dc, meansavrate_dc] = stationary_dist(a, b, n, p, dyncon, f, offset, lmin, lmax); 
[meanass_naive, meanout_naive, meancons_naive, meansavrate_naive] = stationary_dist(a, b, n, p, naive, f, offset, lmin, lmax);

disp(['Model                        Assets    Output  Consump. Sav.Rate']);
disp(['Quasi-hyperbolic         ' sprintf('%10.4f %8.4f %8.4f %8.4f', meanass_qh, meanout_qh, meancons_qh, meansavrate_qh)]);
disp(['Dynamically Consistent   ' sprintf('%10.4f %8.4f %8.4f %8.4f', meanass_dc, meanout_dc, meancons_dc, meansavrate_dc)]);
disp(['Naive                    ' sprintf('%10.4f %8.4f %8.4f %8.4f', meanass_naive, meanout_naive, meancons_naive, meansavrate_naive)]);
disp('All quantities refer to mean value of stationary distribution ');

%% Plotting results
vars = who;
S = cell2struct(cellfun(@(v) evalin('caller', v), vars, 'UniformOutput', false), vars, 1);
figure_stoch_og(S);
clear S;

function [meanass, meanout, meancons, meansavrate] = stationary_dist(a, b, n, p, policy, f, offset, lmin, lmax)

Nsim = 50000;
Tsim = 1000;
burn_in = 500;

b_adj = b;
test_pts = linspace(a, b, 100);
for i = 1:length(test_pts)
    if offset(test_pts(i)) > b
        b_adj = test_pts(i-1);
        break;
    end
end

agrid = linspace(a, b_adj, n+1)';
policy_low = @(x) interp1(agrid, policy, x, 'linear', 'extrap');
policy_high = @(x) policy_low(offset(x));

ygrid = [lmin; lmax];
ydist = [p; 1-p];

rng(2017);
yrand = rand(Nsim, Tsim);
asim = zeros(Nsim, Tsim);
asim(:,1) = 0.5 * ones(Nsim, 1);
yindsim = 1 + (yrand > ydist(1));
ysim = ygrid(yindsim);

for t = 1:(Tsim-1)
    low_mask = (yindsim(:,t) == 1);
    asim(low_mask, t+1) = policy_low(asim(low_mask, t));
    asim(~low_mask, t+1) = policy_high(asim(~low_mask, t));
end

t_keep = (burn_in+1):Tsim;
a_steady = asim(:, t_keep);
y_steady = ysim(:, t_keep);
yind_steady = yindsim(:, t_keep);

outsim = f(a_steady) .* (y_steady / lmin);
csim = zeros(size(a_steady));
csim(:,1:end-1) = outsim(:,1:end-1) - a_steady(:,2:end);

low_mask = (yind_steady(:,end) == 1);
next_assets = zeros(size(a_steady,1),1);
next_assets(low_mask) = policy_low(a_steady(low_mask, end));
next_assets(~low_mask) = policy_high(a_steady(~low_mask, end));
csim(:,end) = outsim(:,end) - next_assets;

meanass = mean(a_steady(:));
meanout = mean(outsim(:));
meancons = mean(csim(:));
meansavrate = 1 - meancons / meanout;

end