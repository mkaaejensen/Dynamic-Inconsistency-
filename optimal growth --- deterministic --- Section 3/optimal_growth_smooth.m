clear;
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
addpath(fullfile(pwd, 'mainscripts'));

%% Initiation, storing and loading
loadoutput  = 0;                    % 0: initiates from horizontal policy; 1: load opts from output.mat with handle specified next
    Handle = 'optimal_growthsmooth';  % Appended to save and load filenames (saves
                                    % in "savedmatfiles" with a separate file called dyncon is saved when beta=1)
figonly     = 0;                    % = 1 draws figures from saved solution(s) [no algorithm call]

%% Preferences
delta        = 0.96;                % Discount factor
calmode      = 0;                   % Calibration mode to compare DC and Markov policies
    deltaDC  = delta;               % Discount factor for DC case
beta         = 0.995;                % Quasi-hyperbolic discount factor
risk_aver    = 2;                   % CRRA parameter (risk aversion)
                                    % computational notes: When different from 1, restarting is extremely effective (run for a while to reduce loss, cancel, run again with loadoutput=1)
                                    %                      Note that there does not appear to exist smooth solutions except when risk_aver=1

%% Production
A        = 5;                       % Technology index
capshare = 0.36;                    % Capital share
depr     = 1;                       % Rate of depreciation

%% Bounds
a           = 0;                    % Lower bound (borrowing limit)
b           = 6;                    % Upper bound for assets (best practice is to compute larger interval than interval of interest)
bspecfig    = b-2;                  % Separate upper bound for Figure 1 (due to potential imprecisions near upper bound)

%% Labor productivity process
lmin = 1;                           % Low productivity
p    = 0;                           % Probability of low state
lmax = (1-p*lmin)/(1-p);           % High productivity (this implies mean productivity =1)

% Production ("Income production") function
f = @(a) A.*a.^(capshare)+(1-depr)*a; % f(k_t)=k_{t+1}+c_t (full depreciation). Adjust as appropriate

% offset = "\hat{x}(x) that uniquely solves  f(\hat{x},\underline{l})=f(x,\overline{l})" (see paper).  
% in deterministic models, insert a dummy (this can be used in all deterministic models)
offset = @(ass) ass;                % dummy



%% Choice of algorithm
algorithm   = 'tsa';           % Options:
                                   % 'default'   %   Default (optimized CMA-ES for monotone policies)
                                   % 'tsa'       %   Hybrid CMA-ES-policy iteration algorithm (sometimes very effective in stochastic models and is preferred in dynamically consistent models)
                                   
                                   %[ 'ssr'       %   State Space recursion                                    
                                                 %   NOT fully implemented
                                                 %   --- to use import
                                                 %   steady states (see folder /steady
                                                 %   states --- Online
                                                 %   Appendix Section 1)]
                                  
                                   
%% Precision parameters
dissaving   = 0;                    % = 1 restricts second-mover to weak dissaving policies in Lossfunction
gridprec    = 150;                   % Grid intervals per unit of assets
dsp         = 2;                    % Second-mover direct search precision (if increased to 5 or 10, finds Krusell-Smith solutions)
innerstep   = 0.1 / sqrt(gridprec*(b-a)); % Initial step length for CMA-ES
tol         = 1e-8;                 % Tolerance (value-function iteration)
tola        = 1e-14;                % CMA-ES loss tolerance (default variant only)
MaxIt       = 5000;                 % Maximum iterations

%% Internal settings (changes from now on might lead to errors)
n = round(gridprec*(b-a),0);        % Number of grid intervals
ssrdc       = 0;                    % if = 1, Loss algorithm is used to compute DC and naive solutions in ssr case
RESTART     = 1;                    % 0: iteration from opts; 1: find initial candidate via supervised policy iteration

if figonly == 1
    cutoff = 300;
    MaxIt = 0;
    tol = 1e-20;
    shortv = 0;
    loadoutput = 1;
end

%% Display options
howkaplan = 0;                     % Shows standard Kaplan aiyagari simulation (if relevant)
%showSM     = 1;                    % Shows second-mover in main plot
shownaive  = 0;                     % Shows naive policy in main plot (requires recoding for optimal growth model)
showDC     = 0;                     % Shows dynamically consistent solution
bspecfigr   = b;                    % Separate b for other figures
bspecfigcv  = b;
bspecfigcss = b;

%% Utility function and discount sequence
if risk_aver == 1
    u = @(c) log(c);
else    
    u = @(c) (c.^(1-risk_aver)-1)./(1-risk_aver);
end

% Discount sequence
cutoff = 150;                       % Cutoff for finite-sum approximations
dis = @(j) beta.*(delta.^j);
dis = dis(1:cutoff);

if exist('r', 'var')
    disp(['Interest rate: ' num2str(100*r) '%, modified golden rule interest rate ' num2str(100*(1/delta-1)) '%, Wage: ' num2str(W) ]);
    disp(['Borrowing limit: ' num2str(a) ', Consumption at borrowing limit ' num2str(r*a+W*lmin) ', Natural borrowing limit: ' num2str(-((W*lmin)/(r)))]);
end

if (p==0 || p==1)
    disp(['Deterministic model (p = ' num2str(p) ')']);
else
    disp(['Probability of low labor endowment: p = ' num2str(p) ' , low / high endowment ' num2str(lmin) ' / ' num2str(lmax)]);
end

%% Retrieve dynamically consistent solution (Kaplan code)
if exist('r', 'var')
    disp('Computing saving functions in dynamically consistent (DC) case with Kaplan EGP algorithm')
    if calmode == 1 && ~strcmp(algorithm, 'ssr')
        [dyncon_kaplan, grid_kaplan] = aiyagari(a, b, n, r, W, lmin, lmax, risk_aver, p, deltaDC);
    else 
        [dyncon_kaplan, grid_kaplan] = aiyagari(a, b, n, r, W, lmin, lmax, risk_aver, p, delta);
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

if loadoutput == 1
    if exist(optssave, 'file') == 2
        load(optssave); 
        disp(['Initiating from opts saved in ' optssave])
    else
        disp(['No ' optssave ' file found. Please make sure loadoutput = 0 and figonly = 0']);
        return;
    end
else
    opts = a*ones(1, n+1);
end

%% Executing algorithm
tic;

if MaxIt>0
if strcmp(algorithm, 'ssr')
    disp('Initiating State Space Recursion Algorithm')
    cct   = 1;  
    dissaving= 0;
    con   = 0;
    branch = 1;
    disp(['Computing Markov policy with beta = ' num2str(beta)]);
    [SS, opts, Br] = algorithm_ssr(a, b, n, u, dis, cct, opts, dissaving, f, offset, p, tol, con, branch);
    ids = transpose(Br);
    %save(optssave,'opts');
    if ssrdc==1
        disp('Computing Markov policy with beta = 1 (dynamically consistent solution)');
        dis = @(j) (delta.^j);
        dis = dis(1:cutoff);
        [~, dyncon, ~] = algorithm_ssr(a, b, n, u, dis, cct, opts, dissaving, f, offset, p, tol, con, branch);
        dyncon = transpose(dyncon);
        disp('Computing Naive solution');
        dis = @(j) beta.*(delta.^j);
        dis = dis(1:cutoff);
        [~, naive] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 0, tol, f, offset,dsp);
        naive = transpose(naive);
    end
else
    if beta ==1
        disp(['Initiating ' algorithm ' CMA-ES and computing Markov policy with beta = ' num2str(beta) ' (dynamic consistency)']);
    else
        disp(['Initiating ' algorithm ' CMA-ES and computing Markov policy with beta = ' num2str(beta) ' (dynamic inconsistency)']);
    end
    opts = algorithm_handle(a, b, n, MaxIt, innerstep, opts, RESTART, dissaving, dis, p, u, tol, f, offset,dsp,tola);
    disp('Computing naive solution')
end

time = round(toc);
disp(['...Total Ego Loss algorithm run time: ' num2str(time) ' seconds.']);

if beta == 1 && figonly == 0
    dyncon = opts;
    save(dynconsave,'dyncon');
elseif figonly == 0
    save(optssave,'opts');
end
end
 
if shownaive==1
    [~, naive] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 0, tol, f, offset,dsp);
    naive = transpose(naive);
else
    [MAN, asimN, ysimN, TsimN]= deal([]);
end


[~, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 0, tol, f, offset,dsp);
ids = transpose(Br);
   
%% Plotting results
% colors
str = '#000080';
nblue = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
ngreen = [0 0.5 0];
str = '#800000';
nred = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255; % better colors

figure(1);
fplot(@(j) j, [a, b], ':', 'linewidth', 1, 'Color', uint8([5 5 5]), 'HandleVisibility', 'off')
hold on

x1 = a:((b-a)/n):b;
t1 = 1:1:n+1;
%if showSM == 1    
    [~, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 0, tol, f, offset,dsp);
    ids = transpose(Br);    
    p4=plot(x1, ids, '-', 'LineWidth', 3, 'color', nred, 'DisplayName', '$g^I$ (Second-mover)');
    hold on
%end      

hold on   

if showDC==1
    displayNameDC = sprintf('Dynamically Consistent ($\\beta = 1$)');
    p1 = plot(x1, dyncon, '-', 'color', ngreen, 'LineWidth', 3, 'DisplayName', displayNameDC);
end
 
displayNameMP = sprintf('$g^E$ (First-mover)');
p3 = plot(x1, opts(t1), '-', 'LineWidth', 3, 'color', nblue, 'DisplayName', displayNameMP);

if risk_aver == 1 && depr == 1  
    G=linspace(a,b,n+1);
    z=(capshare.*beta.*delta.*A)./(1-capshare.*delta.*(1-beta)).*(G).^(capshare); % analytical solution from K-S 2002/2008
    h3=plot(x1,z,'--','LineWidth',2.5,'Color', ngreen,'DisplayName','Symmetric Minimax Equilibrium');
    legend([h3,p3,  p4], 'location', 'northwest','FontSize',16, 'Interpreter', 'latex');     
else
    legend([p3 p4], 'location', 'northwest','FontSize',16, 'Interpreter', 'latex');
end
legend boxoff

xlabel('Assets at date $t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Assets at date $t+1$', 'Interpreter', 'latex', 'FontSize', 16);
xlim([a bspecfig]); 

hold off