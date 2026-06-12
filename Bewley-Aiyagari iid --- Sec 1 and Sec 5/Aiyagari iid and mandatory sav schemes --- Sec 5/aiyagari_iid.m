clear;
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
addpath(fullfile(pwd, 'figurescripts'));
addpath(fullfile(pwd, 'mainscripts'));



%% Saving Mandate (forced values from dyncon which must be saved by running first with beta=1) below sm:
sm = 0;

%% Initiation, storing and loading
loadoutput  = 1;                    % 0: initiates from horizontal policy; 1: load opts from output.mat with handle specified next
    Handle = 'Aiyagari_iid';        % Appended to save and load filenames (saves
                                    % in "savedmatfiles" with a separate file called dyncon is saved when beta=1)
figonly     = 1;                    % = 1 draws figures from saved solution(s) [no algorithm call]
shortv      = 0;                    % = 1 prints mean and breaks after first figure    
%% Preferences
delta        = 0.95;                % Discount factor
calmode      = 1;                   % Calibration mode to compare DC and Markov policies
    deltaDC  = 0.95;                % Discount factor for DC case
beta         = 0.88;                % Quasi-hyperbolic discount factor
risk_aver    = 1;                   % Rate of risk aversion

%% Production
A        = 1.00;                    % Technology index
capshare = 0.36;                    % Capital share
depr     = 0.08;                    % Rate of depreciation
k        = 5.31;                    % Capital-labor ratio

%% Bounds
a           = 0;                    % Lower bound (borrowing limit)
b           = 25;                   % Upper bound for assets
    bspecfig    = 8;                % Separate upper bound for Figure 1

%% Labor productivity process
lmin = 0.15;                        % Low productivity
p    = 0.5;                         % Probability of low state
lmax = (1-p*lmin)/(1-p);           % High productivity <=> p lmin+(1-p)lmax=1 (normalization)

% Prices (from production function)
F  = @(a) A*a^(capshare);          % Production function
DF = @(a) A*capshare*a^(capshare-1); % MPK

r = DF(k)-depr;                     % Interest rate (R-1)
W = F(k) - DF(k)*k;                 % Wage rate

% Income functions
f = @(ass) (1+r)*ass+W*lmin;        % Income in low state
offset = @(ass) ass+(lmax-lmin)*W/(1+r); % Offset function for high state

%% Output/figure settings
showval     = 0;                    % = 1 displays value function among figures
extrads     = 1;                    % Full suite of inequality statistics (Teil, Atkinson, Palma, etc)
extratailm  = 0;                    % Full tail analysis

%% Choice of algorithm
algorithm   = 'default';            

%% Precision parameters
dissaving   = 1;                    % = 1 restricts second-mover to weak dissaving policies in Lossfunction
gridprec    = 20;                   % Grid intervals per unit of assets
dsp         = 4;                    % Second-mover direct search precision multiplier
innerstep   = 1 / sqrt(gridprec*(b-a)); % Initial step length for CMA-ES
tol         = 1e-8;                 % Tolerance (value-function iteration)
MaxIt       = 5000;                 % Maximum iterations
if beta == 1
    tola = 1e-8;
else
    if sm ~=a
        tola = 2.5e-3;                  % very high precision fine tuned for paper's model to ensure stabilization. With powerful CPU or GPU support, set lambda (much) higher in default. The impact will be very marginal however (especially for distributions since these stabilize before means)
    else                           
        tola = 2.0e-3;                  % without mandate (sm=0), 2.0e-3 or lower is recommended
    end
end

%% Internal settings (changes from now on might lead to errors)
n = round(gridprec*(b-a),0);        % Number of grid intervals

renewaldc   = 0;                    % if = 1, Ego Loss algorithm is used to compute DC and naive solution in renewal case
RESTART     = 1;                    % 0: iteration from opts; 1: find initial candidate via supervised policy iteration

if figonly == 1
    cutoff = 300;
    MaxIt = 0;
    tol = 1e-20;
    shortv = 0;
    loadoutput = 1;
end

%% Display options
showkaplan = 0;                     % Shows standard Kaplan aiyagari simulation (if relevant)
shownaive  = 0;                     % Shows naive policy in main plot (only low productivity)
bspecfigr   = 15;                   % Separate b for other figures
bspecfigcv  = 15;
bspecfigcss = 10;

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
    disp(['Capital-output ratio  ' num2str(k/F(k))]);
    disp(['Interest rate: ' num2str(100*r) '%, modified golden rule interest rate ' num2str(100*(1/delta-1)) '%, Wage: ' num2str(W) ]);
    disp(['Borrowing limit: ' num2str(a) ', Consumption at borrowing limit ' num2str(r*a+W*lmin) ', Natural borrowing limit: ' num2str(-((W*lmin)/(r)))]);
end

if (p==0 || p==1)
    disp(['Deterministic model (p = ' num2str(p) ')']);
else
    disp(['Probability of low labor endowment: p = ' num2str(p) ' , low / high endowment ' num2str(lmin) ' / ' num2str(lmax)]);
end

if sm > a
    disp(['Mandatory saving scheme enabled for asset levels below ' num2str(sm)]);
end

%% Retrieve dynamically consistent solution (Kaplan code)
if exist('r', 'var')
    disp('Computing saving functions in dynamically consistent (DC) case with Kaplan EGP algorithm')
    [dyncon_kaplan, grid_kaplan] = aiyagari(a, b, n, r, W, lmin, lmax, risk_aver, p, delta);
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
 
%% Executing algorithm
tic;

if beta == 1 
        dyncon = a*ones(1, n+1);
        naive  = dyncon;
    else
        if strcmp(algorithm, 'renewal')  
            if renewaldc==0
                dyncon = dyncon_kaplan;
                dyncon = interp1(grid_kaplan, dyncon, linspace(a, b, n+1), 'linear');
            end
        else
            if isfile(dynconsave) %&& renewaldc==0
                load(dynconsave);
            elseif renewaldc == 0
                disp('Please run with beta=1 to generate dynamically consistent solution for comparison plots');
                return;
            end
            if renewaldc==0
            dyncon = interp1(linspace(a, b, length(dyncon)), dyncon, linspace(a, b, n+1), 'linear');
            end 
        end
        if renewaldc==0
        [~, Br] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 0, tol, f, offset,dsp,sm);
        naive = transpose(Br);
        end
end


if loadoutput == 1
    if exist(optssave, 'file') == 2
        load(optssave); 
        opts = interp1(linspace(a, b, length(opts)), opts, linspace(a, b, n+1), 'linear');
        disp(['Initiating from opts saved in ' optssave])
    else
        disp(['No ' optssave ' file found. Initiating from Horizontal Policy at borrowing limit']);
        opts = a*ones(1, n+1);
        RESTART = 1;
end    
elseif loadoutput == 0
    opts = a*ones(1, n+1);
    if ~strcmp(algorithm, 'renewal')
    disp('Initiating from Horizontal Policy at borrowing limit')
    end
end
if MaxIt>0     
    if beta ==1
        disp(['Initiating ' algorithm ' CMA-ES and computing Markov policy with beta = ' num2str(beta) ' (dynamic consistency)']);
    else
        disp(['Initiating ' algorithm ' CMA-ES and computing Markov policy with beta = ' num2str(beta) ' (dynamic inconsistency)']);
        %if exist(dynconsave, 'file') == 2
        %load(dynconsave);
        dyncon = dyncon_kaplan;%  mandatory saving scheme (to use EGP). Current uses Loss algorithm, interp1(linspace(a, b, length(dyncon)), dyncon, linspace(a, b, n+1), 'linear');
        %else
         %   disp('Please run with beta=1 to generate dyncon.mat for comparison plot');
        %return;
    end
    opts = algorithm_handle(a, b, n, MaxIt, innerstep, opts, RESTART, dissaving, dis, p, u, tol, f, offset,dsp,dyncon,sm,tola);
    disp('Computing naive solution')
    [~, naive] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 0, tol, f, offset,dsp,sm);
    naive = transpose(naive);
end

time = round(toc);
disp(['...Total Ego Loss algorithm run time: ' num2str(time) ' seconds.']);
 
%% saving

%if isempty(Handle)
 %   optssave = 'output.mat';
 %   dynconsave = 'dyncon.mat';
%else
%    optssave = sprintf('output_%s.mat', Handle);
%    dynconsave = sprintf('dyncon_%s.mat', Handle);
%end



% Save the file
%save(filename, 'myVariable');

if beta == 1 && figonly == 0
    dyncon = opts;
    save(dynconsave,'dyncon');
elseif figonly == 0
    save(optssave,'opts');
end
%end
%dyncon = interp1(grid_kaplan, dyncon_kaplan, linspace(a, b,length(dyncon))', 'linear');
[~, VF]   = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 2, tol, f, offset,dsp,sm); 
[~, VFDC] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 2, tol, f, offset,dsp,sm);
if exist('r', 'var')
    [MA, asim, ysim, Tsim] = stationary_dis(a, b, n, r, W, lmin, lmax, p, opts, f, offset);
    [MAK, asimK, ysimK, TsimK] = stationary_dis(a, b, n, r, W, lmin, lmax, p, dyncon, f, offset);
    [MAKap, asimKap, ysimKap, TsimKap] = stationary_dis(a, b, n, r, W, lmin, lmax, p, interp1(grid_kaplan, dyncon_kaplan, linspace(a, b,length(dyncon))', 'linear'), f, offset);
end

if shownaive==1
    %[~, naive] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 0, tol, f, offset,dsp);
    %naive = transpose(naive);
    [MAN, asimN, ysimN, TsimN] = stationary_dis(a, b, n, r, W, lmin, lmax, p, naive, f, offset);
    disp(['Naive mean is: ' num2str(MAN) ]);
else
    [MAN, asimN, ysimN, TsimN]= deal([]);
end
%dyncon = interp1(grid_kaplan, dyncon_kaplan, linspace(a, b,length(dyncon))', 'linear'); % optional overwrite of dyncon 

%if showSM == 1
    [~, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 0, tol, f, offset,dsp,sm);
    ids = transpose(Br);
%     
%else
%    ids = [];
%end



%% Plotting results
addpath(fullfile(pwd, 'figurescripts'));
% saving workspace in structure for easy export
%% Compute grid and income values

% Translation from offset-based to state-separated policies for figures_random.m
G = linspace(a, b, n+1)';
lsim=ysim;
lsimK=ysimK;
lsimKap=ysimKap;
temp = interp1(G, opts(:), offset(G), 'linear', 'extrap');
opts = [opts(:); temp(:)];
temp = interp1(G, dyncon(:), offset(G), 'linear', 'extrap');
dyncon = [dyncon(:); temp(:)];
if exist('naive', 'var') && ~isempty(naive)
    temp = interp1(G, naive(:), offset(G), 'linear', 'extrap');
    naive = [naive(:); temp(:)];
end
if exist('ids', 'var') && ~isempty(ids)
    temp = interp1(G, ids(:), offset(G), 'linear', 'extrap');
    Br = [ids(:); temp(:)];
end

vars = who;
S = cell2struct(cellfun(@(v) evalin('caller', v), vars, 'UniformOutput', false), vars, 1);

% colors
str = '#000080';
nblue = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
ngreen = [0 0.5 0];
str = '#800000';
nred = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255; % better colors

if p~=1 && p~=0 
    if exist('r','var')%% Plotting random case
    figures_random(S);
    else  
    %% FIGURE 1: Policy functions
    
figure(1); clf;
policy_QH_low = @(x) interp1(linspace(a, b, length(opts)), opts, x, 'linear', 'extrap');
policy_QH_high = @(x) interp1(linspace(a, b, length(opts)), opts, offset(x), 'linear', 'extrap');

policy_A_low = @(x) interp1(linspace(a, b, length(dyncon)), dyncon, x, 'linear', 'extrap');
policy_A_high = @(x) interp1(linspace(a, b, length(dyncon)), dyncon, offset(x), 'linear', 'extrap');

policy_N_low = @(x) interp1(linspace(a, b, length(naive)), naive, x, 'linear', 'extrap');
policy_N_high = @(x) interp1(linspace(a, b, length(naive)), naive, offset(x), 'linear', 'extrap');
displayNameDC = sprintf('Dyn. Cons. ($A_t\\in \\{ \\underline{A},\\overline{A} \\}$)');
displayNameDCH= sprintf('Dynamically Consistent ($l_t=\\overline{l}$)');
displayNameN = sprintf('Na\\"{i}ve Solution ($\\beta = %s$)', num2str(beta));
displayNameMP  = sprintf('Quasi-hyperbolic ($A_t=\\underline{A}$)');
displayNameMPH = sprintf('Quasi-hyperbolic ($A_t=\\overline{A}$)');

x1 = linspace(a,b,n+1);

%glow  = @(x) interp1(x1, opts, x, 'linear', 'extrap');
%ghigh = @(x) glow(offset(x));
optsH = policy_QH_high(x1);
dynconH = policy_A_high(x1);
naiveH = policy_N_high(x1);

fplot(@(j) j, [a, b], ':', 'linewidth', 1, 'Color', uint8([5 5 5]), 'HandleVisibility', 'off')
hold on

p1 = plot(x1, dyncon, '-', 'color', "k", 'LineWidth', 3, 'DisplayName', displayNameDC);

if shownaive == 1
    p2 = plot(x1, naive, '--', 'LineWidth', 3, 'Color', "k", 'DisplayName', displayNameN);
end

p3 = plot(x1, opts, '-', 'LineWidth', 3, 'color', ngreen, 'DisplayName', displayNameMP);

p4 =plot(x1, optsH, '-', 'LineWidth', 3, 'color', nred,'DisplayName', displayNameMPH);% plot(x_shifted, y1, '-', 'LineWidth', 3, 'color', nred,'DisplayName', displayNameMPH);
p5 = plot(x1, dynconH, '-', 'color', "k", 'LineWidth', 3);
set(p5, 'HandleVisibility', 'off');

if shownaive == 0
    legend([p1,p3,p4], 'location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');
else
    legend([p1, p3, p2], 'location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');
end
%if beta==1
%    plot(grid_kaplan, dyncon_kaplan, 'y--', 'LineWidth', 1.5, 'DisplayName', 'Dyn. cons. (Kaplan)')
%end
legend boxoff
xlabel('Assets at date $t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Assets at date $t+1$', 'Interpreter', 'latex', 'FontSize', 16);
xlim([a bspecfig]);
hold off
drawnow
    end
else
    %% plotting deterministic case
    if exist('r', 'var')
        figures_det(a, b, n, opts, dyncon, naive, grid_kaplan, dyncon_kaplan, ...
                    W, r, lmin, lmax, beta, 0, showkaplan, dis, u, p, tol, algorithm, ...
                    VF, VFDC, ids, bspecfig, bspecfigr,bspecfigcv,MA, asim, ysim, Tsim, MAK, asimK, ysimK, TsimK)
        
        %(a, b, n, opts, dyncon, naive, grid_kaplan, dyncon_kaplan, ...
              %      W, r, lmin, lmax, beta, showSM, showkaplan, dis, u, p, tol, algorithm, ...
              %      VF, VFDC, ids, bspecfig, MA, asim, ysim, Tsim, MAK, asimK, ysimK, TsimK)
    else 
        figure(1);
        fplot(@(j) j, [a, b], ':', 'linewidth', 1, 'Color', uint8([5 5 5]), 'HandleVisibility', 'off')
        hold on
        
        
        x1 = a:((b-a)/n):b;
        t1 = 1:1:n+1;
        %if showSM == 1    
        %     
        %    [~, Br] =     Lossfunction(a, b, n, opts, speedy, dis, u, p, 0, tol, f, offset,dsp);
        %    ids = transpose(Br);    
        %    plot(x1, ids, '--', 'LineWidth', 1.5, 'color', nblue, 'DisplayName', 'Second-Mover');
        %    hold on
        %end      
        %y1 = opts(t1);

    if showkaplan==1 && ~strcmp(algorithm, 'renewal')
        plot(grid_kaplan, dyncon_kaplan, 'y--', 'LineWidth', 1.5, 'DisplayName', 'Dynamically consistent (Kaplan)')
    end
    hold on   

    displayNameDC = sprintf('Dynamically Consistent ($\\beta = 1$)');
    p1 = plot(x1, dyncon, '-', 'color', nblue, 'LineWidth', 3, 'DisplayName', displayNameDC);

%displayNameN = sprintf('Na\\"{i}ve Solution ($\\beta = %s$)', num2str(beta));
%p2 = plot(x, naive, '--', 'LineWidth', 3, 'Color', ngreen, 'DisplayName', displayNameN);

    displayNameMP = sprintf('Markov Policy ($\\beta = %s$)', num2str(beta));
    p3 = plot(x1, opts(t1), '-', 'LineWidth', 3, 'color', nred, 'DisplayName', displayNameMP);
    if risk_aver == 1 && isa(F, 'function_handle') && isequal(F(2), A * 2^(capshare))
        G=linspace(a,b,n+1);
        z=(capshare.*beta.*delta.*A)./(1-capshare.*delta.*(1-beta)).*(G).^(capshare); % analytical solution from K-S 2002/2008
        h3=plot(x1,z,'--','LineWidth',2.5,'Color', ngreen,'DisplayName','Symmetric Minimax Equilibrium');
        legend([p1, p3, h3], 'location', 'northwest','FontSize',16, 'Interpreter', 'latex');
    else
        legend([p1, p3], 'location', 'northwest','FontSize',16, 'Interpreter', 'latex');
    end
    legend boxoff

xlabel('Assets at date $t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Assets at date $t+1$', 'Interpreter', 'latex', 'FontSize', 16);
xlim([a bspecfig]); 

hold off
    end
end

clear S;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aiyagari model
% Endogenous Grid Points with IID Income
% Greg Kaplan 2017
% Edited to two-state case by Martin Kaae Jensen 2025
function [dyncons, grid] = aiyagari(a,b,n,r,W,lmin,lmax,risk_aver,p,delta)
  

% asset grids
amax        = b;                       % upper bound
borrow_lim  = a;                       % borrowing limit
na          = n+1;                       % # grid points 
agrid_par   = 1;                           % Linear grid 

% computation
max_iter    = 500;
tol_iter    = 1.0e-6;
Nsim        = 50000;
Tsim        = 400;



% UTILITY FUNCTION
if risk_aver==1
    u = @(c)log(c);
else    
    u = @(c)(c.^(1-risk_aver)-1)./(1-risk_aver);
end    
u1    = @(c) c.^(-risk_aver);
u1inv = @(c) c.^(-1./risk_aver);

% SET UP GRIDS
agrid = linspace(0,1,na)';
agrid = agrid.^(1./agrid_par);
agrid = borrow_lim + (amax-borrow_lim).*agrid;
 
ygrid    = [lmin lmax].';  
ydist    = [p 1-p].';  
ycumdist = cumsum(ydist);
ny       = length(ydist);

% DRAW RANDOM NUMBERS
rng(2017);
yrand = rand(Nsim,Tsim);

% SIMULATE LABOR EFFICIENCY REALIZATIONS
yindsim = zeros(Nsim,Tsim);
for it = 1:Tsim
    yindsim(yrand(:,it)<=ycumdist(1),it) = 1;
    for iy = 2:ny
        yindsim(yrand(:,it)>ycumdist(iy-1) & yrand(:,it)<=ycumdist(iy),it) = iy;
    end
end
ysim = ygrid(yindsim);

iterKL  = 0;
KLdiff  = 1;
R       = 1+r;
yscale  = 1;

conguess = zeros(na,ny);
for iy = 1:ny
    conguess(:,iy) = r.*agrid + W.*yscale.*ygrid(iy);
end
con = conguess;

iter = 0;
cdiff= 1000;
while iter <= max_iter && cdiff>tol_iter
    iter = iter + 1;
    sav  = zeros(na,ny);
    conlast = con;

    emuc    = u1(conlast)*ydist; 
    muc1    = delta.*R.*emuc; 
    con1    = u1inv(muc1);
    ass1    = zeros(na,ny);

    for iy = 1:ny
        ass1(:,iy) = (con1 + agrid - W.*yscale.*ygrid(iy))./R;
        for ia  = 1:na 
            if agrid(ia)<ass1(1,iy)
                sav(ia,iy) = borrow_lim;
            else
                sav(ia,iy) = lininterp1(ass1(:,iy),agrid,agrid(ia));
            end                
        end
        con(:,iy) = R.*agrid + W.*yscale.*ygrid(iy) - sav(:,iy);
    end

    cdiff = max(max(abs(con-conlast)));
end

dyncons = sav(:,1);
grid    = agrid;

end


function [yi] = lininterp1(x,y,xi)
%% Linear interpolation function from Kaplan (2017). Maintained for comparability
%[yi] = lininterp1(x,y,xi)
%x is N x 1
%y is N x 1
%xi is 1 x 1
%yi is 1 x 1
%same as lininterp except is faster and only looks up one function
%extrapolates out of range

placeLow = find(xi<x,1)-1;
if placeLow == 0
    placeLow = 1;
end
if isempty(placeLow)
    placeLow = length(x)-1;
end    
placeHigh = placeLow+1;
xLow    = x(placeLow);
xHigh   = x(placeHigh);
yLow    = y(placeLow);
yHigh   = y(placeHigh);

yi = yLow +(xi-xLow).*(yHigh-yLow)./(xHigh-xLow);
end
function [Mean, asim, ysim, Tsim] = stationary_dis(a, b, n, R, W, lmin, lmax, p, s, f, offset)
% Monte-Carlo simulation of stationary distribution, based on Kaplan (2017)
% Modified 2025 with burn-in period

%% setting up grid and policy
% Determine the adjusted upper bound by testing the offset function
% We need to find b_adj such that offset(b_adj) <= b
% This prevents the policy function from being queried outside its defined range
test_pts = linspace(a, b, 100);
max_safe_val = b;
for i = 1:length(test_pts)
if offset(test_pts(i)) <= b
 max_safe_val = test_pts(i);
else
break; % Stop when first point out of range determined
end
end
b_adj = max_safe_val;
% Recreate grid
agrid_full = a + (b - a) .* linspace(0, 1, n+1)'; % Original full grid
agrid = agrid_full(agrid_full <= b_adj); % Restrict grid to valid points
if isempty(agrid) || agrid(end) < b_adj
 agrid = [agrid; b_adj]; % Add b_adj
end
% Define interpolation functions using the provided policy function s
fmin = @(x) lininterp1(agrid, s, x); % Low state (no shock)
fmax = @(x) fmin(offset(x)); % High state (with shock)
% Precompute savings decisions for both states
for j = 1:length(agrid)
 sav(j, 1) = fmin(agrid(j)); % Low labor state
 sav(j, 2) = fmax(agrid(j)); % High labor state
end
%% The remainder follows Kaplan (2017) but with discrete distribution
% labor realization simulation
burn_in = 1000;
Nsim = 50000;
Tsim = 500;
Tsim_full = Tsim + burn_in;
% Define labor supply grid and distribution
ygrid = [lmin lmax]';
ydist = [p 1-p]';
ycumdist = cumsum(ydist);
ny = length(ydist);
% DRAW RANDOM NUMBERS
rng(2017); % ensures reproducability by always making same draws
yrand = rand(Nsim, Tsim_full);
yindsim = 1 + (yrand > p);
ysim_full = ygrid(yindsim);
asim_full = zeros(Nsim, Tsim_full);
% Create interpolating function
for iy = 1:ny
 savinterp{iy} = griddedInterpolant(agrid, sav(:, iy), 'linear', 'nearest');
end
%% Simulating assets
% Loop over time periods
for it = 1:Tsim_full
if it < Tsim_full
for iy = 1:ny
 asim_full(yindsim(:, it) == iy, it+1) = savinterp{iy}(asim_full(yindsim(:, it) == iy, it));
end
end
end
%% Extract post-burn-in data
asim = asim_full(:, burn_in+1:end);
ysim = ysim_full(:, burn_in+1:end);
%% Mean assets, labor supply
Mean = mean(asim(:, Tsim))/mean(ysim(:, Tsim));
end