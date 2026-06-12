clear;
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
addpath(fullfile(pwd, 'figurescripts'));
addpath(fullfile(pwd, 'mainscripts'));

%% Saving mandate

sm = 1.75;% sm=1.75 in comparison; %set the limit below which households are mandated to save at dynamically consistent level

%% Preferences 
delta     = 0.95;  
beta       = 0.88;   
risk_aver = 1;         

if beta == 1 % for consistency (no mandate imposed on DC solution)
    sm = 0;
end;
%% Initiation

figonly     = 0    ;  % =1 draws figures from precumputed solutions
shortv      = 0;      % breaks after first figure
loadoutput  =  0;     % 0: horizontal policy; 1: load opts from output.mat
Handle = 'intro_fig_panel_B'; % Appended to save and load (can be left empty/%'ed out leading to repeated overwrite)

RESTART     = 0;       % 0: iteration from opts; 1: multiple policy iteration steps; 2: precisely one policy iteration step
calmode      =0;       % calibration mode to compare DC and Markov policies
    deltaDC     = delta; % Discount factor for comparison
algorithm   = 'tsa';
                       % Options are:
                      % 'tsa'       %   tsa (optimized CMA-ES implementation in best response space)            
                      % 'default'   %   default (optimized CMA-ES for monotone policies). 
dissaving = 1;        % if = 1, second-mover is limited to weak dissaving policies in Lossfunction (improves speed and accuracy)
renewaldc   = 0;      % if = 1, Ego Loss algorithm is used to compute DC and naive solution in renewal case. If = 0, DC computed via Kaplan's EGP algorithm
 
if strcmp(algorithm, 'tsa')  
    RESTART = 0;
else 
    RESTART =1;
end
%% Bounds and parameters
a           = 0;  
b           = 10;
bspecfig    = 6;   
bspecfigr   = 8;
bspecfigcv  = 8; 
bspecfigcss = 8;
gridprec    = 20; 
dsp         = 6;


%% Display options

    if figonly == 1 % overwrite for consistency
        shortv = 0;
        loadoutput = 1;
    end
showhisto  = 1;    % =1 shows asset distribution (banded), = 2 uses kernel density estimate
showkaplan = 0;    % shows standard Kaplan aiyagari simulation (if relevant)
showSM     = 1;    % shows second-mover in main plot (only low productivity)
shownaive  = 0;    % shows naive policy in main plot (only low productivity)




% Various deep parameters
n = round(gridprec*(b-a),0);    
innerstep   = 1 / sqrt(gridprec*(b-a)); 
tola=1e-6; % algoritm tolerance (1e-3 already produces very good approximate Markov policies so can be increased)

% Asset grid
t = linspace(0, 1, n+1);
G = a + (b-a)*(t);        % Asset grid
G = G(:);

% Income at each asset level on the grid
%IncomeG_low = f_low(G);                % Income in low state
%IncomeG_high = f_high(G);              % Income in high state
 
if figonly == 1
    cutoff=300;
    MaxIt = 0;
    tol     = 1e-20;
else
    MaxIt   = 5000;
    tol     = 1e-8; 
end
 

%% Utility function and discount sequence
if risk_aver == 1
    u = @(c) log(c);
else    
    u = @(c) (c.^(1-risk_aver)-1)./(1-risk_aver);
end

% Discount sequence
cutoff = 150;                   % cutoff for finite-sum approximations 
dis = @(j) beta.*(delta.^j);
dis = dis(1:cutoff);

%% Labor endowment process
lmin = 0.97;   
p    = 0.5; 
lmax = (1-p*lmin)/(1-p);        % <=> p lmin+(1-p)lmax=1 (normalization)

%% Production and prices
A        = 1.068;               % technology index
capshare = 0.36;                % capital share
depr     = 0.10;                % rate of depreciation    
k        = 5.3;                 % capital-labor ratio

F  = @(a) A*a^(capshare); 
DF = @(a) A*capshare*a^(capshare-1);

disp(['Capital-output ratio  ' num2str(k/F(k))]);
r = DF(k)-depr;
W = F(k) - DF(k)*k;
 

%% optional price overwrite
%r = 0.03; 
%W = 5.674; 


f = @(ass) (1+r)*ass+W*lmin;
offset = @(ass) ass+(lmax-lmin)*W/(1+r); 


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
   %disp(['Note that best practice is to always make a run with beta = 1 to confirm that EGP algorithm produces same results as Ego Loss algorithm in DC case.'])
    if calmode == 1 && ~strcmp(algorithm, 'renewal')
        [dyncon_kaplan, grid_kaplan] = aiyagari(a, b, n, r, W, lmin, lmax, risk_aver, p, deltaDC);
    else 
        [dyncon_kaplan, grid_kaplan] = aiyagari(a, b, n, r, W, lmin, lmax, risk_aver, p, delta);
    end
end

%% save and load handles and folder specifications

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
        [~, Br] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 0, tol, f, offset,dsp);
        naive = transpose(Br);
        end
end


if loadoutput == 1
    if exist(optssave, 'file') == 2
        load(optssave); 
        opts = interp1(linspace(a, b, length(opts)), opts, linspace(a, b, n+1), 'linear');
        disp(['Initiating from opts saved in ' optssave])
    else
        disp(['No ' optssave ' file found. Please make sure loadoutput = 0 and figonly = 0']);
        return;
    end
elseif loadoutput == 2
    if exist(dynconsave, 'file') == 2
        load(dynconsave);
        dyncon = interp1(linspace(a, b, length(dyncon)), dyncon, linspace(a, b, n+1), 'linear');
        %agrid  = linspace(a, b, n+1);
        %disp('Initiating from Full Commitment solution')
    else
        disp('Please run with beta=1 to generate dyncon.mat');
        return;
    end
    opts = dyncon;
    disp('Initiating from Full Commitment solution')
elseif loadoutput == 3
    load(dynconsave);
    dyncon = interp1(agrid, dyncon, linspace(a, b, n+1), 'linear');
    agrid  = linspace(a, b, n+1);
    [~, Br] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 0, tol, f, offset,dsp);
    naive = transpose(Br);
    opts  = naive;
    disp('Initiating from Naive solution')
elseif loadoutput == 0
    opts = a*ones(1, n+1);
    if ~strcmp(algorithm, 'renewal')
    disp('Initiating from Horizontal Policy at borrowing limit')
    end
elseif loadoutput == 4
    load(optssave)
    disp('Initiating from opts saved in output.mat')
    opts = interp1(linspace(a, b, length(opts)), opts, linspace(a, b, n+1), 'linear');
    for j = 1:n
        opts(j+1) = max(opts(j), opts(j+1));
    end
end
if MaxIt>0
if strcmp(algorithm, 'renewal')
    disp('Initiating Renewal-Based Numerical Algorithm')
    cct   = 1; 
    con   = 1; % = 0 uses constant policy for continuations =1 uses opts
    branch = 0; % = 0 uses value-function iteration, = 1 uses summation
    disp(['Computing Markov policy with beta = ' num2str(beta)]);
    [SS, opts, Br] = algorithm_renewal(a, b, n, u, dis, cct, opts, dissaving, f, offset, p, tol, con, branch);
    ids = transpose(Br);
    %save(optssave,'opts');
    if renewaldc==1
        disp('Computing Markov policy with beta = 1 (dynamically consistent solution)');
        disdc = @(j) (delta.^j);
        disdc = disdc(1:cutoff);
        [~, dyncon, ~] = algorithm_renewal(a, b, n, u, disdc, cct, opts, dissaving, f, offset, p, tol, con, branch);
        dyncon = transpose(dyncon);
        disp('Computing Naive solution');
        [~, naive] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 0, tol, f, offset,dsp);%Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 0, tol, f, offset,dsp);
        naive = transpose(naive);                
    end
else
    if beta ==1
        disp(['Initiating ' algorithm ' CMA-ES and computing Markov policy with beta = ' num2str(beta) ' (dynamic consistency)']);
    else
        disp(['Initiating ' algorithm ' CMA-ES and computing Markov policy with beta = ' num2str(beta) ' (dynamic inconsistency)']);
        %if exist(dynconsave, 'file') == 2
        %load(dynconsave);
        dyncon = dyncon_kaplan;%  mandatory saving schemes change from interp1(linspace(a, b, length(dyncon)), dyncon, linspace(a, b, n+1), 'linear');
        %else
         %   disp('Please run with beta=1 to generate dyncon.mat for comparison plot');
        %return;
    end
    opts = algorithm_handle(a, b, n, MaxIt, innerstep, opts, RESTART, dissaving, dis, p, u, tol, f, offset,dsp,dyncon,sm,tola);
    disp('Computing naive solution')
    [~, naive] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 0, tol, f, offset,dsp);
    naive = transpose(naive);
end
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
[~, VF]   = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 2, tol, f, offset,dsp); 
[~, VFDC] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 2, tol, f, offset,dsp);
if exist('r', 'var')
    [MA, asim, ysim, Tsim] = MeanAss_kaplan(a, b, n, r, W, lmin, lmax, p, opts, f, offset);
    [MAK, asimK, ysimK, TsimK] = MeanAss_kaplan(a, b, n, r, W, lmin, lmax, p, dyncon, f, offset);
    [MAKap, asimKap, ysimKap, TsimKap] = MeanAss_kaplan(a, b, n, r, W, lmin, lmax, p, interp1(grid_kaplan, dyncon_kaplan, linspace(a, b,length(dyncon))', 'linear'), f, offset);
end



if shownaive==1
    %[~, naive] = Lossfunction(a, b, n, dyncon, dissaving, dis, u, p, 0, tol, f, offset,dsp);
    %naive = transpose(naive);
    [MAN, asimN, ysimN, TsimN] = MeanAss_kaplan(a, b, n, r, W, lmin, lmax, p, naive, f, offset);
    disp(['Naive mean is: ' num2str(MAN) ]);
else
    [MAN, asimN, ysimN, TsimN]= deal([]);
end
%dyncon = interp1(grid_kaplan, dyncon_kaplan, linspace(a, b,length(dyncon))', 'linear'); % optional overwrite of dyncon 

if showSM == 1
    [~, Br] = Lossfunction(a, b, n, opts, dissaving, dis, u, p, 0, tol, f, offset,dsp);
    ids = transpose(Br);
    idic=0;
    idic = idic + sum(ids(2:n+1) < ids(1:n));
    %disp(['SM strategy decreased on ' num2str(100*idic/n) ' per cent of grid intervals (if this is high, consider increasing number of grid points).'])
else
    ids = [];
end



%% Plotting results
addpath(fullfile(pwd, 'figurescripts'));
% saving workspace in structure for easy export
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
                    W, r, lmin, lmax, beta, showSM, showkaplan, dis, u, p, tol, algorithm, ...
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
        if showSM == 1    
             
            [~, Br] =     Lossfunction(a, b, n, opts, speedy, dis, u, p, 0, tol, f, offset,dsp);
            ids = transpose(Br);    
            plot(x1, ids, '--', 'LineWidth', 1.5, 'color', nblue, 'DisplayName', 'Second-Mover');
            hold on
        end      
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

function [Mean, asim, ysim, Tsim] = MeanAss_kaplan(a, b, n, R, W, lmin, lmax, p, s, f, offset)
% NOTES: Monto-Carlo simulation of stationary distribution following Kaplan (2017)
% Modified 2025 to ensure that policies are never queried where not defined (makes
% no difference if b is sufficiently large).

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
        break;  % Stop when first point out of range determined
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
fmin = @(x) lininterp1(agrid, s, x);        % Low state (no shock)
fmax = @(x) fmin(offset(x));                % High state (with shock)

% Precompute savings decisions for both states
for j = 1:length(agrid)
    sav(j, 1) = fmin(agrid(j)); % Low labor state
    sav(j, 2) = fmax(agrid(j)); % High labor state
end

%% The remainder follows Kaplan (2017) but with discrete distribution
% labor realization simulation
Nsim = 50000;
Tsim = 500;

% Define labor supply grid and distribution
ygrid = [lmin lmax]';
ydist = [p 1-p]';
ycumdist = cumsum(ydist);
ny = length(ydist);

% DRAW RANDOM NUMBERS
rng(2017); % ensures reproducability by always making same draws
yrand = rand(Nsim, Tsim);
yindsim = zeros(Nsim, Tsim);
for it = 1:Tsim
    yindsim(yrand(:, it) <= ycumdist(1), it) = 1;
    for iy = 2:ny
        yindsim(yrand(:, it) > ycumdist(iy-1) & yrand(:, it) <= ycumdist(iy), it) = iy;
    end
end
ysim = ygrid(yindsim);
asim = zeros(Nsim, Tsim);

% Create interpolating function
for iy = 1:ny
    savinterp{iy} = griddedInterpolant(agrid, sav(:, iy), 'linear', 'nearest');
end

%% Simulating assets
% Loop over time periods
for it = 1:Tsim
    if it < Tsim
        for iy = 1:ny
            asim(yindsim(:, it) == iy, it+1) = savinterp{iy}(asim(yindsim(:, it) == iy, it));
        end
    end
end

%% Mean assets/labor supply
Mean = mean(asim(:, Tsim))/mean(ysim(:, Tsim));
end