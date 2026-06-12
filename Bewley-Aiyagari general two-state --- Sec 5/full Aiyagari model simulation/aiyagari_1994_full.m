% Aiyagari (1994) Model Implementation
% Uses Endogenous Grid Points method with Rouwenhorst discretization
% Based on Kaplan (2017) with modifications:
% - Uses theoretical means by default (UseTheoretical = 1)
% - Added controls and statistics for stationary distribution verification
% - Implementation by Martin K Jensen (2025)


clear;
close all;

%% PARAMETERS
% preferences
risk_aver   = 2;            % Rate of risk aversion (Aiyagari uses 1,3,5) 
beta        = 0.955;

% production
A           = 1.00;
deprec      = 0.08;
capshare    = 0.36;

% equilibrium method
UseTheoretical = 1;  % SET TO 1 TO USE THEORETICAL MEAN, 0 FOR SIMULATION MEAN (only use for comparison, the problem with simulated mean is that K/L loop targets small changes and so using this to compute stationary eq. is fragile)

%% LABOR ENDOWMENT PROCESS (six-state Rouwenhorst)
rho         = 0.6;
sigma_eps   = 0.2;
ny          = 2;
[P,logy]    = rouwenhorst(ny,rho,sigma_eps); % Discretizes the AR(1) process: log(z_{t+1}) = ρ log(z_t) + σ_eps (1-ρ^2)^(0.5) ε_t where ε_t are i.i.d. N(0,1) 

% levels normalized to mean=1
ygrid       = exp(logy');
ygrid       = ygrid./mean(ygrid)
P

% stationary distribution
[V,D]       = eig(P');
[~,idx]     = min(abs(diag(D)-1));
stat        = V(:,idx); 
stat        = stat./sum(stat);
ydist       = stat';
ycumdist    = cumsum(ydist);



% asset grids
amax        = 20; 
gridprec    = 20;
na          = amax*gridprec;
borrow_lim  = 0;
agrid_par   = 0.5; %1 for linear, 0 for L-shaped.  

% computation
max_iter    = 1000;
tol_iter    = 1.0e-6;
Nsim        = 50000;  % Number of households in simulation (can be reduced for faster execution)
Tsim        = 200;    % Number of time periods in simulation (can be reduced for faster execution)

maxiter_KL  = 100;
tol_KL      = 1.0e-4;
step_KL     = 0.005;
rguess      = 1./beta-1 - 0.001; % a bit lower than inverse of discount rate
 
KLratioguess = 1.20*((1/beta - 1 + deprec)/(A*capshare))^(1/(capshare-1)); % 20% above mgr as initial guess

%% OPTIONS
Display     = 1;
MakePlots   = 1;

% which function to interpolation 
InterpCon = 0;
InterpEMUC = 1;

% tolerance for non-linear solver
options = optimset('Display','Off','TolX',1.0e-6);

%% UTILITY FUNCTION
if risk_aver==1
    u = @(c)log(c);
else    
    u = @(c)(c.^(1-risk_aver)-1)./(1-risk_aver);
end    

u1 = @(c) c.^(-risk_aver);
u1inv = @(u) u.^(-1./risk_aver);

%% SET UP GRIDS
% assets
agrid = linspace(0,1,na)';
agrid = agrid.^(1./agrid_par);
agrid = borrow_lim + (amax-borrow_lim).*agrid;

%% DRAW RANDOM NUMBERS
rng(2017);
yrand = rand(Nsim,Tsim);

%% SIMULATE LABOR EFFICIENCY REALIZATIONS
if Display >=1
    disp(['Simulating labor efficiency realizations']);
end
yindsim = zeros(Nsim,Tsim);

% Pre-compute cumulative transition matrix for efficiency
cumP = cumsum(P,2);
    
% Initialize first period according to stationary distribution
for i = 1:Nsim
    yindsim(i,1) = find(yrand(i,1) <= ycumdist, 1);
end

% Generate subsequent periods using vectorized 6-state Markov transitions
for it = 2:Tsim
    prev = yindsim(:,it-1);
    U = yrand(:,it);
    for k = 1:ny
        idx = prev==k;
        yindsim(idx,it) = 1 + sum(U(idx) > cumP(k,:),2);
    end
end
    
ysim = ygrid(yindsim);

%% ITERATE OVER KL RATIO
KLratio = KLratioguess;

iterKL = 0;
KLdiff = 1;

while iterKL <= maxiter_KL && abs(KLdiff)>tol_KL
    iterKL = iterKL + 1;

    r   = A.*capshare.*KLratio^(capshare-1) - deprec;
    R   = 1+r;
    wage= A.*(1-capshare).* KLratio^capshare;

    % rescale efficiency units of labor so that output = 1
    yscale = 1;
    
    % initialize consumption function in first iteration only
    if iterKL==1
        conguess = zeros(na,ny);
        for iy = 1:ny
            conguess(:,iy) = r.*agrid + wage.* yscale.*ygrid(iy);
        end
        con = conguess;
     end

    % solve for policy functions with EGP
    iter = 0;
    cdiff = 1000;
    while iter <= max_iter && cdiff>tol_iter
        iter = iter + 1;
        sav = zeros(na,ny);

        conlast = con;

        % Calculate expected marginal utilities for each current state (vectorized)
        emuc = u1(conlast) * P';
        muc1 = beta.*R.*emuc;
        
        % Calculate state-specific consumption
        con1 = zeros(na, ny);
        for iy = 1:ny
            con1(:,iy) = u1inv(muc1(:,iy));
        end

        % loop over income
        for iy = 1:ny
            ass1(:,iy) = (con1(:,iy) + agrid - wage.* yscale.*ygrid(iy))./R;

            % loop over current period assets
            for ia  = 1:na 
                if agrid(ia)<ass1(1,iy) %borrowing constraint binds
                    sav(ia,iy) = borrow_lim;
                else %borrowing constraint does not bind;
                    sav(ia,iy) = lininterp1(ass1(:,iy),agrid,agrid(ia));
                end                
            end
            con(:,iy) = R.*agrid +wage.* yscale.*ygrid(iy) - sav(:,iy);
        end   

        cdiff = max(max(abs(con-conlast)));
        if Display >=2
            disp([' Iteration no. ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end
    end    

    %simulate: start at assets from last iteration
    if iterKL==1
        asim = zeros(Nsim,Tsim);        
    elseif iterKL>1
        asim(:,1) = asim(:,Tsim);
    end
    
    % create interpolating function
    for iy = 1:ny
        savinterp{iy} = griddedInterpolant(agrid,sav(:,iy),'linear');
    end
    
    % loop over time periods
    for it = 1:Tsim
        if Display >=2 && mod(it,100) ==0
            disp([' Simulating, time period ' int2str(it)]);
        end
                
        % asset choice
        if it<Tsim
            for iy = 1:ny
                asim(yindsim(:,it)==iy,it+1) = savinterp{iy}(asim(yindsim(:,it)==iy,it));
            end
        end
    end
    
    % assign actual labor income values
    labincsim = wage.*yscale.*ysim;

    if UseTheoretical
        % Calculate THEORETICAL mean from policy functions
        Q_direct = zeros(na*ny, na*ny);
        
        % Build transition matrix for joint process (assets,state)
        for i = 1:na  % current asset grid point
            for iy = 1:ny  % current income state
                a_next = sav(i,iy);  % policy function gives next period assets
                
                % Find weights for interpolation
                j = find(agrid <= a_next, 1, 'last');
                if isempty(j) || j == 0
                    j = 1;
                end
                
                if j == na
                    wgt = 0;  % No interpolation needed
                else
                    wgt = (a_next - agrid(j)) / (agrid(j+1) - agrid(j));
                end
                
                % Calculate transition to all possible states
                for iy_next = 1:ny
                    prob = P(iy, iy_next);  % Probability of state transition
                    
                    % Current index in the combined state space
                    idx_current = (i-1)*ny + iy;
                    
                    if j == na
                        idx_next = (na-1)*ny + iy_next;
                        Q_direct(idx_current, idx_next) = Q_direct(idx_current, idx_next) + prob;
                    else
                        idx_next_low = (j-1)*ny + iy_next;
                        idx_next_high = j*ny + iy_next;
                        
                        Q_direct(idx_current, idx_next_low) = Q_direct(idx_current, idx_next_low) + prob * (1-wgt);
                        Q_direct(idx_current, idx_next_high) = Q_direct(idx_current, idx_next_high) + prob * wgt;
                    end
                end
            end
        end

        % Solve for stationary distribution of joint process using eigs
        [V,~] = eigs(Q_direct',1,'lm');
        pi_joint = (V / sum(V))';
        if any(imag(pi_joint) ~= 0)
            pi_joint = real(pi_joint);  % Handle numerical precision issues
        end
        
        % Marginalize to get asset distribution
        pi_asset = zeros(1, na);
        for i = 1:na
            for iy = 1:ny
                pi_asset(i) = pi_asset(i) + pi_joint((i-1)*ny + iy);
            end
        end

        % Calculate theoretical mean
        Ea_theory = pi_asset * agrid;
        
        % We still need L from the simulation for consistency
        L = yscale.*mean(ysim(:,Tsim));
        
        % Use theoretical mean instead of simulation mean
        KLrationew = Ea_theory; 
        
        % Also calculate simulation-based KL ratio for comparison
        Ea_sim = mean(asim(:,Tsim));
        KLratio_sim = Ea_sim ./ L;
        
        KLdiff = KLrationew./KLratio - 1;
        if Display >=1
            disp(['Equm iter ' int2str(iterKL), ', r = ',num2str(r), ...
                  ', KL ratio (theory): ',num2str(KLrationew), ...
                  ', KL ratio (sim): ',num2str(KLratio_sim), ...
                  ', KL diff: ',num2str(KLdiff*100) '%']);
        end
    else
        % mean assets and efficiency units using simulation
        Ea_sim = mean(asim(:,Tsim));
        L = yscale.*mean(ysim(:,Tsim));
        
        KLrationew = Ea_sim ./ L;
        
        KLdiff = KLrationew./KLratio - 1;
        if Display >=1
            disp(['Equm iter ' int2str(iterKL), ', r = ',num2str(r), ...
                  ', KL ratio: ',num2str(KLrationew), ...
                  ', KL diff: ',num2str(KLdiff*100) '%']);
        end
    end

    KLratio = (1-step_KL)*KLratio + step_KL*KLrationew; 
end

%% Final results comparison
% Calculate final theoretical mean with higher precision
disp('Computing final theoretical mean with higher precision...');

% Create finer grid for more accurate calculation
na_fine = na;   
agrid_fine = agrid;%linspace(0,1,na_fine)';
%agrid_fine = agrid_fine.^(1./agrid_par);
%agrid_fine = borrow_lim + (amax-borrow_lim).*agrid_fine;

% Interpolate policy functions to finer grid
sav_fine = zeros(na_fine, ny);
for iy = 1:ny
    sav_fine(:,iy) = savinterp{iy}(agrid_fine);
end

% Build transition matrix using direct policy functions
Q_fine = zeros(na_fine*ny, na_fine*ny);

% Build transition matrix for joint process (assets,state) on fine grid
for i = 1:na_fine  % current asset grid point
    for iy = 1:ny  % current income state
        a_next = sav_fine(i,iy);  % policy function gives next period assets
        
        % Find weights for interpolation
        j = find(agrid_fine <= a_next, 1, 'last');
        if isempty(j) || j == 0
            j = 1;
        end
        
        if j == na_fine
            wgt = 0;  % No interpolation needed
        else
            wgt = (a_next - agrid_fine(j)) / (agrid_fine(j+1) - agrid_fine(j));
        end
        
        % Calculate transition to all possible states
        for iy_next = 1:ny
            prob = P(iy, iy_next);  % Probability of state transition
            
            % Current index in the combined state space
            idx_current = (i-1)*ny + iy;
            
            if j == na_fine
                idx_next = (na_fine-1)*ny + iy_next;
                Q_fine(idx_current, idx_next) = Q_fine(idx_current, idx_next) + prob;
            else
                idx_next_low = (j-1)*ny + iy_next;
                idx_next_high = j*ny + iy_next;
                
                Q_fine(idx_current, idx_next_low) = Q_fine(idx_current, idx_next_low) + prob * (1-wgt);
                Q_fine(idx_current, idx_next_high) = Q_fine(idx_current, idx_next_high) + prob * wgt;
            end
        end
    end
end

% Solve for stationary distribution of joint process using eigs
[V,~] = eigs(Q_fine',1,'lm');
pi_joint_fine = (V / sum(V))';
if any(imag(pi_joint_fine) ~= 0)
    pi_joint_fine = real(pi_joint_fine);  % Handle numerical precision issues
end

% Marginalize to get asset distribution
pi_fine = zeros(1, na_fine);
for i = 1:na_fine
    for iy = 1:ny
        pi_fine(i) = pi_fine(i) + pi_joint_fine((i-1)*ny + iy);
    end
end

% Calculate final theoretical mean
final_theory_mean = pi_fine * agrid_fine;
final_sim_mean = mean(asim(:,Tsim));

% Display comparison
disp('FINAL EQUILIBRIUM RESULTS:');
disp(['Interest rate: ' num2str(r*100) '%']);
disp(['Wage: ' num2str(wage)]);
if UseTheoretical
    disp(['Simulation mean: ' num2str(final_sim_mean)]);
    disp(['Theoretical mean (used for equilibrium): ' num2str(final_theory_mean)]);
    disp(['Simulation vs theoretical difference: ' num2str((final_sim_mean-final_theory_mean)/final_theory_mean*100) '%']);
else
    disp(['Simulation mean (used for equilibrium): ' num2str(final_sim_mean)]);
    disp(['Theoretical mean (for comparison): ' num2str(final_theory_mean)]);
    disp(['Theoretical vs simulation difference: ' num2str((final_theory_mean-final_sim_mean)/final_sim_mean*100) '%']);
end

%% Run extended simulation to check convergence to theoretical distribution
disp('Running extended simulation to check convergence to theoretical distribution...');
Tsim_ext = 5000;
asim_ext = zeros(Nsim, 1);
asim_ext(:) = asim(:,Tsim);  % Start from final distribution of main sim

% Extended simulation also needs income state
yindsim_ext = zeros(Nsim, 1);
yindsim_ext(:) = yindsim(:,Tsim);  % Start from final income state

% Track mean during extended simulation
mean_track = zeros(Tsim_ext, 1);
mean_track(1) = mean(asim_ext);

% Extended simulation
for it = 1:Tsim_ext-1
    if mod(it, 1000) == 0
        disp(['Extended simulation period: ' num2str(it)]);
    end
    
    % Generate next period's income using Markov process
    yrand_ext = rand(Nsim, 1);
    
    % Apply vectorized Markov transitions
    prev = yindsim_ext;
    U = yrand_ext;
    for k = 1:ny
        idx = prev==k;
        yindsim_ext(idx) = 1 + sum(U(idx) > cumP(k,:),2);
    end
    
    % Apply policy function based on current state and assets
    for iy = 1:ny
        idx = (yindsim_ext == iy);
        asim_ext(idx) = savinterp{iy}(asim_ext(idx));
    end
    
    % Track mean
    mean_track(it+1) = mean(asim_ext);
end

% Plot convergence path
figure;
plot(1:Tsim_ext, mean_track, 'LineWidth', 1.5);
hold on;
if UseTheoretical
    yline(final_theory_mean, 'r--', 'Theoretical Mean');
    reference_mean = final_theory_mean;
    reference_label = 'theoretical';
else
    yline(final_sim_mean, 'r--', 'Simulation Mean');
    reference_mean = final_sim_mean;
    reference_label = 'simulation';
end
title('Extended Simulation Convergence');
xlabel('Time Periods');
ylabel('Mean Assets');
grid on;
if UseTheoretical
    legend('Extended Simulation', 'Theoretical Mean');
else
    legend('Extended Simulation', 'Simulation Mean');
end

% Calculate final simulation mean
ext_sim_mean = mean_track(end);
disp(['Extended simulation mean after ' num2str(Tsim_ext) ' periods: ' num2str(ext_sim_mean)]);
disp(['Divergence from ' reference_label ' mean: ' num2str((ext_sim_mean-reference_mean)/reference_mean*100) '%']);

%% MAKE PLOTS
if MakePlots == 1 
    figure(2);
    
    % consumption policy function
    subplot(2,4,1);
    plot(agrid,con(:,1),'b-',agrid,con(:,ny),'r-','LineWidth',1);
    grid;
    xlim([0 amax]);
    title('Consumption Policy Function');
    legend('Lowest income state','Highest income state');

    % savings policy function
    subplot(2,4,2);
    plot(agrid,sav(:,1)-agrid,'b-',agrid,sav(:,ny)-agrid,'r-','LineWidth',1);
    hold on;
    plot(agrid,zeros(na,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([0 amax]);
    title('Savings Policy Function (a''-a)');
    
    % income distribution
    subplot(2,4,5);
    bar(ygrid, ydist);
    ylabel('Probability')
    title('Income distribution');
    
    % asset distributions comparison
    subplot(2,4,6:7);
    % Plot both simulation and theoretical distributions
    histogram(asim(:,Tsim), 50, 'Normalization', 'probability', 'FaceColor', [.7 .7 .7], 'EdgeColor', 'black');
    hold on;
    
    % Create theoretical histogram for plotting
    edges = linspace(min(agrid_fine), max(agrid_fine), 51);
    th_hist = zeros(length(edges)-1, 1);
    for i = 1:length(edges)-1
        for j = 1:na_fine
            if agrid_fine(j) >= edges(i) && agrid_fine(j) < edges(i+1)
                th_hist(i) = th_hist(i) + pi_fine(j);
            end
        end
    end
    
    % Plot theoretical histogram
    bar(edges(1:end-1) + diff(edges)/2, th_hist, 'FaceColor', 'none', 'EdgeColor', 'r', 'LineWidth', 1.5);
    
    ylabel('')
    title('Asset Distribution Comparison');
    legend('Simulation', 'Theoretical');
    
    % Mean comparison
    subplot(2,4,8);
    bar([1, 2, 3], [final_sim_mean, ext_sim_mean, final_theory_mean]);
    set(gca, 'XTickLabel', {'Simulation', 'Extended', 'Theoretical'});
    title('Mean Comparison');
    ylabel('Mean Assets');
    grid on;
end

%% CALCULATE GINI COEFFICIENTS AND BORROWING CONSTRAINTS

% Simulation-based Gini coefficients
aysim = asim(:,Tsim);
asset_gini_sim   = calculate_gini(aysim);

% Calculate both gross and net income for simulation
labinc_final = labincsim(:,Tsim);  % Labor income: w*l
gross_incsim = (1+r) * aysim + labinc_final;  % Gross: (1+r)*a + w*l
net_incsim = r * aysim + labinc_final;        % Net: r*a + w*l

% Calculate Gini coefficients for gross and net income
gross_gini_sim  = calculate_gini(gross_incsim);
net_gini_sim    = calculate_gini(net_incsim);

% Fraction of agents at the borrowing constraint (simulation)
frac_constrained = mean(asim(:,Tsim) <= borrow_lim) * 100;

% Theoretical Gini coefficients
asset_gini_theory = calculate_gini_from_distribution(agrid_fine, pi_fine);

% Calculate theoretical income Gini coefficients using joint distribution
% We already have the joint distribution pi_joint_fine and the fine grids
gross_inc_values = [];
net_inc_values = [];
joint_probs = [];

for i = 1:na_fine
    for iy = 1:ny
        idx = (i-1)*ny + iy;
        if pi_joint_fine(idx) > 1e-12  % Only include states with meaningful probability
            asset_val = agrid_fine(i);
            labor_inc = wage * yscale * ygrid(iy);
            
            gross_inc_values = [gross_inc_values; (1+r) * asset_val + labor_inc];
            net_inc_values = [net_inc_values; r * asset_val + labor_inc];
            joint_probs = [joint_probs; pi_joint_fine(idx)];
        end
    end
end

% Normalize probabilities
joint_probs = joint_probs / sum(joint_probs);

% Calculate theoretical Gini coefficients
gross_gini_theory = calculate_gini_from_distribution(gross_inc_values, joint_probs);
net_gini_theory = calculate_gini_from_distribution(net_inc_values, joint_probs);

% Theoretical borrowing constraint fraction
% Count mass at borrowing constraint in theoretical distribution
frac_constrained_theory = 0;
for i = 1:na_fine
    if agrid_fine(i) <= borrow_lim + 1e-6  % Same tolerance as simulation
        for iy = 1:ny
            idx = (i-1)*ny + iy;
            frac_constrained_theory = frac_constrained_theory + pi_joint_fine(idx);
        end
    end
end
frac_constrained_theory = frac_constrained_theory * 100;

% Display results
disp(' ');
disp('INEQUALITY AND BORROWING CONSTRAINT MEASURES:');
fprintf('  Asset Gini          Theoretical: %.4f | Simulation: %.4f\n', asset_gini_theory, asset_gini_sim);
%fprintf('  Gross Income Gini   Theoretical: %.4f | Simulation: %.4f\n', gross_gini_theory, gross_gini_sim);
fprintf('  Income Gini     Theoretical: %.4f | Simulation: %.4f\n', net_gini_theory, net_gini_sim);
fprintf('  Borrowing Constrained Theoretical: %.2f%% | Simulation: %.2f%%\n', frac_constrained_theory, frac_constrained);



%% ROUWENHORST FUNCTION
function [P, logy] = rouwenhorst(n, rho, sigma_eps)
% ROUWENHORST Discretizes an AR(1) process using Rouwenhorst method
%   [P, logy] = rouwenhorst(n, rho, sigma_eps)
%   
%   Inputs:
%     n         : number of grid points
%     rho       : persistence parameter
%     sigma_eps : standard deviation of innovation
%   
%   Outputs:
%     P         : n x n transition matrix
%     logy      : n x 1 grid for log income

% Step size
psi = sqrt(n-1) * sigma_eps / sqrt(1-rho^2);

% Grid
logy = linspace(-psi, psi, n)';

% Base case for n=2
p = (1+rho)/2;
P = [p, 1-p; 1-p, p];

% Recursively build larger transition matrices
for i = 3:n
    P_old = P;
    P = zeros(i,i);
    
    P(1:i-1, 1:i-1) = P(1:i-1, 1:i-1) + p * P_old;
    P(1:i-1, 2:i)   = P(1:i-1, 2:i)   + (1-p) * P_old;
    P(2:i, 1:i-1)   = P(2:i, 1:i-1)   + (1-p) * P_old;
    P(2:i, 2:i)     = P(2:i, 2:i)     + p * P_old;
    
    P(2:i-1, :) = P(2:i-1, :) / 2;
end
end

%% Helper function from original code
function [yi] = lininterp1(x,y,xi)
%% Linear interpolation function from Kaplan (2017)
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

function gini = calculate_gini(values)
    % calculate_gini  Gini coefficient from a raw sample vector
    %   gini = calculate_gini(values) computes the standard Gini index
    %   for the vector values, ignoring NaNs.
    
    % Remove NaNs and sort
    values = values(~isnan(values));
    values = sort(values);
    n = numel(values);
    if n == 0
        gini = NaN;
        return;
    end
    
    % Cumulative approach
    idx = (1:n)';
    total = sum(values);
    gini = (2*sum(idx .* values) - (n+1)*total) / (n * total);
end

function gini = calculate_gini_from_distribution(values, probabilities)
    % calculate_gini_from_distribution  Gini from discrete distribution
    %   gini = calculate_gini_from_distribution(values, probabilities)
    %   where values is a vector of outcomes and probabilities sums to 1.
    
    % Force column vectors and sort by value
    v = values(:);
    p = probabilities(:);
    [v, ix] = sort(v);
    p = p(ix);
    p = p / sum(p);  % Normalize just in case
    
    % Calculate Gini using standard discrete formula
    n = length(v);
    total_income = sum(v .* p);
    
    if total_income == 0
        gini = 0;
        return;
    end
    
    % Standard discrete Gini formula
    gini = 0;
    for i = 1:n
        for j = 1:n
            gini = gini + p(i) * p(j) * abs(v(i) - v(j));
        end
    end
    gini = gini / (2 * total_income);
end