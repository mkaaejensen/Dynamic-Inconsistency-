function figure_det(S)
% FIGURE_DET plots all figures for the deterministic case (p==0 or p==1)
% using the same code as your original deterministic branch,
% but every instance of income computed as (1+S.r)*a+S.W (or variants)
% has been replaced by S.f(a).

% Compute parameters from S.dis:
delta = S.dis(2)/S.dis(1);
beta = S.dis(1)/delta;  % local beta (to avoid conflict with S.beta)
cutoff = 300;


%% Define local colors
str = '#000080';
nblue = sscanf(str(2:end),'%2x%2x%2x',[1 3]) / 255;
ngreen = [0 0.5 0];
str = '#800000';
nred = sscanf(str(2:end),'%2x%2x%2x',[1 3]) / 255;

%% FIGURE 1: Policy functions
figure(1);
fplot(@(j) j, [S.a, S.b], ':', 'linewidth', 1, 'Color', uint8([5 5 5]), 'HandleVisibility', 'off');
hold on;
x = S.a : ((S.b-S.a)/S.n) : S.b;
if S.showSM == 1
    if ~strcmp(S.algorithm, 'renewal')
        [~, Br] = Lossfunction(S.a, S.b, S.n, S.dyncon, S.dissaving, S.dis, S.u, S.p, 0, S.tol, S.f, S.offset, S.dsp);
        ids = transpose(Br);
    end
    plot(x, ids, '--', 'LineWidth', 1.5, 'color', nblue, 'DisplayName', 'Second-Mover');
    hold on;
end

x1 = S.a : ((S.b-S.a)/S.n) : S.b;
t1 = 1:1:S.n+1;
y1 = S.opts(t1);

if S.showkaplan == 1 && ~strcmp(S.algorithm, 'renewal')
    plot(S.grid_kaplan, S.dyncon_kaplan, 'y--', 'LineWidth', 1.5, 'DisplayName', 'Dynamically consistent (Kaplan)');
end
hold on;

displayNameDC = sprintf('Dynamically Consistent ($\\beta = 1$)');
p1 = plot(x, S.dyncon, '-', 'color', nblue, 'LineWidth', 3, 'DisplayName', displayNameDC);

displayNameN = sprintf('Na\\"{i}ve Solution ($\\beta = %s$)', num2str(beta));
p2 = plot(x, S.naive, '--', 'LineWidth', 3, 'Color', ngreen, 'DisplayName', displayNameN);

displayNameMP = sprintf('Markov Policy ($\\beta = %s$)', num2str(beta));
p3 = plot(x1, y1, '-', 'LineWidth', 3, 'color', nred, 'DisplayName', displayNameMP);

legend([p1, p3, p2], 'location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');  
legend boxoff;
xlabel('Assets at date $t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Assets at date $t+1$', 'Interpreter', 'latex', 'FontSize', 16);
xlim([S.a, S.bspecfig]); 
hold off;

%% ====== FIGURE 2: Lifetime Utility Ratios vs. a0 ======
% We define policy functions on matching grids.

grid_opts   = linspace(S.a, S.b, length(S.opts));
grid_dyncon = linspace(S.a, S.b, length(S.dyncon));
grid_naive  = linspace(S.a, S.b, length(S.naive));

g_markov = @(a_val) interp1(grid_opts, S.opts, a_val, 'linear', 'extrap');
g_dyn    = @(a_val) interp1(grid_dyncon, S.dyncon, a_val, 'linear', 'extrap');
g_naive  = @(a_val) interp1(grid_naive, S.naive, a_val, 'linear', 'extrap');

n_a0 = S.n + 1;
a0_vals = linspace(S.a, S.b, n_a0);

% Pre-allocate arrays for lifetime utility computations.
perc_markov = zeros(1, n_a0);
perc_naive  = zeros(1, n_a0);
perc_hybrid = zeros(1, n_a0);

V_markov_arr = zeros(1, n_a0); 
V_dc_arr     = zeros(1, n_a0); 
V_naive_arr  = zeros(1, n_a0);
U_markov_arr = zeros(1, n_a0);
U_dc_arr     = zeros(1, n_a0);
U_naive_arr  = zeros(1, n_a0);
U_hybrid_arr = zeros(1, n_a0);

max_iter = 150;


for i = 1:n_a0
    a0 = a0_vals(i);
     
    % --- Markov Policy ---
    a1 = g_markov(a0);
    c0 = S.f(a0) - a1;  % instead of (1+S.r)*a0 + S.W - a1
    U0_markov = S.u(c0);
    U_markov = U0_markov;
    a_cur = a1;
    for iter=1:max_iter
        a_next = g_markov(a_cur);
        c = S.f(a_cur) - a_next;
        term = beta * (delta^iter) * S.u(c);
        U_markov = U_markov + term;
        a_cur = a_next;        
    end
    U_markov_arr(i) = U_markov;
    V_markov = (U_markov - U0_markov) / (beta * delta);
    V_markov_arr(i) = V_markov;

    % --- DC Policy ---
    a1 = g_dyn(a0);
    c0 = S.f(a0) - a1;
    U0_dc = S.u(c0);
    U_dc = U0_dc;
    a_cur = a1;
    for iter=1:max_iter
        a_next = g_dyn(a_cur);
        c = S.f(a_cur) - a_next;
        term = beta * (delta^iter) * S.u(c);
        U_dc = U_dc + term;
        a_cur = a_next; 
    end
    U_dc_arr(i) = U_dc;
    V_dc = (U_dc - U0_dc) / (beta * delta);
    V_dc_arr(i) = V_dc;

    % --- Naive Policy ---
    a1 = g_naive(a0);
    c0 = S.f(a0) - a1;
    U0_naive = S.u(c0);
    U_naive = U0_naive;
    a_cur = a1;
    for iter=1:max_iter
        a_next = g_naive(a_cur);
        c = S.f(a_cur) - a_next;
        term = beta * (delta^iter) * S.u(c);
        U_naive = U_naive + term;
        a_cur = a_next; 
    end
    U_naive_arr(i) = U_naive;
    V_naive = (U_naive - U0_naive) / (beta * delta);
    V_naive_arr(i) = V_naive;

    % --- Hybrid Policy: one-step naive, then DC ---
    a1 = g_naive(a0);
    cons0 = S.f(a0) - a1;
    U_hybrid = S.u(cons0);
    a_cur = a1;
    iter = 0;
    while iter <= max_iter
        a_next = g_dyn(a_cur);
        c = S.f(a_cur) - a_next;
        term = beta * (delta^iter) * S.u(c);  
        U_hybrid = U_hybrid + term;
        a_cur = a_next;
        iter = iter + 1;
    end
    U_hybrid_arr(i) = U_hybrid;
    
    % Compare versus the hybrid (full commitment) solution
    perc_markov(i) = 100 * (U_markov / U_hybrid);
    perc_naive(i)  = 100 * (U_naive  / U_hybrid);
    perc_hybrid(i) = 100 * (U_hybrid / U_hybrid);
end
 
%% Plot Figure 2
figure(2); clf;
displayNameDC = sprintf('Dynamically Consistent');
displayNameN = sprintf('Na\\"{i}ve Solution');
displayNameMP = sprintf('Markov Equilibrium');
displayNameFC = sprintf('Full Commitment');
plot(a0_vals, perc_markov, '-',  'LineWidth',3, 'Color', nred);  hold on;
plot(a0_vals, perc_naive,  '--', 'LineWidth',3, 'Color', ngreen);
xlabel('Initial Assets $x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Disc. Util. in \% of Full Commitment Sol.', 'Interpreter', 'latex', 'FontSize', 13);
legend(displayNameMP, displayNameN); 
legend('location', 'southwest','FontSize',16, 'Interpreter', 'latex');
legend boxoff;
xlim([S.a, S.bspecfigr]);
drawnow;

%% FIGURE 3: Absolute Utility Values
figure(3); clf;
p_dc_abs = plot(a0_vals, U_dc_arr, '-', 'LineWidth',3, 'Color', nblue, 'DisplayName', 'Dynamically Consistent');
hold on;
p_markov_abs = plot(a0_vals, U_markov_arr, '-', 'LineWidth',3, 'Color', nred, 'DisplayName', 'Markov Equilibrium');
p_naive_abs = plot(a0_vals, U_naive_arr, '--', 'LineWidth',3, 'Color', ngreen, 'DisplayName', sprintf('Na\\"{i}ve Solution'));
hold off;
xlabel('Initial Assets', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Utility', 'Interpreter', 'latex', 'FontSize', 16);
legend([p_dc_abs, p_markov_abs, p_naive_abs], 'Location','southeast','FontSize',16, 'Interpreter','latex');
legend boxoff;
xlim([S.a, S.bspecfigr]);

%% FIGURE 4: Consumption Sequences vs. Time
max_tau = 300;
tol_steady = 1e-8;
steady_periods = 3;
a_markov = S.b;
a_dyn = S.b;
a_naive = S.b;
a_hybrid = S.b;
cons_markov = zeros(1, max_tau);
cons_dyn = zeros(1, max_tau);
cons_naive = zeros(1, max_tau);
cons_hybrid = zeros(1, max_tau);

for tau = 1:max_tau
    next_a = g_markov(a_markov);
    cons_markov(tau) = S.f(a_markov) - next_a;
    a_markov = next_a;
    
    next_a = g_dyn(a_dyn);
    cons_dyn(tau) = S.f(a_dyn) - next_a;
    a_dyn = next_a;
    
    next_a = g_naive(a_naive);
    cons_naive(tau) = S.f(a_naive) - next_a;
    a_naive = next_a;
    
    if tau == 1
        next_a = g_naive(a_hybrid);
    else
        next_a = g_dyn(a_hybrid);
    end
    cons_hybrid(tau) = S.f(a_hybrid) - next_a;
    a_hybrid = next_a;
    
    if tau >= steady_periods
        steady_markov = max(cons_markov(tau-steady_periods+1:tau)) - min(cons_markov(tau-steady_periods+1:tau)) < tol_steady;
        steady_dyn    = max(cons_dyn(tau-steady_periods+1:tau))    - min(cons_dyn(tau-steady_periods+1:tau))    < tol_steady;
        steady_naive  = max(cons_naive(tau-steady_periods+1:tau))  - min(cons_naive(tau-steady_periods+1:tau))  < tol_steady;
        steady_hybrid = max(cons_hybrid(tau-steady_periods+1:tau)) - min(cons_hybrid(tau-steady_periods+1:tau)) < tol_steady;
        
        if steady_markov && steady_dyn && steady_naive && steady_hybrid
            break;
        end
    end
end

time_shifted = 0:(tau-1);
cons_markov = cons_markov(1:tau);
cons_dyn    = cons_dyn(1:tau);
cons_naive  = cons_naive(1:tau);
cons_hybrid = cons_hybrid(1:tau);

figure(4); clf;
plot(time_shifted, cons_markov, '-',  'LineWidth',2.5, 'Color', nred);   hold on;
plot(time_shifted, cons_dyn,    '-', 'LineWidth',2.5, 'Color', nblue);
plot(time_shifted, cons_naive,  '--', 'LineWidth',2.5, 'Color', ngreen);
plot(time_shifted, cons_hybrid, ':',  'LineWidth',2.5, 'Color', [0 0 0]);
xlabel('Time','FontSize',14);
ylabel('Consumption','FontSize',14);
legend('Markov Policy', 'Dynamically Consistent', 'Naïve Policy', 'Full Commitment','Location','Best');
grid on;
legend boxoff;

%% ====== FIGURE 5:  
displayNameDC = sprintf('Dynamically Consistent');
displayNameN = sprintf('Na\\"{i}ve Solution');
displayNameMP = sprintf('Markov Policy');
N_val = length(S.opts);
grid_opts = linspace(S.a, S.b, N_val);
grid_dyncon = linspace(S.a, S.b, N_val);
grid_naive = linspace(S.a, S.b, N_val);
g_markov = @(a_val) interp1(grid_opts, S.opts, a_val, 'linear', 'extrap');
g_dyn = @(a_val) interp1(grid_dyncon, S.dyncon, a_val, 'linear', 'extrap');
g_naive = @(a_val) interp1(grid_naive, S.naive, a_val, 'linear', 'extrap');
n_points = 100;
a_vals = linspace(S.a, S.b, n_points);
V_markov_sum = zeros(size(a_vals));
V_naive_sum  = zeros(size(a_vals));
V_dyn_sum    = zeros(size(a_vals));

for i = 1:length(a_vals)
    a0 = a_vals(i);
    a_current = a0;
    sum_val = 0;
    for t = 0:cutoff
        a_next = g_markov(a_current);
        c_t = S.f(a_current) - a_next;
        sum_val = sum_val + (delta^t)*S.u(c_t);
        a_current = a_next;
    end
    V_markov_sum(i) = sum_val;
    
    a_current = a0;
    sum_val = 0;
    for t = 0:cutoff
        a_next = g_naive(a_current);
        c_t = S.f(a_current) - a_next;
        sum_val = sum_val + (delta^t)*S.u(c_t);
        a_current = a_next;
    end
    V_naive_sum(i) = sum_val;
    
    a_current = a0;
    sum_val = 0;
    for t = 0:cutoff
        a_next = g_dyn(a_current);
        c_t = S.f(a_current) - a_next;
        sum_val = sum_val + (delta^t)*S.u(c_t);
        a_current = a_next;
    end
    V_dyn_sum(i) = sum_val;
end
 
figure(5); clf;
plot(a_vals, V_dyn_sum, 'LineWidth',3, 'Color', nblue, 'DisplayName', displayNameDC);
hold on;
plot(a_vals, V_markov_sum, 'LineWidth',3, 'Color', nred, 'DisplayName', displayNameMP);
hold off;
xlabel('Saving', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Value Function ($W$)', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location','southeast','FontSize',16,'Interpreter','latex'); 
legend boxoff;
xlim([S.a, S.bspecfigcv]);
end
