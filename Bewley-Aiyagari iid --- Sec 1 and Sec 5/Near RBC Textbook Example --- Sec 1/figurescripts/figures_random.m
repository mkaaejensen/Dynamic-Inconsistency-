function figures_random(S)
% All workspace variables are contained in the structure S.
% Every reference uses S.<variable>.

%% Define local colors
str = '#000080';
nblue = sscanf(str(2:end),'%2x%2x%2x',[1 3]) / 255;
ngreen = [0 0.5 0];
str = '#800000';
nred = sscanf(str(2:end),'%2x%2x%2x',[1 3]) / 255;

%% Various additional calculations

% Policy functions for forward simulations
policy_QH_low = @(x) interp1(linspace(S.a, S.b, length(S.opts)), S.opts, x, 'linear', 'extrap');
policy_QH_high = @(x) interp1(linspace(S.a, S.b, length(S.opts)), S.opts, S.offset(x), 'linear', 'extrap');

policy_A_low = @(x) interp1(linspace(S.a, S.b, length(S.dyncon)), S.dyncon, x, 'linear', 'extrap');
policy_A_high = @(x) interp1(linspace(S.a, S.b, length(S.dyncon)), S.dyncon, S.offset(x), 'linear', 'extrap');

policy_N_low = @(x) interp1(linspace(S.a, S.b, length(S.naive)), S.naive, x, 'linear', 'extrap');
policy_N_high = @(x) interp1(linspace(S.a, S.b, length(S.naive)), S.naive, S.offset(x), 'linear', 'extrap');

Nsim = 50000;
aysim  = S.asim(:, S.Tsim)  ./ mean(S.ysim(:, S.Tsim));
aysimK = S.asimK(:, S.TsimK) ./ mean(S.ysimK(:, S.TsimK));
aysimKap = S.asimKap(:, S.TsimKap) ./ mean(S.ysimK(:, S.TsimKap));
quantf  = @(p) quantile(aysim, p);
quantfK = @(p) quantile(aysimK, p);
MEAN  = integral(quantf, 0, 1);
MEANK = integral(quantfK, 0, 1);
Lorenz  = @(p) (integral(quantf, 0, p)) ./ MEAN;
LorenzK = @(p) (integral(quantfK, 0, p)) ./ MEANK;
j = 1;
for t1 = 0:0.01:1
    LP(j)  = Lorenz(t1);
    LPK(j) = LorenzK(t1);
    j = j + 1;
end
 
 

pgrid = 0:0.01:1;

% Gini coefficients via the trapezoidal rule
Gini_QH = 1 - 2 * trapz(pgrid, LP);
Gini_A  = 1 - 2 * trapz(pgrid, LPK);

%% Display statistics
 

drawnow

%% Rounding for figures
S.MA = round(S.MA, 1);
S.MAK = round(S.MAK, 1);

%% FIGURE 1: Policy functions
figure(1); clf;

displayNameDC = sprintf('Dynamically consistent ($\\beta=1$)');
displayNameDCH = sprintf('Dynamically Consistent ($l_t=\\overline{l}$)');
displayNameN = sprintf('Na\\"{i}ve Solution ($\\beta = %s$)', num2str(S.beta));
displayNameMP = sprintf('Quasi-hyperbolic ($\\beta = %s$)', num2str(S.beta));
displayNameMPH = sprintf('Quasi-hyperbolic ($l_t=\\overline{l}$)');
displayNameSM = sprintf('Q-H, Second-mover ($\\beta = %s$)', num2str(S.beta));

x1 = linspace(S.a, S.b, S.n+1);

optsH = policy_QH_high(x1);
dynconH = policy_A_high(x1);
naiveH = policy_N_high(x1);

fplot(@(j) j, [S.a, S.b], ':', 'linewidth', 1, 'Color', uint8([5 5 5]), 'HandleVisibility', 'off')
hold on

p1 = plot(x1, S.dyncon, '-', 'color', nblue, 'LineWidth', 3, 'DisplayName', displayNameDC);

if S.shownaive == 1
    p2 = plot(x1, S.naive, '--', 'LineWidth', 3, 'Color', "k", 'DisplayName', displayNameN);
end

p3 = plot(x1, S.opts, '-', 'LineWidth', 3, 'color', nred, 'DisplayName', displayNameMP);
%p4 = plot(x1, optsH, '-', 'LineWidth', 3, 'color', nred, 'DisplayName', displayNameMPH);
%p5 = plot(x1, dynconH, '-', 'color', "k", 'LineWidth', 3);
%set(p5, 'HandleVisibility', 'off');

if S.shownaive == 0
    legend([p1, p3], 'location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');
else
    legend([p1, p3, p2], 'location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');
end
if S.beta == 1
    plot(S.grid_kaplan, S.dyncon_kaplan, 'y--', 'LineWidth', 1.5, 'DisplayName', 'Dyn. Con.(EGP Algor.)')
end
legend boxoff
xlabel('Assets at date $t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Assets at date $t+1$', 'Interpreter', 'latex', 'FontSize', 16);
xlim([S.a, S.bspecfig]);
hold off
drawnow

if S.shortv == 1
    disp(['Mean assets QH: ' num2str(S.MA)]);
    disp(['Mean assets Aiyagari: ' num2str(S.MAK)]);
    return;
end
 

figure(2); clf;
fplot(@(j) j, [S.a, S.b], ':', 'linewidth', 0.5, 'Color', uint8([5 5 5]), 'HandleVisibility', 'off')
hold on
x = S.a : ((S.b-S.a)/S.n) : S.b;
[~, Br] = Lossfunction(S.a, S.b, S.n, S.opts, S.dissaving, S.dis, S.u, S.p, 0, S.tol, S.f, S.offset, S.dsp);
ids = transpose(Br);
plot(x, ids, '-', 'LineWidth', 2.5, 'color', nblue, 'DisplayName', 'Second-Mover');
hold on
x1 = S.a : ((S.b-S.a)/S.n) : S.b;
t1 = 1 : 1 : S.n+1;
y1 = S.opts(t1);
plot(x1, y1, '-', 'LineWidth', 2.5, 'color', nred, 'DisplayName', 'Markov Policy')
plot(x, S.dyncon, '-', 'color', nblue, 'LineWidth', 2.5, 'DisplayName', 'Dynamically consistent');
plot(x, S.naive, '--', 'LineWidth', 2.5, 'color', ngreen, 'DisplayName', 'Naive solution');
if S.showkaplan == 1
    plot(S.grid_kaplan, S.dyncon_kaplan, 'y--', 'LineWidth', 1.5, 'DisplayName', 'Dynamically consistent (EGP Algor.)')
end
xlim([S.a S.b]);
legend('location', 'northwest')
legend boxoff
hold off
 
%% FIGURE 3: Policy functions (alt)
figure(3); clf;
 



fplot(@(j) j, [S.a, S.b], ':', 'linewidth', 1, 'Color', uint8([5 5 5]), 'HandleVisibility', 'off')
hold on

p1 = plot(x1, S.dyncon, '-', 'color', nblue, 'LineWidth', 3, 'DisplayName', displayNameDC);

if S.shownaive == 1
    p2 = plot(x1, S.naive, '--', 'LineWidth', 3, 'Color', "k", 'DisplayName', displayNameN);
end

p3 = plot(x, ids, '-', 'LineWidth', 3, 'color', nred, 'DisplayName', displayNameSM);
%p4 = plot(x1, optsH, '-', 'LineWidth', 3, 'color', nred, 'DisplayName', displayNameMPH);
%p5 = plot(x1, dynconH, '-', 'color', "k", 'LineWidth', 3);
%set(p5, 'HandleVisibility', 'off');

if S.shownaive == 0
    legend([p1, p3], 'location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');
else
    legend([p1, p3, p2], 'location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');
end
if S.beta == 1
    plot(S.grid_kaplan, S.dyncon_kaplan, 'y--', 'LineWidth', 1.5, 'DisplayName', 'Dyn. cons. (EGP Algor.)')
end
legend boxoff
xlabel('Assets at date $t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Assets at date $t+1$', 'Interpreter', 'latex', 'FontSize', 16);
xlim([S.a, S.bspecfig]);
hold off
drawnow