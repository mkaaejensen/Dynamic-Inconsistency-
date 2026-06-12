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
% Create appropriate grid based on power_param
FigGrid = S.G;

% Extract policies for low and high states
opts_low = S.opts(1:S.n+1);
opts_high = S.opts(S.n+2:end);
dyncon_low = S.dyncon(1:S.n+1);
dyncon_high = S.dyncon(S.n+2:end);

% Create policy function interpolation for both states
policy_QH_low = @(x) interp1(FigGrid, opts_low, x, 'linear', 'extrap');
policy_QH_high = @(x) interp1(FigGrid, opts_high, x, 'linear', 'extrap');

policy_A_low = @(x) interp1(FigGrid, dyncon_low, x, 'linear', 'extrap');
policy_A_high = @(x) interp1(FigGrid, dyncon_high, x, 'linear', 'extrap');

% Also handle naive policy if available
if isfield(S, 'naive') && ~isempty(S.naive)
    naive_low = S.naive(1:S.n+1);
    naive_high = S.naive(S.n+2:end);
    policy_N_low = @(x) interp1(FigGrid, naive_low, x, 'linear', 'extrap');
    policy_N_high = @(x) interp1(FigGrid, naive_high, x, 'linear', 'extrap');
end

if S.shortv == 1
    disp('Key asset distribution statistics (means only because shortv = 1)');
else
    disp('Key asset distribution statistics');
end

% Means
MEAN  = S.MA;
MEANK = S.MAK;
 
disp(['Mean assets QH: ' num2str(S.MA) ' | Aiyagari: ' num2str(S.MAK)]);
 
 
%% Rounding for figures
S.MA = round(S.MA, 1);
S.MAK = round(S.MAK, 1);
 
%% FIGURE 1: Policy functions
figure(1); clf;

displayNameDC = sprintf('Aiyagari ($l_t\\in \\{ \\underline{l},\\overline{l} \\}$)');
displayNameDCH = sprintf('Dynamically Consistent ($l_t=\\overline{l}$)');
displayNameN = sprintf('Na\\"{i}ve Solution ($\\beta = %s$)', num2str(S.beta));
displayNameMP = sprintf('Quasi-hyperbolic ($l_t=\\underline{l}$)');
displayNameMPH = sprintf('Quasi-hyperbolic ($l_t=\\overline{l}$)');

fplot(@(j) j, [S.a, S.b], ':', 'linewidth', 1, 'Color', uint8([5 5 5]), 'HandleVisibility', 'off')
hold on

p1 = plot(FigGrid, dyncon_low, '-', 'color', "k", 'LineWidth', 3, 'DisplayName', displayNameDC);
p5 = plot(FigGrid, dyncon_high, '-', 'color', "k", 'LineWidth', 3);
set(p5, 'HandleVisibility', 'off');

if S.shownaive == 1 && isfield(S, 'naive') && ~isempty(S.naive)
    p2 = plot(FigGrid, naive_low, '--', 'LineWidth', 3, 'Color', "k", 'DisplayName', displayNameN);
end

p3 = plot(FigGrid, opts_low, '-', 'LineWidth', 3, 'color', ngreen, 'DisplayName', displayNameMP);
p4 = plot(FigGrid, opts_high, '-', 'LineWidth', 3, 'color', nred, 'DisplayName', displayNameMPH);
 
if S.shownaive == 0
    legend([p1, p3, p4], 'location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');
else
    legend([p1, p3, p2], 'location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');
end

if S.beta == 1 && isfield(S, 'dyncon_kaplan') && isfield(S, 'grid_kaplan')
    plot(S.grid_kaplan, S.dyncon_kaplan(1:length(S.grid_kaplan)), 'y--', 'LineWidth', 1.5, 'DisplayName', 'Dyn. cons. (Kaplan)')
end

legend boxoff
xlabel('Assets at date $t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Assets at date $t+1$', 'Interpreter', 'latex', 'FontSize', 16);
xlim([S.a, S.bspecfig]);
hold off
drawnow

if S.shortv == 1
    return;
end

Nsim = 50000;
if isfield(S, 'asimKap') && isfield(S, 'lsimKap') && isfield(S, 'TsimKap')
    asimKap = S.asimKap;
end

 

% Use last period cross section for distributional statistics
asim = S.asim(:,S.Tsim);
lsim = S.lsim(:,S.Tsim);
asimK  = S.asimK(:, S.TsimK);
lsimK = S.lsimK(:, S.TsimK);
quantf  = @(p) quantile(asim, p);
quantfK = @(p) quantile(asimK, p);

 

Lorenz  = @(p) (integral(quantf, 0, p)) ./ mean(asim(:));
LorenzK = @(p) (integral(quantfK, 0, p)) ./ mean(asimK(:));
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

%% Display full statistics
disp(['Median assets QH: ' num2str(quantile(asim, .5)) ' | Aiyagari: ' num2str(quantile(asimK, .5))]);
disp(['Skewness QH: ' num2str(skewness(asim)) ' | Aiyagari: ' num2str(skewness(asimK))]);
disp(['Gini coefficient QH: ' num2str(Gini_QH) ' | Aiyagari: ' num2str(Gini_A)]);
disp(['Fraction borrowing constrained QH: ' num2str(sum(asim <= S.a) ./ Nsim * 100) '% | Aiyagari: ' num2str(sum(asimK <= S.a) ./ Nsim * 100) '%']);

borrow_lim=S.a;
if S.extrads == 1
    % -------------------------------------------------------------------------
% Wealth-distribution statistics: QH vs. Aiyagari benchmark
% -------------------------------------------------------------------------

sorted_QH = sort(asim);
sorted_A  = sort(asimK);
total_QH  = sum(sorted_QH);
total_A   = sum(sorted_A);

% ---------------------------------------------------------------
% Cumulative population & wealth shares
% ---------------------------------------------------------------
cum_pop       = (1:Nsim)'/Nsim;          % share of population
cum_wealth_QH = cumsum(sorted_QH)/total_QH;
cum_wealth_A  = cumsum(sorted_A) /total_A;

% ---------------------------------------------------------------
% Inverse Lorenz shares: population share whose cumulative wealth <= X%% of total
% (computed by inverting pre‑computed Lorenz curves LP, LPK evaluated on pgrid)
% ---------------------------------------------------------------
wealth_cutoffs = [0.001 0.01 0.02 0.05 0.10 0.20 0.30 0.40]; % 0.1,1,2,5,10,20,30,40 percent

% Ensure strictly increasing Lorenz ordinates before inversion
[uLP , idxLP ]  = unique(LP ,'stable');
[uLPK, idxLPK]  = unique(LPK,'stable');

invLorenz_QH = interp1(uLP ,  pgrid(idxLP ),  wealth_cutoffs, 'linear', 'extrap');
invLorenz_A  = interp1(uLPK,  pgrid(idxLPK),  wealth_cutoffs, 'linear', 'extrap');

% ---------------------------------------------------------------
% Top & bottom shares
% ---------------------------------------------------------------
bottom_shares_QH = interp1(cum_pop,cum_wealth_QH,[0.01 0.05 0.10 0.20 0.40 0.50]);
bottom_shares_A  = interp1(cum_pop,cum_wealth_A ,[0.01 0.05 0.10 0.20 0.40 0.50]);
top_shares_QH    = 1 - interp1(cum_pop,cum_wealth_QH,[0.999 0.99 0.95 0.90]);
top_shares_A     = 1 - interp1(cum_pop,cum_wealth_A ,[0.999 0.99 0.95 0.90]);

% ---------------------------------------------------------------
% Decile shares
% ---------------------------------------------------------------
Decile_QH = zeros(1,10);
Decile_A  = zeros(1,10);
for d = 1:10
    i1 = floor((d-1)*Nsim/10)+1;
    i2 = floor( d   *Nsim/10);
    Decile_QH(d) = sum(sorted_QH(i1:i2))/total_QH*100;
    Decile_A(d)  = sum(sorted_A (i1:i2))/total_A *100;
end

% ---------------------------------------------------------------
% Basic inequality indices
% ---------------------------------------------------------------
mQ    = mean(asim);    mA  = mean(asimK);
maskQ = asim  > 0;     maskA = asimK > 0;

Theil_QH = mean((asim(maskQ)./mQ).*log(asim(maskQ)./mQ));
Theil_A  = mean((asimK(maskA)./mA).*log(asimK(maskA)./mA));

epsA  = 0.5;
termQ = mean((asim ./mQ).^(1-epsA));
A_QH  = 1 - termQ^(1/(1-epsA));
termA = mean((asimK./mA).^(1-epsA));
A_Aiyagari = 1 - termA^(1/(1-epsA));

medQ = median(asim);   medA = median(asimK);
med_mean_QH = medQ/mQ;
med_mean_A  = medA/mA;

cv_QH = std(asim)/mQ;
cv_A  = std(asimK)/mA;

% ---------------------------------------------------------------
% Polarization & concentration
% ---------------------------------------------------------------
Palma_QH = top_shares_QH(4)/bottom_shares_QH(5);   % Top10 / Bot40
Palma_A  = top_shares_A (4)/bottom_shares_A (5);

top01_QH = top_shares_QH(1)*100;   % Top 0.1%
top01_A  = top_shares_A (1)*100;

constraint_mass_QH = sum(asim <= S.a + 1e-6)/Nsim*100;
constraint_mass_A  = sum(asimK<= S.a + 1e-6)/Nsim*100;

% ---------------------------------------------------------------
% === ORGANIZED OUTPUT ===
% ---------------------------------------------------------------
fprintf('\n=== WEALTH DISTRIBUTION STATISTICS ===\n');

% Basic inequality measures
fprintf('\n--- Basic Measures ---\n');
fprintf('Gini coefficient        QH: %.3f | A: %.3f\n', Gini_QH, Gini_A);
fprintf('Theil index            QH: %.3f | A: %.3f\n', Theil_QH, Theil_A);
fprintf('Atkinson (eps=0.5)     QH: %.3f | A: %.3f\n', A_QH, A_Aiyagari);
fprintf('Median/Mean ratio      QH: %.3f | A: %.3f\n', med_mean_QH, med_mean_A);
fprintf('Coeff. of variation    QH: %.3f | A: %.3f\n', cv_QH, cv_A);

% Wealth shares
fprintf('\n--- Wealth Shares ---\n');
fprintf('Top 0.1%%               QH: %4.1f%% | A: %4.1f%%\n',top01_QH,top01_A);
fprintf('Top 1%%                 QH: %4.1f%% | A: %4.1f%%\n',100*top_shares_QH(2),100*top_shares_A(2));
fprintf('Top 5%%                 QH: %4.1f%% | A: %4.1f%%\n',100*top_shares_QH(3),100*top_shares_A(3));
fprintf('Top 10%%                QH: %4.1f%% | A: %4.1f%%\n',100*top_shares_QH(4),100*top_shares_A(4));
fprintf('Bottom 50%%             QH: %4.1f%% | A: %4.1f%%\n',100*bottom_shares_QH(6),100*bottom_shares_A(6));
fprintf('Bottom 40%%             QH: %4.1f%% | A: %4.1f%%\n',100*bottom_shares_QH(5),100*bottom_shares_A(5));

% Decile shares (compact, two lines each)
fprintf('\n--- Decile Shares ---\n');
fprintf('QH:   ');
for d = 1:10
    fprintf('D%d: %4.1f%%  ',d,Decile_QH(d));
    if d==5, fprintf('\n      '); end
end
fprintf('\nA:    ');
for d = 1:10
    fprintf('D%d: %4.1f%%  ',d,Decile_A(d));
    if d==5, fprintf('\n      '); end
end
 

% Tail measures
fprintf('\n--- Tail Measures ---\n');
fprintf('Palma ratio            QH: %.2f | A: %.2f\n',Palma_QH,Palma_A);
fprintf('Borrowing constrained  QH: %.1f%% | A: %.1f%%\n',constraint_mass_QH,constraint_mass_A);

end


drawnow

%% FIGURE 2: Asset distributions
if isrow(S.opts)
    S.opts = S.opts'; % transform to improve speed
end

%% FIGURE 2 (asset distributions)

figure(2); clf;
bin_width = 1;

% Define bins
edges = S.a : bin_width : S.b;

% Get raw histcounts (probabilities) including the first bin at a=0
histQH = histcounts(S.asim(:, S.Tsim),  edges, 'Normalization','probability');
histA  = histcounts(S.asimK(:, S.TsimK), edges, 'Normalization','probability');

% Extract the point‐mass at the borrowing limit (first bin)
massQH0 = histQH(1);
massA0  = histA(1);

hold on;



% 2) “Stairs” histogram for the continuous part (assets > a)
h1 = histogram(...
    S.asim(S.asim(:,S.Tsim)>=S.a, S.Tsim), ...
    edges, 'Normalization','probability', 'DisplayStyle','stairs');
set(h1, 'EdgeColor', nred,  'LineWidth', 3);

h2 = histogram(...
    S.asimK(S.asimK(:,S.TsimK)>=S.a, S.TsimK), ...
    edges, 'Normalization','probability', 'DisplayStyle','stairs');
set(h2, 'EdgeColor', nblue, 'LineWidth', 3);

% Legend dummies
dummy1 = plot(nan, nan, 'Color', nred,  'LineWidth', 2.5);
dummy2 = plot(nan, nan, 'Color', nblue, 'LineWidth', 2.5);
legend([dummy1,dummy2], ...
       {['Quasi-hyp. (mean ' num2str(S.MA) ')'], ['Aiyagari (mean ' num2str(S.MAK) ')']}, ...
       'Location','northeast', 'FontSize', 16, 'Interpreter','latex');
legend boxoff;

% Axes limits
x_max = max([ max(h1.BinEdges(h1.Values>1e-4)), max(h2.BinEdges(h2.Values>1e-4)) ]);
x_min = min([ S.asim(:); S.asimK(:) ]);
if isnan(x_max) || x_max <= x_min
    x_max = x_min + 1;
end
xlim([x_min, x_max]);
ylim([0, max([h1.Values, h2.Values]) * 1.1]);

% Labels & formatting
xlabel('Assets', 'Interpreter','latex', 'FontSize', 16);
ylabel('Fraction of households', 'Interpreter','latex', 'FontSize', 16);
box off;
hold off;

 
%% FIGURE 3 (asset distributions with discrete mass at borrowing limit)
figure(3); clf;
bw   = 0.25;           % bin width
bw_k = 1;              % KDE bandwidth (increase for more smoothness)
a    = S.a;

dQ  = S.asim(:, S.Tsim);   dA  = S.asimK(:, S.TsimK);
mQ  = mean(dQ <= a);       mA  = mean(dA <= a);
cQ  = dQ(dQ > a);          cA  = dA(dA > a);

% create xi starting exactly at a
xi  = linspace(a, max([cQ; cA]), 2000);
fQ  = ksdensity(cQ, xi, 'Bandwidth', bw_k) * (1 - mQ) * bw;
fA  = ksdensity(cA, xi, 'Bandwidth', bw_k) * (1 - mA) * bw;
% clamp the first point to the spike mass
fQ(1) = mQ;
fA(1) = mA;

hold on
plot(xi, fQ, 'Color', nred,  'LineWidth', 3);
plot(xi, fA, 'Color', nblue, 'LineWidth', 3);
% center dots at a
plot(a, mQ, 'o', 'MarkerSize', 8, 'MarkerFaceColor', nred,  'MarkerEdgeColor', nred);
plot(a, mA, 'o', 'MarkerSize', 8, 'MarkerFaceColor', nblue, 'MarkerEdgeColor', nblue);

xlim([a, xi(end)]);
ylim([0, max([mQ, mA, fQ, fA]) * 1.1]);

xlabel('Assets','Interpreter','latex','FontSize',16);
ylabel('Fraction of households','Interpreter','latex','FontSize',16);

d1 = plot(nan,nan,'Color',nred,'LineWidth',2.5);
d2 = plot(nan,nan,'Color',nblue,'LineWidth',2.5);
legend([d1,d2], ...
    {['Quasi-hyp. (mean ' num2str(S.MA) ')'], ['Aiyagari (mean ' num2str(S.MAK) ')']}, ...
    'Location','northeast','FontSize',16,'Interpreter','latex');
legend boxoff; box off; hold off;


% === FIGURE 4: kernel estimators
figure(4); clf;
bw   = 0.25;        % bin-width factor  
bw_k = 1;           % KDE bandwidth
a    = S.a;         % borrowing limit

dQ = S.asim (:, S.Tsim );          dA  = S.asimK(:, S.TsimK);
cQ = dQ(dQ >= a);                  cA  = dA(dA >= a);

xi = linspace(a, max([cQ; cA]), 2000);

fQ = ksdensity(cQ, xi, 'Bandwidth', bw_k, 'BoundaryCorrection','reflection') * bw;
fA = ksdensity(cA, xi, 'Bandwidth', bw_k, 'BoundaryCorrection','reflection') * bw;

hold on
plot(xi, fQ, 'Color', nred , 'LineWidth', 3);
plot(xi, fA, 'Color', nblue, 'LineWidth', 3);

xlim([a, xi(end)]);
ylim([0, max([fQ, fA])*1.1]);

xlabel('Assets',                'Interpreter','latex','FontSize',16);
ylabel('Density','Interpreter','latex','FontSize',16);

d1 = plot(nan,nan,'Color',nred , 'LineWidth',2.5);
d2 = plot(nan,nan,'Color',nblue, 'LineWidth',2.5);
legend([d1,d2], ...
       {['Quasi-hyp. (Mean ' num2str(S.MA ) ')'], ...
        ['Aiyagari (Mean ' num2str(S.MAK) ')']}, ...
       'Location','northeast','FontSize',16,'Interpreter','latex');
legend boxoff; box off; hold off;

%% main figure 4 (comment out to obtain kernel estimate instead)
% figure(4); clf;
% bin_width=1;
%edges = S.a : bin_width : S.b;
%countsQH = histcounts(S.asim(:,S.Tsim), edges, 'Normalization','probability');
%countsA = histcounts(S.asimK(:,S.TsimK), edges, 'Normalization','probability');
%xLeft = edges(1:end-1);
%threshold = 0.003;
%mask = (countsQH >= threshold) | (countsA >= threshold);
%lastIdx = find(mask, 1, 'last');
%xMax = xLeft(lastIdx);
%windowSize = 3;
%countsQH_smooth = movmean(countsQH, windowSize);
%countsA_smooth = movmean(countsA, windowSize);

%xx = linspace(xLeft(1), xMax, 2000);
%yyQH = interp1(xLeft, countsQH_smooth, xx, 'spline');
%yyA = interp1(xLeft, countsA_smooth, xx, 'spline');
%plot(xx, yyQH, '-', 'Color', nred, 'LineWidth', 3); hold on;
%plot(xx, yyA, '-', 'Color', nblue, 'LineWidth', 3);
%xlim([S.a, xMax]);
%ylim([0, max([yyQH, yyA]) * 1.1]);
%xlabel('Assets', 'Interpreter', 'latex', 'FontSize', 16);
%ylabel('Fraction of households', 'Interpreter', 'latex', 'FontSize', 16);
%dummy1 = plot(nan, nan, 'Color', nred, 'LineWidth', 2.5);
%dummy2 = plot(nan, nan, 'Color', nblue, 'LineWidth', 2.5);
%legend([dummy1, dummy2], {['Quasi-hyp. (mean ' num2str(S.MA) ')'], ['Aiyagari (mean ' num2str(S.MAK) ')']}, 'Location', 'northeast', 'FontSize', 16, 'Interpreter', 'latex');
%legend boxoff;
%box off;



figure(5); clf;
fplot(@(x) x, [S.a, S.b], ':', 'LineWidth', 1, 'Color', uint8([5 5 5]), 'HandleVisibility', 'off')
hold on
h1 = plot(pgrid, LP, '-', 'color', nred, 'LineWidth', 2.5);
h2 = plot(pgrid, LPK, '-', 'color', nblue, 'LineWidth', 2.5);
label1 = sprintf('Quasi-Hyperbolic (Gini: %.2f)', Gini_QH);
label2 = sprintf('Aiyagari (Gini: %.2f)', Gini_A);
legend([h1, h2], label1, label2, 'Interpreter', 'latex', 'fontsize', 16);
xlim([0 1]);
xlabel('Cumulative Share of Individuals', 'Interpreter', 'latex', 'fontsize', 14)
ylabel('Cumulative Share of Assets', 'Interpreter', 'latex', 'fontsize', 14)
legend('location', 'northwest')
legend boxoff
hold off


figure(6); clf;
fplot(@(j) j, [S.a, S.b], ':', 'linewidth', 0.5, 'Color', uint8([5 5 5]), 'HandleVisibility', 'off')
hold on

% Get best response - modified for two-state Markov
Br=S.Br;
Br_low = Br(1:S.n+1);
Br_high = Br(S.n+2:end);

plot(FigGrid, Br_low, '--', 'LineWidth', 1.5, 'color', nblue, 'DisplayName', 'Second-Mover (low state)');
plot(FigGrid, Br_high, '--', 'LineWidth', 1.5, 'color', ngreen, 'DisplayName', 'Second-Mover (high state)');
hold on
plot(FigGrid, opts_low, '-', 'LineWidth', 2.5, 'color', nred, 'DisplayName', 'Markov Policy (low state)')
plot(FigGrid, opts_high, '-', 'LineWidth', 2.5, 'color', [0.8, 0.4, 0], 'DisplayName', 'Markov Policy (high state)')
plot(FigGrid, dyncon_low, '-', 'color', nblue, 'LineWidth', 2.5, 'DisplayName', 'Dyn. consistent (low state)');
plot(FigGrid, dyncon_high, '-', 'color', [0, 0.7, 0.7], 'LineWidth', 2.5, 'DisplayName', 'Dyn. consistent (high state)');

if isfield(S, 'naive') && ~isempty(S.naive)
    plot(FigGrid, naive_low, '--', 'LineWidth', 2.5, 'color', ngreen, 'DisplayName', 'Naive solution (low state)');
    plot(FigGrid, naive_high, '--', 'LineWidth', 2.5, 'color', [0.5, 0.5, 0], 'DisplayName', 'Naive solution (high state)');
end

if isfield(S, 'showkaplan') && S.showkaplan == 1 && isfield(S, 'dyncon_kaplan') && isfield(S, 'grid_kaplan')
    plot(S.grid_kaplan, S.dyncon_kaplan(1:length(S.grid_kaplan)), 'y--', 'LineWidth', 1.5, 'DisplayName', 'Dyn. cons. (Kaplan)')
end
xlim([S.a S.b]);
legend('location', 'northwest')
legend boxoff
hold off



 % === FIGURE 6: Asset Decile Shares Comparison ===
% Assumes Decile_QH, Decile_A (in %), Gini_QH, Gini_A, nred, nblue are in workspace

figure(7); clf;
deciles = 1:10;

% bar width so pairs just touch
w      = 0.3;
pos_QH = deciles - w/2;
pos_A  = deciles + w/2;

% Quasi-hyp. bars (no edge)
bar(pos_QH, Decile_QH, w, ...
    'FaceColor', nred,  ...
    'EdgeColor', 'none');
hold on

% Aiyagari bars (no edge)
bar(pos_A,  Decile_A, w, ...
    'FaceColor', nblue, ...
    'EdgeColor', 'none');

% Axes
xticks(deciles);
xticklabels(arrayfun(@(d) sprintf('D%d',d), deciles, 'UniformOutput', false));
xlabel('Population Decile', 'Interpreter','latex','FontSize',16);
ylabel('Share of Assets',   'Interpreter','latex','FontSize',16);
xlim([0.5, 10.5]);
ylim([0, max([Decile_QH, Decile_A]) * 1.1]);

% Legend with Ginis to 2 decimals, no swatch borders
lg = legend( ...
    sprintf('Quasi-hyp. (Palma Ratio: %.2f)', Palma_QH), ...
    sprintf('Aiyagari   (Palma Ratio: %.2f)', Palma_A), ...
    'Location','northwest', ...
    'Interpreter','latex', ...
    'FontSize',16 );
legend boxoff
set(findobj(lg,'Type','Patch'),'EdgeColor','none');

box off
drawnow



%% Income distribution computation
% Extract simulated labor productivities from final period
l_actual = S.lsim(:, S.Tsim);
lK_actual = S.lsimK(:, S.TsimK);

% Compute net income: interest on assets plus wage income
income_net = S.r * S.asim(:, S.Tsim) + S.W * l_actual;
incomeK_net = S.r * S.asimK(:, S.TsimK) + S.W * lK_actual;

%% Compute Gini coefficients
% Net income
quantf_net = @(p) quantile(income_net, p);
quantfK_net = @(p) quantile(incomeK_net, p);
MEAN_net = integral(quantf_net, 0, 1);
MEANK_net = integral(quantfK_net, 0, 1);
Lorenz_net = @(p) (integral(quantf_net, 0, p)) ./ MEAN_net;
LorenzK_net = @(p) (integral(quantfK_net, 0, p)) ./ MEANK_net;

j = 1;
for t = 0:0.01:1
    LP_net(j) = Lorenz_net(t);
    LPK_net(j) = LorenzK_net(t);
    j = j + 1;
end

pgrid = 0:0.01:1;
Gini_QH_net = 1 - 2 * trapz(pgrid, LP_net);
Gini_A_net = 1 - 2 * trapz(pgrid, LPK_net);

disp(['Income Gini QH: ' num2str(Gini_QH_net, '%.3f') ' | Aiyagari: ' num2str(Gini_A_net, '%.3f')]);

%% Figure 8: Income Distributions 
figure(8); clf;
bin_width = 0.1;
h1 = histogram(income_net, 'BinWidth', bin_width, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
set(h1, 'EdgeColor', nred, 'LineWidth', 3);
hold on;
h2 = histogram(incomeK_net, 'BinWidth', bin_width, 'Normalization', 'probability', 'DisplayStyle', 'stairs');
set(h2, 'EdgeColor', nblue, 'LineWidth', 3);
dummy1 = plot(nan, nan, 'Color', nred, 'LineWidth', 2.5);
dummy2 = plot(nan, nan, 'Color', nblue, 'LineWidth', 2.5);
legend([dummy1, dummy2], {['Quasi-hyp. (mean ' num2str(mean(income_net), '%.1f') ')'], ['Aiyagari (mean ' num2str(mean(incomeK_net), '%.1f') ')']}, 'Location', 'northeast', 'FontSize', 16, 'Interpreter', 'latex');
legend boxoff;
x_max = max([max(h1.BinEdges(h1.Values > 1e-4)), max(h2.BinEdges(h2.Values > 1e-4))]);
xlim([min([income_net(:); incomeK_net(:)]), x_max]);
ylim([0, max([h1.Values, h2.Values]) * 1.1]);
xlabel('Income', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Fraction of households', 'Interpreter', 'latex', 'FontSize', 16);
box off;
hold off;

%% Figure 9: Income distributions (smoothed)
figure(9); clf;
binedges = min([income_net; incomeK_net]) : bin_width : max([income_net; incomeK_net]);
bincenters = (binedges(1:end-1) + binedges(2:end)) / 2;
h1_smooth = ksdensity(income_net, bincenters, 'Function', 'pdf', 'Bandwidth', 0.8);
h2_smooth = ksdensity(incomeK_net, bincenters, 'Function', 'pdf', 'Bandwidth', 0.8);
plot(bincenters, h1_smooth, 'Color', nred, 'LineWidth', 3);
hold on;
plot(bincenters, h2_smooth, 'Color', nblue, 'LineWidth', 3);
dummy1 = plot(nan, nan, 'Color', nred, 'LineWidth', 2.5);
dummy2 = plot(nan, nan, 'Color', nblue, 'LineWidth', 2.5);
legend([dummy1, dummy2], {['Quasi-hyp. (mean ' num2str(mean(income_net), '%.1f') ')'], ['Aiyagari (mean ' num2str(mean(incomeK_net), '%.1f') ')']}, 'Location', 'northeast', 'FontSize', 16, 'Interpreter', 'latex');
legend boxoff;
idx_last1 = find(h1_smooth > 1e-4, 1, 'last');
idx_last2 = find(h2_smooth > 1e-4, 1, 'last');
idx_last = max(idx_last1, idx_last2);
x_max = bincenters(idx_last);
xlim([min(bincenters), x_max]);
ylim([0, max([h1_smooth, h2_smooth]) * 1.1]);
xlabel('Income', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Fraction of households', 'Interpreter', 'latex', 'FontSize', 16);
box off;
hold off;
drawnow;

%% Figure 10: Lorenz Curves and Gini Coefficients
figure(10); clf;
fplot(@(x) x, [0,1], ':', 'LineWidth', 1, 'Color', uint8([5 5 5]), 'HandleVisibility', 'off');
hold on;
h1_line = plot(pgrid, LP_net, '-', 'Color', nred, 'LineWidth', 2.5);
h2_line = plot(pgrid, LPK_net, '-', 'Color', nblue, 'LineWidth', 2.5);
legend([h1_line, h2_line], {['Quasi-Hyperbolic (Gini: ' num2str(Gini_QH_net, '%.2f') ')'], ['Aiyagari (Gini: ' num2str(Gini_A_net, '%.2f') ')']}, 'Interpreter', 'latex', 'FontSize', 16);
xlim([0 1]);
xlabel('Cumulative Share of Individuals', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Cumulative Share of Net Income', 'Interpreter', 'latex', 'FontSize', 14);
legend('location', 'northwest');
legend boxoff;
hold off;

if S.extratailm == 1
%% === TAIL INEQUALITY MEASURES ===

% 1. TAIL RATIOS
% Ratio of top percentiles to median
p99_p50_QH = prctile(sorted_QH, 99) / prctile(sorted_QH, 50);
p99_p50_A = prctile(sorted_A, 99) / prctile(sorted_A, 50);

% Ratio of top percentile to mean
p99_mean_QH = prctile(sorted_QH, 99) / mean(sorted_QH);
p99_mean_A = prctile(sorted_A, 99) / mean(sorted_A);

% 2. AVERAGE WEALTH IN TAIL
% Mean wealth of top 5% relative to overall mean
top5_mean_QH = mean(sorted_QH(round(0.95*Nsim):end)) / mean(sorted_QH);
top5_mean_A = mean(sorted_A(round(0.95*Nsim):end)) / mean(sorted_A);

%% === DISPLAY RESULTS ===
fprintf('\n=== TAIL INEQUALITY COMPARISON ===\n');
fprintf('P99/P50 ratio:        QH = %.2f, A = %.2f\n', p99_p50_QH, p99_p50_A);
fprintf('P99/Mean ratio:       QH = %.2f, A = %.2f\n', p99_mean_QH, p99_mean_A);
fprintf('Top5%% mean/pop mean:  QH = %.2f, A = %.2f\n', top5_mean_QH, top5_mean_A);

%% === VISUAL COMPARISON ===
% Log-log plot of tail
figure;
top10_idx = round(0.90 * Nsim):Nsim;
loglog(1:length(top10_idx), sorted_QH(top10_idx), 'b-', 'LineWidth', 2);
hold on;
loglog(1:length(top10_idx), sorted_A(top10_idx), 'r-', 'LineWidth', 2);
xlabel('Rank from top 10%');
ylabel('Wealth');
title('Tail Comparison (log-log)');
legend('QH', 'A', 'Location', 'southwest');
grid on;  
end

if S.showval==0
    return;
end

%% Value Function Plots
warning('off', 'MATLAB:fplot:NotVectorized');
% Create function handles for each state
VF_low_func = @(z) subsref(S.VF(z), struct('type', '()', 'subs', {{':',1}}));
VF_high_func = @(z) subsref(S.VF(z), struct('type', '()', 'subs', {{':',2}}));
VFDC_low_func = @(z) subsref(S.VFDC(z), struct('type', '()', 'subs', {{':',1}}));
VFDC_high_func = @(z) subsref(S.VFDC(z), struct('type', '()', 'subs', {{':',2}}));

% Figure for low state value function
figure(11); clf;
fplot(VF_low_func, [S.a, S.b], '-', 'linewidth', 2.5, 'Color', nblue);
hold on;
fplot(VFDC_low_func, [S.a, S.b], '-', 'linewidth', 2.5, 'Color', nred);
xlabel('Assets', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Value Function (Low State)', 'Interpreter', 'latex', 'FontSize', 16);
xlim([S.a, S.b]);
hold off;

% Figure for high state value function
figure(12); clf;
fplot(VF_high_func, [S.a, S.b], '-', 'linewidth', 2.5, 'Color', nblue);
hold on;
fplot(VFDC_high_func, [S.a, S.b], '-', 'linewidth', 2.5, 'Color', nred);
xlabel('Assets', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Value Function (High State)', 'Interpreter', 'latex', 'FontSize', 16);
xlim([S.a, S.b]);
hold off;