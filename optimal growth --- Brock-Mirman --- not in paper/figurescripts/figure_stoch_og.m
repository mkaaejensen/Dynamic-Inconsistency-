function figure_stoch_og(S)
% All workspace variables are contained in the structure S.
% Every reference uses S.<variable>.

%% Define local colors
str = '#000080';
nblue = sscanf(str(2:end),'%2x%2x%2x',[1 3]) / 255;
ngreen = [0 0.5 0];
str = '#800000';
nred = sscanf(str(2:end),'%2x%2x%2x',[1 3]) / 255;

%% FIGURE 1: Policy functions
figure(1); clf;

displayNameDC = sprintf('Dynamically Consistent ($l_t=\\underline{l}$)');
displayNameDCH = sprintf('Dynamically Consistent ($l_t=\\overline{l}$)');
displayNameMP = sprintf('Quasi-hyperbolic ($l_t=\\underline{l}$)');
displayNameMPH = sprintf('Quasi-hyperbolic ($l_t=\\overline{l}$)');

x1 = linspace(S.a, S.b, S.n+1);

% Compute high productivity policies
optsH = interp1(x1, S.opts, S.offset(x1), 'linear', 'extrap');
dynconH = interp1(x1, S.dyncon, S.offset(x1), 'linear', 'extrap');

% 45 degree line
fplot(@(j) j, [S.a, S.b], ':', 'linewidth', 1, 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off')
hold on

% Plot policies
p1 = plot(x1, S.dyncon, '--', 'color', "k", 'LineWidth', 3, 'DisplayName', displayNameDC);
p2 = plot(x1, dynconH, ':', 'color', "k", 'LineWidth', 3, 'DisplayName', displayNameDCH);
p3 = plot(x1, S.opts, '-', 'LineWidth', 3, 'color', ngreen, 'DisplayName', displayNameMP);
p4 = plot(x1, optsH, '-', 'LineWidth', 3, 'color', nred, 'DisplayName', displayNameMPH);

legend([p1, p2, p3, p4], 'location', 'northwest', 'FontSize', 16, 'Interpreter', 'latex');
legend boxoff
xlabel('Assets at date $t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Assets at date $t+1$', 'Interpreter', 'latex', 'FontSize', 16);
xlim([S.a, S.bspecfig]);
hold off
drawnow

end