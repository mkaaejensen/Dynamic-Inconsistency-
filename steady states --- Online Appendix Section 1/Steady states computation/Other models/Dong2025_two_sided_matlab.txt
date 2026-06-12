clc;
close all;
clear all;

%% Main parameters and functional forms for Dong model
delta=0.9;
beta=0.95;
dis=@(j) beta.*delta.^j;        % Discount function
cutoff=600;                     % Cut-off for computations
dis=dis(1:1:cutoff);            % Discount sequence

rl=-1000; rh=2.12;              % Not used for ylim if rl=-1000
epsilon=0.01;                   % Step length in search for zero

%% Dong (2025) quadratic payoff w(a,s)
% Parameters (can be adjusted)
alpha  = 10;
b      = 1;
lambda = 0.5;
nu     = 4;
zeta   = 6;
rho    = 0.3;
y      = 1;                     % Upper bound on period-1 consumption of addictive good
sbar   = y/(1-rho);             % Upper-bound on stock

disp('Computing steady states in Dong (2025) assuming:');
disp(['  beta = ' num2str(beta) ', delta = ' num2str(delta) ', rho = ' num2str(rho)]);
disp(['  alpha = ' num2str(alpha) ', b = ' num2str(b) ', lambda = ' num2str(lambda) ...
      ', nu = ' num2str(nu) ', zeta = ' num2str(zeta) ', y = ' num2str(y)]);
disp('  w(a,s) = -(alpha/2)*a^2 + b*a - lambda*s - (nu/2)*s^2 + zeta*a*s');
disp(' ');

%% Precision parameters
opts = optimset('fminbnd');
opts.TolX = 1.e-32;
gridc=300;

%% Derivatives of w(a,s)
% w(a,s) = -(alpha/2)a^2 + b a - lambda s - (nu/2)s^2 + zeta a s
w_a = @(a,s) -alpha.*a + b + zeta.*s;
w_s = @(a,s) -lambda - nu.*s + zeta.*a;

% F(s,z) = w(z - rho*s, s)
% State = x (= s), choice = y (= z)
DCX=@(x,y) -rho.*w_a(y - rho.*x, x) + w_s(y - rho.*x, x);   % F_x
DCY=@(x,y) w_a(y - rho.*x, x);                              % F_y

DX = DCX;
DY = DCY;

%% Optional parameters, search and plot intervals
aspec=0; bspec=10;              % Plot interval [aspec,bspec]
amin=aspec; amax=bspec;         % Search interval equal to plot interval

%% Compute isoquants and steady state (right-hand)

tt=0.02:(1/gridc):0.98;
num=numel(tt);
scs=zeros(1,num);     % isoquant of (2), RH
scsc=zeros(1,num);    % isoquant of (1), RH
flag=0;
prev=0;

for jj=1:1:num
    fun=@(x) abs(V(x,tt(jj),epsilon,DCY,DX,DY,dis,cutoff));
    func=@(x) abs(W(x,tt(jj),DCY,DX,DY,dis,cutoff));
    scs(jj)=double(fminbnd(fun,amin,amax));
    scsc(jj)=double(fminbnd(func,amin,amax));
    marg(jj)=V(scs(jj),tt(jj),epsilon,DCY,DX,DY,dis,cutoff);
    margc(jj)=W(scsc(jj),tt(jj),DCY,DX,DY,dis,cutoff);
    curr=scs(jj)-scsc(jj);
    if flag==0 && jj>1 && prev*curr<=0
        sol=scs(jj);
        slope=tt(jj);
        flag=1;
        solind=jj;
    end
    prev=curr;
end

disp('--- Right-hand stable steady state ---');
disp(['  Steady-state stock s* = ' num2str(sol)]);
disp(['  Local slope of policy = ' num2str(slope)]);

%% Compute isoquants and steady state (left-hand)

scsL=zeros(1,num);    % isoquant of (2), LH
scscL=zeros(1,num);   % isoquant of (1), LH (same expression as RH)
flagL=0;
prevL=0;

for jj=1:1:num
    funL=@(x) abs(V_left(x,tt(jj),epsilon,DCY,DX,DY,dis,cutoff));
    funcL=@(x) abs(W(x,tt(jj),DCY,DX,DY,dis,cutoff));   % W expression is the same
    scsL(jj)=double(fminbnd(funL,amin,amax));
    scscL(jj)=double(fminbnd(funcL,amin,amax));
    currL=scsL(jj)-scscL(jj);
    if flagL==0 && jj>1 && prevL*currL<=0
        solL=scsL(jj);
        slopeL=tt(jj);
        flagL=1;
        solindL=jj;
    end
    prevL=currL;
end

if flagL==1
    disp(' ');
    disp('--- Left-hand stable steady state ---');
    disp(['  Steady-state stock s* = ' num2str(solL)]);
    disp(['  Local slope of policy = ' num2str(slopeL)]);
else
    disp(' ');
    disp('--- Left-hand stable steady state ---');
    disp('  No left-hand steady-state candidate found on this grid.');
end

%% Two-sided stability check (differentiability proxy)

two_sided_ok = false;
if flagL==1
    tolx = 1e-3;
    tols = 1e-3;
    dx = abs(sol - solL);
    ds = abs(slope - slopeL);
    if dx <= tolx && ds <= tols
        two_sided_ok = true;
        disp(['  Two-sided stability test: PASS (|Δs*| = ' num2str(dx) ', |Δslope| = ' num2str(ds) ').']);
    else
        disp(['  Two-sided stability test: FAIL (|Δs*| = ' num2str(dx) ', |Δslope| = ' num2str(ds) ').']);
    end
end

%% Time-consistent and effective-discount steady states (s_delta and s_delta,beta)

gamma_fun = @(s,delta_hat) ...
    ( w_a((1-rho)*s,s) - delta_hat*(rho*w_a((1-rho)*s,s) - w_s((1-rho)*s,s)) );

% Time-consistent consumer with discount factor delta: s_delta
g_tc = @(s) abs(gamma_fun(s,delta));
s_delta = fminbnd(g_tc,0,sbar);

% Effective-discount consumer (delta-hat) for given (beta,delta): s_delta_beta
delta_hat = (beta*delta)/(1 - delta + beta*delta);
g_eff = @(s) abs(gamma_fun(s,delta_hat));
s_delta_beta = fminbnd(g_eff,0,sbar);

disp(' ');
disp('--- Benchmarks for comparison ---');
disp(['  Time-consistent steady-state stock s_delta         = ' num2str(s_delta)]);
disp(['  Effective-discount steady-state stock s_delta_beta = ' num2str(s_delta_beta)]);

if sol > s_delta_beta
    disp('  => The steady state is an excessive-consumption trap (s* > s_delta_beta).');
else
    disp('  => The steady state is not an excessive-consumption trap (s* <= s_delta_beta).');
end

%% Common y-axis limits based on isoquants (so both figures are comparable)

all_y = [scs scsc scsL scscL];
ymin_plot = min(all_y);
ymax_plot = max(all_y);
yrng = ymax_plot - ymin_plot;
if yrng == 0
    yrng = 0.1;
end
pad = 0.05*yrng;
ymin_plot = ymin_plot - pad;
ymax_plot = ymax_plot + pad;

%% FIGURE 1: Original picture (one (1), one (2), RH steady state)

figure('Name', 'Dong model');
h1 = plot(tt,scsc,'-','color','#800000','Linewidth',3);  % Isoquant of (1)
hold on;
h2 = plot(tt,scs,'-','color','#000080','Linewidth',3);  % Isoquant of (2) (RH)
hold on;
h3 = plot(slope,sol,'.','Markersize',22,'color',[0 0.5 0]); % Steady state (RH)
hold on;

legend([h1 h2 h3], ...
       'Isoquant of (1)', ...
       'Isoquant of (2)', ...
       'Steady state', ...
       'Fontsize',14,'location','best');

xlabel('Slope','Interpreter','latex','fontsize',14)
ylabel('Addictive-good stock','Interpreter','latex','fontsize',14)
ylim([ymin_plot ymax_plot]);
xlim([tt(1) tt(end)]);
legend boxoff
hold off

%% FIGURE 2: Dong model --- full stability analysis
% (1) plus (2)_RH and (2)_LH, with RH/LH inequality regions indicated
% by hatching:  RH: 45-degree lines where V <= 0; LH: -45-degree lines where V_left >= 0.

figure('Name','Dong model --- full stability analysis');
hold on;

% Hatching parameters
ns_hatch = 40;
na_hatch = 40;
s_h = linspace(tt(1),tt(end),ns_hatch);
a_h = linspace(ymin_plot,ymax_plot,na_hatch);

ds = 0.015*(tt(end) - tt(1));          % half-length in slope direction
da = 0.015*(ymax_plot - ymin_plot);    % half-length in state direction

tolV  = 1e-4;   % avoid drawing exactly on isoquants RH
tolVL = 1e-4;   % avoid drawing exactly on isoquants LH

% RH hatching (V <= 0): 45-degree lines (/)
for jj = 1:ns_hatch
    s_val = s_h(jj);
    for ii = 1:na_hatch
        a_val = a_h(ii);
        v_val = V(a_val,s_val,epsilon,DCY,DX,DY,dis,cutoff);
        if v_val <= 0 && abs(v_val) > tolV
            x1 = s_val - ds;
            x2 = s_val + ds;
            y1 = a_val - da;
            y2 = a_val + da;
            plot([x1 x2],[y1 y2],'Color',[0.7 0.7 0.7],'LineWidth',0.5);
        end
    end
end

% LH hatching (V_left >= 0): -45-degree lines (\)
for jj = 1:ns_hatch
    s_val = s_h(jj);
    for ii = 1:na_hatch
        a_val = a_h(ii);
        vl_val = V_left(a_val,s_val,epsilon,DCY,DX,DY,dis,cutoff);
        if vl_val >= 0 && abs(vl_val) > tolVL
            x1 = s_val - ds;
            x2 = s_val + ds;
            y1 = a_val + da;
            y2 = a_val - da;
            plot([x1 x2],[y1 y2],'Color',[0.5 0.5 0.5],'LineWidth',0.5);
        end
    end
end

% Isoquant of (1)
hW = plot(tt,scsc,'-','color','#000000','Linewidth',2);

% Isoquant of (2) RH
hVRH = plot(tt,scs,'-','color','#000080','Linewidth',2);

% Isoquant of (2) LH
hVLH = plot(tt,scsL,'--','color','#000080','Linewidth',1.5);

% Steady states
hSRH = plot(slope,sol,'.','Markersize',22,'color',[0 0.5 0]);
if flagL==1
    hSLH = plot(slopeL,solL,'x','Markersize',10,'color',[0 0.5 0]);
end

if flagL==1
    legend([hW hVRH hVLH hSRH hSLH], ...
           'Isoquant of (1)', ...
           'Isoquant of (2) (RH)', ...
           'Isoquant of (2) (LH)', ...
           'Steady state (RH)', ...
           'Steady state (LH)', ...
           'Fontsize',14,'location','best');
else
    legend([hW hVRH hSRH], ...
           'Isoquant of (1)', ...
           'Isoquant of (2)', ...
           'Steady state', ...
           'Fontsize',14,'location','best');
end

xlabel('Slope','Interpreter','latex','fontsize',14)
ylabel('Addictive-good stock','Interpreter','latex','fontsize',14)
ylim([ymin_plot ymax_plot]);          % same vertical scale as Figure 1
xlim([tt(1) tt(end)]);
legend boxoff
hold off

%% Local functions

function valw = W(a,s,DCY,DX,DY,dis,cutoff)
    valw=DCY(a,a);
    for j=1:1:cutoff
        valw=valw + (dis(j))*s^(j-1)*(DX(a,a)+DY(a,a)*s);
    end
end

function val = V(a,s,epsilon,DCY,DX,DY,dis,cutoff)
    e=epsilon;
    val=DCY(a+e,a+s*e);
    for jv=1:1:cutoff
        val=val + (dis(jv))*(s^(jv-1)*(DX(a+(s^jv)*e,a+(s^(jv+1)*e)) ...
                   + s*DY(a+(s^jv)*e,a+(s^(jv+1)*e))));
    end
end

function val = V_left(a,s,epsilon,DCY,DX,DY,dis,cutoff)
    e=epsilon;
    val=DCY(a-e,a-s*e);
    for jv=1:1:cutoff
        val=val + (dis(jv))*(s^(jv-1)*(DX(a-(s^jv)*e,a-(s^(jv+1)*e)) ...
                   + s*DY(a-(s^jv)*e,a-(s^(jv+1)*e))));
    end
end
