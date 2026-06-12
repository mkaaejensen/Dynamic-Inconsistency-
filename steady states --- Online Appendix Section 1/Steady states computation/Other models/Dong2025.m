clc;
close all;
clear all;

%% Main parameters and functional forms for Panel A (Dong model)
delta=0.9;
beta=0.95;
dis=@(j) beta.*delta.^j;        % Discount function
cutoff=600;                     % Cut-off for computations
dis=dis(1:1:cutoff);            % Discount sequence

rl=-1000; rh=2.12;              % Range in figures (set rl=-1000 to disable y-lim enforcement)
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

disp('Computing right-hand stable steady state(s) in Dong (2025) assuming:');
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
showslope=1;                    % =1 shows computed slope in figure (unused but kept)
forceshowks=0;                  % (unused here, kept for consistency)
aspec=0; bspec=10;              % Plot interval [aspec,bspec]
min=aspec; max=bspec;           % By default search interval equal to plot interval

%% Compute isoquants and steady state (Panel A)
tic;
tt=0.02:(1/gridc):0.98;
num=numel(tt);
scs=zeros(1,num);
scsc=zeros(1,num);
flag=0;

for jj=1:1:num
    fun=@(x) abs(V(x,tt(jj),epsilon,DCY,DX,DY,dis,cutoff));
    func=@(x) abs(W(x,tt(jj),DCY,DX,DY,dis,cutoff));
    scs(jj)=double(fminbnd(fun,min,max));
    scsc(jj)=double(fminbnd(func,min,max));
    marg(jj)=V(scs(jj),tt(jj),epsilon,DCY,DX,DY,dis,cutoff);
    margc(jj)=W(scsc(jj),tt(jj),DCY,DX,DY,dis,cutoff);
    curr=scs(jj)-scsc(jj);
    if flag==0 & jj>1 & prev*curr<=0
        sol=scs(jj);
        slope=tt(jj);
        flag=1;
        solind=jj;
    end
    prev=curr;
end

disp('--- Steady state from Panel-A computation (MPE, right-hand stable) ---');
disp(['  Steady-state addiction stock s* = ' num2str(sol)]);
disp(['  Local slope of policy at s*     = ' num2str(slope)]);

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
    disp('  => The computed MPE steady state is an excessive-consumption trap (s* > s_delta_beta).');
else
    disp('  => The computed MPE steady state is not an excessive-consumption trap (s* <= s_delta_beta).');
end

%% Plot Panel A
figure('Name', 'Panel A (Dong model)');
h1=plot(tt,scsc,'-','color','#800000','Linewidth',3);
hold on;
h2=plot(tt,scs,'-','color','#000080','Linewidth',3);
hold on;
h3=plot(slope,sol,'.','Markersize',22,'color',[0 0.5 0]);
hold on;
legend([h1 h2 h3],'Isoquant of (1)','Isoquant of (2)','Steady State', ...
       'Fontsize',14,'location','best');
xlabel('Slope','Interpreter','latex','fontsize',14)
ylabel('Addictive-good stock','Interpreter','latex','fontsize',14)
if rl>-1000
    ylim([rl rh]);
end
xlim([0 1]);
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
