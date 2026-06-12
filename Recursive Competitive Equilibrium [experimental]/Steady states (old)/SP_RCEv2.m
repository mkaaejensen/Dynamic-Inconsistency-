clc;
close all;
clear all;

%% Main parameters and functional forms for Panel A
delta=0.9;
beta=0.95;
dis=@(j) beta.*delta.^j;        % Discount function
cutoff=600;                     % Cut-off for computations
dis=dis(1:1:cutoff);            % Discount sequence
RRA=2;                          % Panel A Rate of Risk Aversion
%rl=1.85; rh=2.12;               % Range in figures (imposed for easy comparison). To disable set rl=-1000;
rl=-1000;
epsilon=0.01;                   % Step length in search for zero

%% Precision parameters
opts = optimset('fminbnd');
opts.TolX = 1.e-32;
gridc=300;

%% Utility and technology
if RRA==1
    u=@(c) log(c);
    Du=@(c) 1./c;
else
    u=@(c) c.^(1-RRA)./(1-RRA);
    Du=@(c) c.^(-RRA);
end
v=u;
Dv=Du;
A=2; caps=0.7;                  % Productivity parameter and capital share
dep=1;                        % Depreciation rate η
f=@(a) A*a.^(caps);
Df=@(a) caps*A*a.^(caps-1);

%% Optional parameters, search and plot intervals
showslope=1;                    % =1 shows computed slope in figure
forceshowks=0;                  % (unused here, kept for consistency)
aspec=0; bspec=10;              % Plot interval [aspec,bspec]
min=aspec; max=bspec;           % By default search interval equal to plot interval

%% -------- SOCIAL PLANNER: isoquants and steady state (Panel A, SP) --------
DCX=@(x,y) (Df(x)+(1-dep)).*Du(f(x)+(1-dep)*x-y);
DCY=@(x,y) -Du(f(x)+(1-dep)*x-y);
DX=@(x,y) (Df(x)+(1-dep)).*Dv(f(x)+(1-dep)*x-y);
DY=@(x,y) -Dv(f(x)+(1-dep)*x-y);

tic;
tt=0.02:(1/gridc):0.98;
num=numel(tt);
scs=zeros(1,num);
scsc=zeros(1,num);
flag=0;

for jj=1:1:num
    fun=@(x) abs(V_SP(x,tt(jj),epsilon,DCY,DX,DY,dis,cutoff));
    func=@(x) abs(W_SP(x,tt(jj),DCY,DX,DY,dis,cutoff));
    scs(jj)=double(fminbnd(fun,min,max));
    scsc(jj)=double(fminbnd(func,min,max));
    marg(jj)=V_SP(scs(jj),tt(jj),epsilon,DCY,DX,DY,dis,cutoff);
    margc(jj)=W_SP(scsc(jj),tt(jj),DCY,DX,DY,dis,cutoff);
    curr=scs(jj)-scsc(jj);
    if flag==0 && jj>1 && prev*curr<=0
        sol=scs(jj);
        slope=tt(jj);
        flag=1;
        solind=jj;
    end
    prev=curr;
end

solSP=sol;
slopeSP=slope;

%% Plot Panel A (SP)
figure('Name', 'Panel A (SP)');
h1=plot(tt,scsc,'-','color','#800000','Linewidth',3);
hold on;
h2=plot(tt,scs,'-','color','#000080','Linewidth',3);
hold on;
h3=plot(slope,sol,'.','Markersize',22,'color',[0 0.5 0]);
hold on;
legend([h1 h2 h3],'Isoquant of (1)','Isoquant of (2)','Steady State','Fontsize',14,'location','southeast');
xlabel('Slope ($s$)','Interpreter','latex','fontsize',14)
ylabel('State ($x^*$)','Interpreter','latex','fontsize',14)
if rl>-1000
    ylim([rl rh]);
end
xlim([0 1]);
legend boxoff
hold off

%% -------- RCE: isoquants and steady state (Panel A, RCE) --------
% RCE primitives (no taxes/transfers)
w=@(k) f(k)-Df(k).*k;
R=@(k) 1+Df(k)-dep;

ttR  = 0.02:(1/gridc):1.00;
numR = numel(ttR);

scsR   = zeros(1,numR);
scscR  = zeros(1,numR);
flagR  = 0;
for jj=1:1:numR
    sCE = ttR(jj);
    funR  = @(x) abs(V_R(x,sCE,epsilon,Du,f,Df,dep,dis,cutoff));
    funcR = @(x) abs(W_R(x,sCE,Du,f,Df,dep,dis,cutoff));

    scsR(jj)   = double(fminbnd(funR,min,max));
    scscR(jj)  = double(fminbnd(funcR,min,max));
    margR(jj)  = V_R(scsR(jj),sCE,epsilon,Du,f,Df,dep,dis,cutoff);
    margcR(jj) = W_R(scscR(jj),sCE,Du,f,Df,dep,dis,cutoff);

    currR = scsR(jj) - scscR(jj);
    if flagR==0 && jj>1 && prevR*currR<=0
        solR   = scsR(jj);
        slopeR = sCE;
        flagR  = 1;
        solindR = jj;
    end
    prevR = currR;
end

%% Print steady states  
fprintf('Planner steady state: k^* = %.10f, s = %.10f\n',solSP,slopeSP);
fprintf('RCE steady state:     k^* = %.10f, s = %.10f\n',solR,slopeR);

%% Plot Full range RCE
figure('Name', 'RCE');
h1R=plot(ttR,scscR,'-','color','#800000','Linewidth',3);
hold on;
h2R=plot(ttR,scsR,'-','color','#000080','Linewidth',3);
hold on;
h3R=plot(slopeR,solR,'.','Markersize',22,'color',[0 0.5 0]);
hold on;
legend([h1R h2R h3R],'Isoquant of (1) RCE','Isoquant of (2) RCE','RCE Steady State','Fontsize',14,'location','southeast');
xlabel('Slope ($s$)','Interpreter','latex','fontsize',14)
ylabel('State ($x^*$)','Interpreter','latex','fontsize',14)
if rl>-1000
    ylim([rl rh]);
end
xlim([0 1]);
legend boxoff
hold off

%% Print steady states
%fprintf('Planner steady state: k^* = %.10f, s = %.10f\n',solSP,slopeSP);
%fprintf('RCE steady state:     k^* = %.10f, s = %.10f\n',solR,slopeR);

%% -------- Local functions: SP --------
function valw = W_SP(a,s,DCY,DX,DY,dis,cutoff)
    valw=DCY(a,a);
    for j=1:1:cutoff
        valw=valw + (dis(j))*s^(j-1)*(DX(a,a)+DY(a,a)*s);
    end
end

function val = V_SP(a,s,epsilon,DCY,DX,DY,dis,cutoff)
    e=epsilon;
    val=DCY(a+e,a+s*e);
    for jv=1:1:cutoff
        val=val + (dis(jv))*(s^(jv-1)*(DX(a+(s^jv)*e,a+(s^(jv+1)*e)) + s*DY(a+(s^jv)*e,a+(s^(jv+1)*e))));
    end
end

%% -------- Local functions: RCE (Theorem 5 primitives) --------
function valwR = W_R(a,s,Du,f,Df,dep,dis,cutoff)
    % F_y(k,y) = -U'(R(k)k + w(k) - y), F_x(k,y) = R(k) U'(R(k)k + w(k) - y)
    Rloc=@(k) 1+Df(k)-dep;
    wloc=@(k) f(k)-Df(k).*k;
    c0 = Rloc(a).*a + wloc(a) - a;        % y = k at k^*
    Fy = -Du(c0);
    Fx = Rloc(a).*Du(c0);
    valwR = Fy;
    for j=1:1:cutoff
        valwR = valwR + dis(j)*s^(j-1)*(Fx + Fy*s);
    end
end

 function valR = V_R(a,s,epsilon,Du,f,Df,dep,dis,cutoff)
    % ---- CE2 with FIXED PRICES at k* = a ----
    Rloc = @(k) 1 + Df(k) - dep;
    wloc = @(k) f(k) - Df(k).*k;

    e    = epsilon;
    % Freeze prices at k* = a
    Rbar = Rloc(a);
    wbar = wloc(a);

    % F_y at (k^*+ε, k^*+sε) with prices evaluated at k^*
    c0  = Rbar.*(a + e) + wbar - (a + s*e);
    Fy0 = -Du(c0);
    valR = Fy0;

    % Continuation terms, again with prices fixed at k^*
    for jv = 1:cutoff
        k_state = a + (s^jv)*e;        % individual capital path
        y_next  = a + (s^(jv+1))*e;    % next-period capital
        c = Rbar.*k_state + wbar - y_next;

        Fx = Rbar.*Du(c);              % ∂/∂k term at fixed prices
        Fy = -Du(c);                   % ∂/∂y term at fixed prices

        valR = valR + dis(jv)*s^(jv-1)*(Fx + Fy*s);
    end
end

