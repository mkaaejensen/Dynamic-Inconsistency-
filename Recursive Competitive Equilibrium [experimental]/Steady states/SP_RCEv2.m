clc;
close all;
%clear all;

diagnostic = 1;                % =1 print extra-candidate diagnostics, =0 suppress

%% Main parameters and functional forms 
delta=0.9;
beta=0.85;
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
%gridc=300;
gridc = round( min(600, max(150, 3/epsilon)) );

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
dep=0.8;                          % Depreciation rate
f=@(a) A*a.^(caps);
Df=@(a) caps*A*a.^(caps-1);

%% Optional parameters, search and plot intervals
showslope=1;                    % =1 shows computed slope in figure
forceshowks=0;                  % (unused here, kept for consistency)
aspec=0; bspec=10;              % Plot interval [aspec,bspec]
xmin=aspec; xmax=bspec;         % Search interval

%% -------- SOCIAL PLANNER: isoquants and steady state (Panel A, SP) --------
DCX=@(x,y) (Df(x)+(1-dep)).*Du(f(x)+(1-dep)*x-y);
DCY=@(x,y) -Du(f(x)+(1-dep)*x-y);
DX=@(x,y) (Df(x)+(1-dep)).*Dv(f(x)+(1-dep)*x-y);
DY=@(x,y) -Dv(f(x)+(1-dep)*x-y);

tic;
tt=0.02:(1/gridc):0.98;
num=numel(tt);
scs   = zeros(1,num);
scsc  = zeros(1,num);
marg  = zeros(1,num);
margc = zeros(1,num);

for jj=1:1:num
    s = tt(jj);
    fun  = @(x) abs(V_SP(x,s,epsilon,DCY,DX,DY,dis,cutoff));
    func = @(x) abs(W_SP(x,s,DCY,DX,DY,dis,cutoff));
    scs(jj)   = double(fminbnd(fun,xmin,xmax));
    scsc(jj)  = double(fminbnd(func,xmin,xmax));
    marg(jj)  = V_SP(scs(jj),s,epsilon,DCY,DX,DY,dis,cutoff);
    margc(jj) = W_SP(scsc(jj),s,DCY,DX,DY,dis,cutoff);
end

phiSP = scs - scsc;

% detect all SP candidates via sign changes
cand_s_SP   = [];
cand_k_SP   = [];
cand_phi_SP = [];

for jj = 2:num
    if phiSP(jj-1)*phiSP(jj) <= 0
        if abs(phiSP(jj-1)) <= abs(phiSP(jj))
            idx = jj-1;
        else
            idx = jj;
        end
        cand_s_SP(end+1)   = tt(idx);      %#ok<AGROW>
        cand_k_SP(end+1)   = scs(idx);     %#ok<AGROW>
        cand_phi_SP(end+1) = phiSP(idx);   %#ok<AGROW>
    end
end

if isempty(cand_s_SP)
    absVals = abs(phiSP);
    minVal  = min(absVals);
    j0      = find(absVals==minVal,1,'first');
    solSP   = scs(j0);
    slopeSP = tt(j0);
else
    absVals_c = abs(cand_phi_SP);
    minVal_c  = min(absVals_c);
    iBestSP   = find(absVals_c==minVal_c,1,'first');
    solSP     = cand_k_SP(iBestSP);
    slopeSP   = cand_s_SP(iBestSP);
end

%% Plot Panel A (SP)
figure('Name', 'Panel A (SP)');
h1=plot(tt,scsc,'-','color','#800000','Linewidth',3);
hold on;
h2=plot(tt,scs,'-','color','#000080','Linewidth',3);
hold on;
h3=plot(slopeSP,solSP,'.','Markersize',22,'color',[0 0.5 0]);
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
margR  = zeros(1,numR);
margcR = zeros(1,numR);

for jj=1:1:numR
    sCE  = ttR(jj);
    funR  = @(x) abs(V_R(x,sCE,epsilon,Du,f,Df,dep,dis,cutoff));
    funcR = @(x) abs(W_R(x,sCE,Du,f,Df,dep,dis,cutoff));

    scsR(jj)   = double(fminbnd(funR,xmin,xmax));
    scscR(jj)  = double(fminbnd(funcR,xmin,xmax));
    margR(jj)  = V_R(scsR(jj),sCE,epsilon,Du,f,Df,dep,dis,cutoff);
    margcR(jj) = W_R(scscR(jj),sCE,Du,f,Df,dep,dis,cutoff);
end

phiR = scsR - scscR;

% detect all RCE candidates via sign changes
cand_s_R   = [];
cand_k_R   = [];
cand_phi_R = [];

for jj = 2:numR
    if phiR(jj-1)*phiR(jj) <= 0
        if abs(phiR(jj-1)) <= abs(phiR(jj))
            idx = jj-1;
        else
            idx = jj;
        end
        cand_s_R(end+1)   = ttR(idx);     %#ok<AGROW>
        cand_k_R(end+1)   = scsR(idx);    %#ok<AGROW>
        cand_phi_R(end+1) = phiR(idx);    %#ok<AGROW>
    end
end

if isempty(cand_s_R)
    absValsR = abs(phiR);
    minValR  = min(absValsR);
    j0R      = find(absValsR==minValR,1,'first');
    solR     = scsR(j0R);
    slopeR   = ttR(j0R);
else
    absValsR_c = abs(cand_phi_R);
    minValR_c  = min(absValsR_c);
    iBestR     = find(absValsR_c==minValR_c,1,'first');
    solR       = cand_k_R(iBestR);
    slopeR     = cand_s_R(iBestR);
end

%% Print steady states  
fprintf('Planner steady state: k^* = %.10f, s = %.10f\n',solSP,slopeSP);
fprintf('RCE steady state:     k^* = %.10f, s (private)* = %.10f\n',solR,slopeR);
fprintf('* Typically of limited relevance, slope of aggregate law of motion should be computed separately in RCE');
% Report additional SP candidates (if any, and diagnostics on)
if diagnostic && numel(cand_s_SP) > 1
    fprintf('Warning: additional SP steady-state candidates detected:\n');
    for i = 1:numel(cand_s_SP)
        if abs(cand_s_SP(i) - slopeSP) > 1e-12 || abs(cand_k_SP(i) - solSP) > 1e-12
            fprintf('  SP candidate: k^* = %.10f, s = %.10f\n',cand_k_SP(i),cand_s_SP(i));
        end
    end
end

% Report additional RCE candidates (if any, and diagnostics on)
if diagnostic && numel(cand_s_R) > 1
    fprintf('Warning: additional RCE steady-state candidates detected:\n');
    for i = 1:numel(cand_s_R)
        if abs(cand_s_R(i) - slopeR) > 1e-12 || abs(cand_k_R(i) - solR) > 1e-12
            fprintf('  RCE candidate: k^* = %.10f, s = %.10f\n',cand_k_R(i),cand_s_R(i));
        end
    end
end

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
        valwR = valwR + (dis(j))*s^(j-1)*(Fx + Fy*s);
    end
end

function valR = V_R(a,s,epsilon,Du,f,Df,dep,dis,cutoff)
    % ---- CE2 with FIXED PRICES at k* = a ----
    Rloc = @(k) 1 + Df(k)-dep;
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
