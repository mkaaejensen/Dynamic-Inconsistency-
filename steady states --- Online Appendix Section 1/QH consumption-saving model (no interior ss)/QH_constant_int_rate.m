% Computes steady state in deterministic programming problems (recursive
% and non-recursive). The recursive case is uniquely given by the discount
% function dis=@(j)% delta.^j .

function parent
%global delta
clc;
close all;
clear ;
%% Main parameters and functional forms
delta=0.95;
beta=0.9;
dis=@(j) beta.*delta.^j; %discount function
cutoff=800;
dis=dis(1:1:cutoff); % discount sequence

figure('Name', 'Discount Function'); 
plot([0:1:20],[1 dis(1:1:20)],'-','color','r');
xlabel('Time');
hold off
drawnow
u=@(c) log(c); % current utility
Du=@(c) 1/c; % future utility
v=@(c) log(c); % future utility
Dv=@(c) 1/c;   % future uility
A=8.2; caps=0.36; % productivity parameter and capital share
f=@(a) 1.04*a+2;
Df=@(a) 1.04;
DCX=@(x,y) Df(x)*Du(f(x)-y);
DCY=@(x,y) -Du(f(x)-y);
DX=@(x,y) Df(x)*Dv(f(x)-y);
DY=@(x,y) -Dv(f(x)-y);
DCY(0.1,0)

%% optional parameters, search and plot intervals
showslope=1; % =1 shows computed slope in figure
forceshowks=0; % =1 to force showing K-S regardless of utility function (by default only shown with log utility)
aspec=0; bspec=12; % plot interval [aspec,bspec]
minb=aspec; maxb=bspec; % by default search interval equal to plot interval
ks=@(x) ((caps*beta*delta*A)/(1-caps*delta*(1-beta)))*(x)^(caps); % Analytical solution used for plot
%kss=((caps*beta*delta*A)/(1-caps*delta*(1-beta)))^(1/(1-caps));

%n=@(a) (1-delta*beta*Df(a))/(delta*(1-beta)); % slope of strategy at steady state a in K-S model
%n=@(a) (1+delta*beta*DX(a,a)/DY(a,a))/(delta*(1-beta)); % slope of strategy at steady state in general quasi-hyperbolic model

opts = optimset('fminbnd');
opts.TolX = 1.e-32;
gridc=200;
tic;
epsilon=0.1;
tt=0.1:(1/gridc):1;
num=numel(tt);
scs=zeros(1,num);
scsc=zeros(1,num);
marg=zeros(1,num);
margc=zeros(1,num);

for jj=1:1:num
     fun=@(x) abs(V(x,tt(jj)));
     func=@(x) abs(W(x,tt(jj)));
     scs(jj)=double(fminbnd(fun,minb,maxb));
     scsc(jj)=double(fminbnd(func,minb,maxb));
     marg(jj)=V(scs(jj),tt(jj));
     margc(jj)=W(scsc(jj),tt(jj));
     curr=scs(jj)-scsc(jj);
     if jj>1 && prev*curr<=0
         sol=scs(jj);
         slope=tt(jj);         
      end
        prev=curr;
end
toc

figure('Name', 'Marginal Payoff Diagnostics (must lie around 0)');
plot(tt,margc,'.','color','r','DisplayName','Current Self at Steady State');
hold on
plot(tt,marg,'.','color','b','DisplayName',['Current Self at Steady State + ' num2str(epsilon)]);
xlabel('Slope of Markov Eq. Policy');
ylabel('Marginal Payoff');
legend('location','best') 
legend boxoff
hold off

figure('Name', 'Inner Conflict Diagram'); 
plot(tt,scsc,'-','color','k','DisplayName','Current Self at Steady State');
hold on;
plot(tt,scs,'-','color','r','DisplayName',['Current Self at Steady State + ' num2str(epsilon)]);
hold on;
plot(slope,sol,'o','Markersize',12,'color','blue','DisplayName','Steady State')
xlabel('Slope of Markov Eq. Policy');
ylabel('Supported Steady State');
legend('location','best') 
legend boxoff
hold off


figure('Name', 'Inner Conflict Diagram (Zoomed)'); 
plot(tt,scsc,'-','color','k','DisplayName','Current Self at Steady State');
hold on;
plot(tt,scs,'-','color','r','DisplayName',['Current Self at Steady State + ' num2str(epsilon)]);
xlabel('Slope of Markov Eq. Policy');
ylabel('Supported Steady State');
legend('location','best') 
xlim([(slope-(3/gridc)) (slope+(3/gridc))])
legend boxoff
hold off
 
 
 
 
 
%% Printing output to command window


if slope<0;
    display(['Slope is negative (inadmissible steady state) --- adjust search interval or']);
    display(['perform exhaustive search.']);
end;
if slope>1;
    display(['Slope is greater than 1 (inadmissible steady state) --- adjust search interval']);
    display(['perform exhaustive search.']);
end; 
display(['LHS of (1) at (0.8,2.1) is: ' num2str(W(2.2,0.8))]);
display(['LHS of (1) at (0.2,1) is: ' num2str(W(1,0.2))]);
display(['LHS of (1) at (0,0) is: ' num2str(W(0,0))]);
display(['LHS of (2) at (0.8,2.1) is: ' num2str(V(2.1,0.8))]);
display(['LHS of (2) at (0.2,1) is: ' num2str(V(1,0.2))]);
display(['LHS of (2) at (0,0) is: ' num2str(V(0,0))]);
display(['Main point is that steady state conditions']);
display(['can never hold. This is also clear from the Marginal Payoff diagnostics.']);
display(['Further, we see that at any steady state condidate, the current self']);
display(['as well as the "neighboring selves" have an incentive to reduce assets.']);
display(['It follows that the origin is the only stable steady state (any other steady state is either']);
display(['semi-stable or unstable).']);


display(['Steady state is ' num2str(sol) ' with slope ' num2str(slope)])
%display(['At steady state tolerance was ' num2str(V(sol,n(sol))) ])
%% plotting figure
figure('Name', 'Quasi-Hyperbolic Planner Model');  
if sol>=bspec
    bspec=sol+1'
end;
if sol <=aspec
    aspec=sol-1;
end;
 
fplot(@(j) j,[aspec,bspec],':','linewidth',0.5,'Color', uint8([5 5 5]),'HandleVisibility','off')
hold on
n=1000;
for t=1:1:n+1;
    x(t)=aspec+(t-1)*(bspec-aspec)/n;
    z(t)=ks(x(t));% Analytical solution from K-S 2002/2008
    slopef(t)=sol+slope*(x(t)-sol);
end;
plot(sol,sol,'ok','Color','k','Markersize',8,'MarkerFaceColor', 'red','DisplayName','Steady State (Computed)') 
hold on

flag=0;
if u(1)==log(1) & u(2)==log(2) & u(500)==log(500);
    plot(x,z,'-','Color', 'blue','LineWidth',2,'DisplayName','Krusell and Smith Analytical Solution')
    flat=1;
    hold on
end;
if flag==0 & forceshowks==1;
    plot(x,z,'-','Color', 'blue','LineWidth',2,'DisplayName','Krusell and Smith Analytical Solution')
    hold on
end;   

if showslope==1;
    plot(x,slopef,'--','Color', 'red','LineWidth',1,'DisplayName','Slope')
end;
legend('location','southeast') 
legend boxoff
hold off

function val = V(a,s)
e=epsilon; % to shorten expressions
val=zeros(1,1);
val=DCY(a+e,a+s*e);
    for j=1:1:cutoff;
        val=val+(dis(j))*(s^(j-1)*DX(a+(s^j)*e,a+(s^(j+1)*e))+(s^j)*DY(a+(s^j)*e,a+(s^(j+1)*e)));
    end;
end

function val = W(a,s)
e=epsilon; % to shorten expressions
val=DCY(a,a);
    for j=1:1:cutoff;
        val=val+(dis(j))*(s^(j-1)*DX(a,a)+(s^j)*DY(a,a));
    end;
end

end


