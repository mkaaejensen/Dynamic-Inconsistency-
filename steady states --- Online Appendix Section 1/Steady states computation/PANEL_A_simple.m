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
rl=1.85; rh=2.12;               % Range in figures (imposed for easy comparison). To disable set rl=-1000;
epsilon=0.01;                   % Step length in search for zero

%% Precision parameters
opts = optimset('fminbnd');
opts.TolX = 1.e-32;
gridc=300;

%% Utility and technology
if RRA==1
    u=@(c) log(c);
    Du=@(c) 1/c;
else
    u=@(c) c^(1-RRA)/(1-RRA);
    Du=@(c) c^(-RRA);
end
v=u;
Dv=Du;
A=2; caps=0.7;                  % Productivity parameter and capital share
f=@(a) A*a^(caps);
Df=@(a) caps*A*a^(caps-1);
DCX=@(x,y) Df(x)*Du(f(x)-y);
DCY=@(x,y) -Du(f(x)-y);
DX=@(x,y) Df(x)*Dv(f(x)-y);
DY=@(x,y) -Dv(f(x)-y);

%% Optional parameters, search and plot intervals
showslope=1;                    % =1 shows computed slope in figure
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
    if flag==0 && jj>1 && prev*curr<=0
        sol=scs(jj);
        slope=tt(jj);
        flag=1;
        solind=jj;
    end
    prev=curr;
end

%% Plot Panel A
figure('Name', 'Panel A');
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
        val=val + (dis(jv))*(s^(jv-1)*(DX(a+(s^jv)*e,a+(s^(jv+1)*e)) + s*DY(a+(s^jv)*e,a+(s^(jv+1)*e))));
    end
end
