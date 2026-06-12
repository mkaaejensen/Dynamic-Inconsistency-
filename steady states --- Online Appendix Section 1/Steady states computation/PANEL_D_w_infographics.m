% Computes steady state in deterministic decision problems for any discount
% sequence

function parent
%global delta
clc;
close all;
clear all;
%% Main parameters and functional forms
alpha=2;
gamma=1;
dis=@(j) (1+alpha.*j).^(-gamma./alpha); %discount function
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
f=@(a) A*a^(caps);
Df=@(a) caps*A*a^(caps-1);
DCX=@(x,y) Df(x)*Du(f(x)-y);
DCY=@(x,y) -Du(f(x)-y);
DX=@(x,y) Df(x)*Dv(f(x)-y);
DY=@(x,y) -Dv(f(x)-y);


%% optional parameters, search and plot intervals
gridc=100;
epsilon=0.2;
rl=2.2; rh=4.5;
showslope=1; % =1 shows computed slope in figure
forceshowks=-1; % =-1 always disable, =1 force showing K-S regardless of utility function (by default only shown with log utility)
aspec=0; bspec=12; % plot interval [aspec,bspec]
min=aspec; max=bspec; % by default search interval equal to plot interval
%ks=@(x) ((caps*beta*delta*A)/(1-caps*delta*(1-beta)))*(x)^(caps); % Analytical solution used for plot
%kss=((caps*beta*delta*A)/(1-caps*delta*(1-beta)))^(1/(1-caps));

%n=@(a) (1-delta*beta*Df(a))/(delta*(1-beta)); % slope of strategy at steady state a in K-S model
n=@(a) (1+delta*beta*DX(a,a)/DY(a,a))/(delta*(1-beta)); % slope of strategy at steady state in general quasi-hyperbolic model

opts = optimset('fminbnd');
opts.TolX = 1.e-32;
 

tt=0.01:(1/gridc):1;
num=numel(tt);
scs=zeros(1,num);
scsc=zeros(1,num);
 
for jj=1:1:num;
     fun=@(x) abs(V(x,tt(jj)));
     func=@(x) abs(W(x,tt(jj)));
     scs(jj)=double(fminbnd(fun,min,max));
     scsc(jj)=double(fminbnd(func,min,max));
     marg(jj)=V(scs(jj),tt(jj));
     margc(jj)=W(scsc(jj),tt(jj));
     curr=scs(jj)-scsc(jj);
     if jj>1 & prev*curr<=0;
         sol=scs(jj);
         slope=tt(jj);         
      end;
        prev=curr;
end;

figure('Name', 'Marginal Payoff Diagnostics (must lie around 0)');
plot(tt,margc,'.','color','r','DisplayName','Current Self at Steady State');
hold on
plot(tt,marg,'.','color','b','DisplayName',['Current Self at Steady State + ' num2str(epsilon)]);
xlabel('Slope of Markov Eq. Policy');
ylabel('Marginal Payoff');
legend('location','best') 
legend boxoff
hold off


figure('Name', 'Panel D'); 
h1=plot(tt,scsc,'-','color','#800000','Linewidth',3);  
hold on;
h2=plot(tt,scs,'-','color','#000080','Linewidth',3); 
hold on;
h3=plot(slope,sol,'.','Markersize',22,'color',[0 0.5 0]);
legend([h1 h2 h3],'Isoquant of (1)','Isoquant of (2)','Steady State','Fontsize',14,'location','northwest');
xlabel('Slope ($s$)','Interpreter','latex','fontsize',14)
ylabel('State ($x^*$)','Interpreter','latex','fontsize',14)
if rl>-1000;
    ylim([rl rh]);
end;
xlim([0 1]);
legend('location','best') 
legend boxoff
hold off


 
%funs=@(x,y) W(x,y);
%syms x y
%sl=@(x) fminbnd(fun,min,max);

%x=sym('x');
%y=sum('y');
%E=[V(x,y)==0,W(x,y)==0];
%range=[aspec,bspec;0,1];
%[sol,slope] = vpasolve(E,[x,y],range);
%funks=@(x) abs(ks(x)-x);
%ini=fminbnd(funks,min,max,opts)
%sol=fminbnd(fun,min,max,opts);

 
 
%% Printing output to command window

 

display(['Steady state is ' num2str(sol) ' with slope ' num2str(slope)])
%display(['At steady state tolerance was ' num2str(V(sol,n(sol))) ])
 

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


