% Computes steady state in deterministic programming problems whether they are recursive
% and non-recursive. The recursive case is uniquely given by the discount
% function dis=@(j)% delta.^j .

function parent
%global delta
clc;
close all;
clear all;
%% Precision parameters
opts = optimset('fminbnd');
opts.TolX = 1.e-32;
gridc=200;
epsilon=0.1;

%% Main parameters and functional forms
delta=0.9;
beta=0.6;
dis=@(j) beta.*delta.^j; %discount function
cutoff=800;
dis=dis(1:1:cutoff); % discount sequence

figure('Name', 'Discount Function'); 
plot([0:1:20],[1 dis(1:1:20)],'-','color','r');
hold off
u=@(c) log(c); % current utility
Du=@(c) 1/c; % future utility
v=u;
Dv=Du;
A=3; caps=0.6; % productivity parameter and capital share
f=@(a) A*a^(caps);
Df=@(a) caps*A*a^(caps-1);
DCX=@(x,y) Df(x)*Du(f(x)-y);
DCY=@(x,y) -Du(f(x)-y);
DX=@(x,y) Df(x)*Dv(f(x)-y);
DY=@(x,y) -Dv(f(x)-y);


%% optional parameters, search and plot intervals
showslope=1; % =1 shows computed slope in figure
forceshowks=0; % =1 to force showing K-S regardless of utility function (by default only shown with log utility)
aspec=0; bspec=10; % plot interval [aspec,bspec]
min=aspec; max=bspec; % by default search interval equal to plot interval
ks=@(x) ((caps*beta*delta*A)/(1-caps*delta*(1-beta)))*(x)^(caps); % Analytical solution used for plot
kss=((caps*beta*delta*A)/(1-caps*delta*(1-beta)))^(1/(1-caps))

%n=@(a) (1-delta*beta*Df(a))/(delta*(1-beta)); % slope of strategy at steady state a in K-S model
%n=@(a) (1+delta*beta*DX(a,a)/DY(a,a))/(delta*(1-beta)); % slope of strategy at steady state in general quasi-hyperbolic model


tic;
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
 

figure('Name', 'Inner Conflict Diagram'); 
plot(tt,scsc,'-','color','k','DisplayName','Current Self at Steady State');
hold on;
plot(tt,scs,'-','color','r','DisplayName',['Current Self at Steady State + ' num2str(epsilon)]);
%plot(slope,sol,'x','Markersize',14,'color','blue','DisplayName','Steady State')
plot(0.731576498490508,2,'.','Markersize',18,'color','blue','HandleVisibility','off')
plot(0.65,2,'.','Markersize',18,'color','red','Displayname',['Marginal utility = ' num2str(W(2,0.65))])
xlabel('Slope of Markov Eq. Policy');
ylabel('Supported Steady State');
legend('location','best') 
legend boxoff
hold off


figure('Name', 'First figure in SS section'); 
plot(tt,scsc,'-','color','#800000','Linewidth',2,'DisplayName','(1)');
hold on;
plot(tt,scs,'-','color','#000080','Linewidth',2,'DisplayName','(2)');
hold on;
plot(slope,sol,'.','Markersize',18,'color',[0 0.5 0],'DisplayName','Steady State')
%plot(0.731576498490508,2,'.','Markersize',18,'color','blue','HandleVisibility','off')
%plot(0.65,2,'.','Markersize',18,'color','red','Displayname',['Marginal utility = ' num2str(W(2,0.65))])
xlabel('$s$','Interpreter','latex','fontsize',12)
ylabel('$x^*$','Interpreter','latex','fontsize',12)
ylim([1 3]);
xlim([0 1]);
legend('location','best') 
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

display(['Steady state is ' num2str(sol) ' with slope ' num2str(slope)])
%display(['At steady state tolerance was ' num2str(V(sol,n(sol))) ])

 
suppertcurrentattwo=double(fminbnd(@(x) abs(W(2,x)),0,1));
%=0.731553306179005 
supportaeattwo=double(fminbnd(@(x) abs(V(2,x)),0,1));
%=0.739996090809652
% so these slopes are 0.73 and 0.74 if we round to two decimal places (to
% keep things nice).
     
%W(2,suppertcurrentattwo)
display(['At the blue dot, a=2 is supported by 0.73 as a steady state. However, the RoC 1-0.73=0.27, is too'])
display(['high from the point of view of 2+epsilon (it is to the left of the red curve). Given the RoC 0.27,  '])
display(['the marginal payoff of a+epsilon is ' num2str(V(2,suppertcurrentattwo)) ' so there is an incentive to throw a party.']);
display(['Another thing we see is that a=2 could be supported "from the left" by the a+e slope ']);
display(['0.74. Since this slope is greater than the a slope (0.73), RoC supporting a+e, 0.26, is lower than the RoC supporting a, and so if a is']);
display(['subjected to the a+e RoC then a will deviate upwards. Specifically, the marginal payoff of a is at this slope/RoC: '  num2str(W(2,supportaeattwo))]);
display(['Note that this is true for any a above the symmetric steady state: Any left-supported kinked ss is above the symmetric ss.']);
display([' ']);
display(['The in-text exercise in the paper. SS with u(c)=-c to -1 is 2.0709 with slope 0.76. Given log utility the marginal utility of a+e is here: ' num2str(V(2.0709,0.76))]);
display(['The actual support slope of a+e is ' num2str(double(fminbnd(@(x) abs(V(2.0709,x)),0,1)))]);

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


