% Computes steady state in deterministic programming problems whether they are recursive
% and non-recursive. The recursive case is uniquely given by the discount
% function dis=@(j)% delta.^j .

function parent
clc;
%close all;
clear all;
 

%% Main parameters and functional forms
delta=0.9;
beta=0.95; 
dis=@(j) beta.*delta.^j;        % Discount function
cutoff=600;                     % Cut-off for computations (anything above 100 leads to very high precision)
dis=dis(1:1:cutoff);            % Discount sequence
RRA=1;                          % Panel B Rate of Risk Aversion
rl=1.8; rh=2.3;                 % Range in figures (imposed for easy comparison). To disable set r1=-1000;
epsilon=0.02;                    % step length in search for zero


if RRA==1;
    u=@(c) log(c);
    Du=@(c) 1/c;
else
    u=@(c) c^(1-RRA)/(1-RRA);
    Du=@(c) c^(-RRA);
end;
v=u;
Dv=Du;
A=2; caps=0.7; % productivity parameter and capital share
f=@(a) A*a^(caps);
Df=@(a) caps*A*a^(caps-1);
DCX=@(x,y) Df(x)*Du(f(x)-y);
DCY=@(x,y) -Du(f(x)-y);
DX=@(x,y) Df(x)*Dv(f(x)-y);
DY=@(x,y) -Dv(f(x)-y);


%% Precision parameters (note that these should ideally match with parent.m)
opts = optimset('fminbnd');
opts.TolX = 1.e-32;
gridc=200;


%% optional parameters, search and plot intervals
showslope=1; % =1 shows computed slope in figure
forceshowks=0; % =1 to force showing K-S regardless of utility function (by default only shown with log utility)
aspec=0; bspec=10; % plot interval [aspec,bspec]
min=aspec; max=bspec; % by default search interval equal to plot interval

%n=@(a) (1+delta*beta*DX(a,a)/DY(a,a))/(delta*(1-beta)); % slope of strategy at steady state in general quasi-hyperbolic model


tic;
tt=0.02:(1/gridc):0.98;
num=numel(tt);
scs=zeros(1,num);
scsc=zeros(1,num);
flag=0;
for jj=1:1:num;
     fun=@(x) abs(V(x,tt(jj)));
     func=@(x) abs(W(x,tt(jj)));%double(abs(n(x)-tt(jj)));
     scs(jj)=double(fminbnd(fun,min,max));
     scsc(jj)=double(fminbnd(func,min,max));
     marg(jj)=V(scs(jj),tt(jj));
     margc(jj)=W(scsc(jj),tt(jj));
     curr=scs(jj)-scsc(jj);
     if flag==0 & jj>1 & prev*curr<=0;
         sol=scs(jj);
         slope=tt(jj);
         flag=1;
     end;
        prev=curr;
end;
 
 
figure('Name', 'Panel B'); 
h1=plot(tt,scsc,'-','color','#800000','Linewidth',3); 
hold on;
h2=plot(tt,scs,'-','color','#000080','Linewidth',3); 
hold on;
h3=plot(slope,sol,'.','Markersize',22,'color',[0 0.5 0]);
legend([h1 h2 h3],'Isoquant of (1)','Isoquant of (2)','Steady State','Fontsize',14,'location','northeast');
xlabel('Slope ($s$)','Interpreter','latex','fontsize',14)
ylabel('State ($x^*$)','Interpreter','latex','fontsize',14)
if rl>-1000;
    ylim([rl rh]);
end;
xlim([0 1]);

legend boxoff
hold off
 
 
 
%% Printing output to command window and plotting Panel B
 
display(['Panel B steady state (RRA = ' num2str(RRA) ') is ' num2str(sol) ' with slope ' num2str(slope)])
display(['LHS of (1) at (0.8,2.1) is: ' num2str(W(2.1,0.8))]);
display(['LHS of (2) at (0.8,2.1) is: ' num2str(V(2.1,0.8))]);
display(['LHS of (1) at (0.1,0.01) is: ' num2str(W(0.01,0.1))]);
display(['LHS of (2) at (0.1,0.01) is: ' num2str(V(0.01,0.1))]);
display(['Approaching the origin (for any slope) leads to positive left-hand sides']);
display(['This is as we would expect given location of steady state and direction of upper contour sets.']);
   
 
%% Nested functions

function valw = W(a,s)
e=epsilon; % to shorten expressions
valw=DCY(a,a);
    for j=1:1:cutoff;
        valw=valw+(dis(j))*s^(j-1)*(DX(a,a)+DY(a,a)*s);
    end;
end


function val = V(a,s)
e=epsilon; % to shorten expressions
val=DCY(a+e,a+s*e);
    for jv=1:1:cutoff;
        val=val+(dis(jv))*(s^(jv-1)*(DX(a+(s^jv)*e,a+(s^(jv+1)*e))+s*DY(a+(s^jv)*e,a+(s^(jv+1)*e))));
    end;
end

end


