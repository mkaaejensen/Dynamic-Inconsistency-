 

function PANEL_C
clc; close all; %clear all;

%% Notes 

% The dynamically consistent case is computationally 'fragile' as explained
% below when epsilon is assigned. The computation is therefore high precision and slightly
% slower than in the other cases.

%% Comparative statics parameters

RRA=1;                          % Rate of Risk Aversion A
PaRRA=2;                        % Rate of Risk Aversion B

%% Main parameters and functional forms
delta=0.9;
beta=1; 
dis=@(j) beta.*delta.^j;        % Discount function
cutoff=600;                     % Cut-off for computations (anything above 100 leads to very high precision)
dis=dis(1:1:cutoff);            % Discount sequence
rl=2.155; rh=2.2;               % Range in figures (imposed for easy comparison). To disable set r1=-1000;
%rl=-1000;
epsilon=0.005;                  % DC case (beta=1) is fragile numerically and requires very high precision
                                % The fragility is because terms cancel out in V and W
                                % (defined below)
[comss, comslope, compstrat]=parent(PaRRA,dis,cutoff,rl,rh,epsilon); % 2 is the comparison RRA
close all
if RRA==1
    u=@(c) log(c);
    Du=@(c) 1/c;
else
    u=@(c) c^(1-RRA)/(1-RRA);
    Du=@(c) c^(-RRA);
end
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
gridc=300;


%% optional parameters, search and plot intervals
aspec=0; bspec=10; % plot interval [aspec,bspec]
min=aspec; max=bspec; % by default search interval equal to plot interval


tic;
tt=0.02:(1/gridc):0.98;
num=numel(tt);
scs=zeros(1,num);
scsc=zeros(1,num);
flag=0;
for jj=1:1:num
     fun=@(x) abs(V(x,tt(jj)));
     func=@(x) abs(W(x,tt(jj)));%double(abs(n(x)-tt(jj)));
     scs(jj)=double(fminbnd(fun,min,max));
     scsc(jj)=double(fminbnd(func,min,max));
     marg(jj)=V(scs(jj),tt(jj));
     margc(jj)=W(scsc(jj),tt(jj));
     curr=scs(jj)-scsc(jj);
     if flag==0 & jj>1 && prev*curr<=0
         sol=scs(jj);
         slope=tt(jj);
         flag=1;
     end
        prev=curr;
end

 
flat=ones(num,1)*comss;

figure('Name', 'Panel C'); 
h1=plot(tt,scsc,'-','color','#800000','Linewidth',3); 
hold on;
h2=plot(tt,scs,'-','color','#000080','Linewidth',3); 
hold on;
h3=plot(tt,compstrat,'--','color','#000080','Linewidth',3); 
hold on;
h4=plot(slope,sol,'.','Markersize',22,'color',[0 0.5 0]); 
hold on
plot(comslope,comss,'.','Markersize',22,'color',[0 0.5 0],'HandleVisibility','off');
legend([h1 h2 h3 h4],'Isoquant of (1)','Isoquant of (2), RRA=1','Isoquant of (2), RRA=2','Steady State','Fontsize',14,'location','northeast');
xlabel('Slope ($s$)','Interpreter','latex','fontsize',14)
ylabel('State ($x^*$)','Interpreter','latex','fontsize',14)
if rl>-1000
    ylim([rl rh]);
end
xlim([0 1]);
legend boxoff
hold off
 
 
 
%% Printing output to command window 
 
display(['Panel A steady state (RRA = ' num2str(PaRRA) ') is ' num2str(comss) ' with slope ' num2str(comslope)]);
display(['Panal B steady state (RRA = ' num2str(RRA) ') is ' num2str(sol) ' with slope ' num2str(slope)])
if RRA==1
    %ks=@(x) ((caps*beta*delta*A)/(1-caps*delta*(1-beta)))*(x)^(caps);   %Analytical solution (used to calculate ss and slope following next)
    kss=((caps*beta*delta*A)/(1-caps*delta*(1-beta)))^(1/(1-caps)); % analytical steady state
    %n=@(a) (1+delta*beta*DX(a,a)/DY(a,a))/(delta*(1-beta)); % slope of strategy at steady state in general quasi-hyperbolic model
    display(['For the case RRA = 1 we know an analytical solution which for comparison is ' num2str(kss)]);
end
 

    %% Nested functions

    function val = V(a,s)
        e=epsilon; % to shorten expressions
        val=DCY(a+e,a+s*e);
        for jv=1:1:cutoff
            val=val+(dis(jv))*(s^(jv-1)*DX(a+(s^jv)*e,a+(s^(jv+1)*e))+(s^jv)*DY(a+(s^jv)*e,a+(s^(jv+1)*e)));
        end
    end

    function valw = W(a,s)
    valw=DCY(a,a);
    for j=1:1:cutoff
        valw=valw+(dis(j))*(s^(j-1)*DX(a,a)+(s^j)*DY(a,a));
    end
    end

end


