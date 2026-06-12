function [dyncons, grid] = aiyagari(a, b, n, r, W, lmin, lmax, risk_aver, P, delta, power_param)
% Aiyagari model
% Endogenous Grid Points with Markov Income
% Modified to properly handle state-dependent policies with Markov transitions

% asset grids
amax = b; % upper bound
borrow_lim = a; % borrowing limit
na = n+1; % # grid points

% computation
max_iter = 500;
tol_iter = 1.0e-6;

% UTILITY FUNCTION
if risk_aver==1
    u = @(c)log(c);
else
    u = @(c)(c.^(1-risk_aver)-1)./(1-risk_aver);
end
u1 = @(c) c.^(-risk_aver);
u1inv = @(c) c.^(-1./risk_aver);

% SET UP POWER GRID
t = linspace(0,1,na);
agrid = borrow_lim + (amax-borrow_lim)*(t.^power_param)';

% Labor endowment grid
ygrid = [lmin lmax].';
ny = length(ygrid);

R = 1+r;
yscale = 1;

% Initialize consumption guess
conguess = zeros(na,ny);
for iy = 1:ny
    conguess(:,iy) = r.*agrid + W.*yscale.*ygrid(iy);
end
con = conguess;

% Policy iteration
iter = 0;
cdiff = 1000;
while iter <= max_iter && cdiff>tol_iter
    iter = iter + 1;
    sav = zeros(na,ny);
    conlast = con;
    
    % FIXED: Compute state-specific expected marginal utilities
    emuc = zeros(na, ny); % Now state-specific (ny columns)
    for iy = 1:ny
        for iy_next = 1:ny
            % No longer using ydist - just transition probabilities
            emuc(:, iy) = emuc(:, iy) + P(iy, iy_next) * u1(conlast(:, iy_next));
        end
    end
    
    % State-specific marginal utilities and consumption
    muc1 = delta.*R.*emuc;
    con1 = zeros(na, ny);
    for iy = 1:ny
        con1(:, iy) = u1inv(muc1(:, iy));
    end
    
    % Endogenous grid method - now state-specific
    ass1 = zeros(na,ny);
    for iy = 1:ny
        ass1(:,iy) = (con1(:, iy) + agrid - W.*yscale.*ygrid(iy))./R;
        
        for ia = 1:na
            if agrid(ia)<ass1(1,iy)
                sav(ia,iy) = borrow_lim;
            else
                sav(ia,iy) = lininterp1(ass1(:,iy),agrid,agrid(ia));
                
            end
        end
        
        con(:,iy) = R.*agrid + W.*yscale.*ygrid(iy) - sav(:,iy);
    end
    
    cdiff = max(max(abs(con-conlast)));
end

% Return policy functions for both states in one vector
dyncons = [sav(:,1); sav(:,2)];
grid = agrid;

end

function [yi] = lininterp1(x,y,xi)
%% Linear interpolation function
% Same as before
placeLow = find(xi<x,1)-1;
if placeLow == 0
    placeLow = 1;
end
if isempty(placeLow)
    placeLow = length(x)-1;
end
placeHigh = placeLow+1;
xLow = x(placeLow);
xHigh = x(placeHigh);
yLow = y(placeLow);
yHigh = y(placeHigh);
yi = yLow +(xi-xLow).*(yHigh-yLow)./(xHigh-xLow);
end