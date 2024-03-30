function [wgp, xgp] = gauss(ngp)
if ngp == 1
    wgp = 2;
    xgp = 0;
elseif ngp == 2
    wgp = [1; 1];
    xgp = [-1/sqrt(3); 1/sqrt(3)];
elseif ngp == 3
    wgp = [5/9; 8/9; 5/9];
    xgp = [-sqrt(3/5); 0; sqrt(3/5)];
elseif ngp == 4
    wgp = [(18-sqrt(30))/36; (18+sqrt(30))/36; (18+sqrt(30))/36; ...
        (18-sqrt(30))/36];
    xgp = [-sqrt(3/7+2/7*sqrt(6/5)); -sqrt(3/7-2/7*sqrt(6/5)); ...
        sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7+2/7*sqrt(6/5))];
else
    error('Number of Gauss Points specified not supported')
end