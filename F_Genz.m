function output = F_Genz(y)
%-------------------------------------------------------------------------%
% Description: function for compute the Genz product peak
% Programer: Juan Manuel Cardenas
% Date: May 10-2019 / Last modification : May 24 - 2021
%-------------------------------------------------------------------------%

% y = y';
% [d,m] = size(y);
% ind = [1:d]';
% aux = ind + 1;
% v_1 = (y + (( (-1).^aux )./aux ) ).^2 ; 
% v_2 = (d/4) + v_1 ; 
% v_3 = (d/4)./(v_2);
% f = prod(v_3);
% f = f';

[m,d] = size(y);
ind    = [1:d];
term_1 = (y+((-1).^(ind+1))./(ind+1)).^2;
term_2 = (d/4)./((d/4)+term_1);
output = prod(term_2,2);

end