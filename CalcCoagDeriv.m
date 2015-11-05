function dvdt = CalcCoagDeriv(t, vcon, p2)
% 
% CalcCoagDeriv calculates the derivates in the population balance
% equations 
%
% Note, use a new, single structure to transfer all the parameters
%
% USAGE:
%
%
% HISTORY:
%  04-05-09: First cut - based heavily on GAJ's code
%
% Adrian Burd, University of Georgia, 2009

n_sections = length(vcon);

vcon_r = vcon';      % vcon passed as column vector - make a row

vcon_shift = [0 vcon_r(1:n_sections-1)];

term1 = vcon_r * p2.b25;
term1 = vcon_r .* term1;

term2 = vcon_r * p2.b1;
term2 = term2 .* vcon_shift;

term3 = p2.linear * vcon;

dvdt = (term1 + term2)' + term3;

