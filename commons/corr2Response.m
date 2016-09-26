function z = corr2Response(x, mu_0, si_0, mu_1, si_1)
% z = corr2Response(x, mu_0, si_0, mu_1, si_1)
%
% Maps correlation responses to tampering probabilities.
% -------------------------------------------------------------------------
% This function is a part of multi-scale analysis toolkit available from:
% https://github.com/pkorus/multiscale-prnu
% The code is provided without any warranty or support for educational and 
% research purposes only. See readme.md for more details.
% -------------------------------------------------------------------------
% Written by Paweł Korus, Shenzhen University and AGH University of Science 
%   and Technology
% Version: September 2016
% Contact: pkorus [at] agh [dot] edu [dot] pl
% -------------------------------------------------------------------------

    z = 1 / (1 + exp(-log(si_1/si_0) - (x - mu_1)^2 / (2*si_1^2) + (x - mu_0)^2 / (2 * si_0^2) ));

    z(isnan(z)) = 1;
    z(z > 1) = 1;
    
end