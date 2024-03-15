% function [v] = function_v(d, lambda, fsd, alpha_value)
%     f = @(y) 2 * (exp(-(y+lambda/sqrt(2)/fsd)^2) * (1+alpha_value*lambda/2/d+alpha_value*fsd/sqrt(2)/d *y) - exp(-y.^2) * (1+alpha_value*fsd/sqrt(2)/d *y)).^2 ./ (erfc(y) - erfc(y + lambda / sqrt(2) / fsd));
%     v = integral(f, -Inf, Inf);
% end
function [v] = function_v(lambda, fsd, d, alpha_value)


%     syms y;
% 
%     expr1 = exp(-(y+lambda/sqrt(2)/fsd)^2);
% 
%     expr2 = 1+alpha_value*lambda/2/d+alpha_value*fsd/sqrt(2)/d *y;
% 
%     expr3 = exp(-y*y);
% 
%     expr4 = 1+alpha_value*fsd/sqrt(2)/d *y;
% 
%     f =  2 * (expr1 * expr2 -expr3 * expr4)^2 / ...
%         (erfc(y) - erfc(y + lambda / sqrt(2) / fsd));
%     
    
    f = @(y) 2 * (exp(-(y+lambda/sqrt(2)/fsd).^2) .* ...
        (1+alpha_value*lambda/2/d+alpha_value*fsd/sqrt(2)/d *y) ...
        - exp(-y.^2) .* (1+alpha_value*fsd/sqrt(2)/d *y)).^2 ./ ...
    (erfc(y) - erfc(y + lambda / sqrt(2) / fsd));
    v = integral(f, 0 , 9.999);
end
