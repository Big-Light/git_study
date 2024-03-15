function [result] = function_sigma(di,alpha_value)

sigma_0=1e-4;     %西格玛0 
    d0=1;           %d0的设置
    result = sqrt(sigma_0 * sigma_0 *power(di/d0,alpha_value));
end