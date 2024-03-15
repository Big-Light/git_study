function [num] = function_g(t_x , t_y , cs_x , cs_y , fk , v_lightspeed , alpha_value , sigma_sk)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
num = power((4*pi*fk*sqrt((t_x - cs_x)^2 + (t_y - cs_y)^2))/v_lightspeed,-alpha_value)*...
            exp(-sigma_sk^2/(2*(10/log(10))^2));
end