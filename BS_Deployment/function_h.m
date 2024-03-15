function h = function_h(lambda, fsd, d, alpha_value)
    % 定义符号表达式
    h_expression = @(y) sqrt(2*pi)*(power(exp(1),-(y+lambda/fsd/sqrt(2)).^2)...
        .*(1+alpha_value*lambda/2/d+alpha_value.*y*fsd/d/sqrt(2)) ...
        -power(exp(1),-y^2)*(1+alpha_value*fsd.*y/d/sqrt(2))).^2....
        /integral(@(t) exp(-t.^2), sqrt(2)*y , sqrt(2)*y + lambda/fsd);

    % 返回符号表达式
    h = h_expression;
end
