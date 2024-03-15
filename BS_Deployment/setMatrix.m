function [A] = setMatrix(matrix)
%判断一行中最大的元素，如果这个元素大于0.5，则把该元素设置为1,其余元素设置为0

% 对每一行进行操作
for i = 1:size(matrix, 1)
    % 找到当前行的最大值和对应的列索引
    [max_val, max_idx] = max(matrix(i, :));
    
    % 如果最大值大于0.5，则将该行的最大值设置为1，其余元素设置为0
    if max_val > 0.5
        matrix(i, :) = zeros(1, size(matrix, 2)); % 将整行设置为0
        matrix(i, max_idx) = 1; % 将最大值所在位置设置为1
    else
        matrix(i, :) = zeros(1, size(matrix, 2)); % 将整行设置为0
    end
end
A = matrix;

end