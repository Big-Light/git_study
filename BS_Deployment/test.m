v_x = [1,1,1,0,0];
M_X = v_x * v_x';
M_V = diag(v_x);
l = [1,1,1,1,1];
disp(l * v_x');

disp(sqrt(trace(M_X .* M_V)));