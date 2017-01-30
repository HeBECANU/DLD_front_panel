function y = fit_2d_gaussian(Param,binned_output)

y = 0;
size_vec = size(binned_output);
max_y = size_vec(2);
max_x = size_vec(1);

for xpos = 1:max_x
    for ypos = 1:max_y
        y = y + (Param(1)*exp(-0.5*((xpos-Param(5)).^2./(Param(2)).^2))*exp(-0.5*((ypos-Param(6)).^2./(Param(3)).^2)) + Param(4) - binned_output(xpos,ypos)).^2;
    end
end
end