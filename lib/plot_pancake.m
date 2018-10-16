function [] = plot_pancake(filename_no_ext)

three_channel_output_combined_sorted = [];

for filenum = 1:10

[number_hits,three_channel_output_combined_sorted_temp] = dld_read_5channels_reconst_multi_a(filename_no_ext,1,1,0);

%three_channel_output_combined_sorted = [three_channel_output_combined_sorted , three_channel_output_combined_sorted_temp] ;
three_channel_output_combined_sorted = cat(1, three_channel_output_combined_sorted, three_channel_output_combined_sorted_temp);

end

rot_angle = 0.61;%+3.142/4;
sin_theta = sin(rot_angle);
cos_theta = cos(rot_angle);

three_channel_output_sorted_rot(:,1) = three_channel_output_combined_sorted(:,1);
three_channel_output_sorted_rot(:,2) = three_channel_output_combined_sorted(:,2)*cos_theta - three_channel_output_combined_sorted(:,3)*sin_theta;
three_channel_output_sorted_rot(:,3) = three_channel_output_combined_sorted(:,2)*sin_theta + three_channel_output_combined_sorted(:,3)*cos_theta;

t_min = 0.0;
t_max = 0.1;
x_min = -0.02;
x_max = 0.00;
y_min = -0.03;
y_max = 0.03;


rowstokeep = three_channel_output_sorted_rot(:,1)<t_max;
three_channel_output_sorted_rot = three_channel_output_sorted_rot(rowstokeep,:);

rowstokeep = three_channel_output_sorted_rot(:,1)>t_min;
three_channel_output_sorted_rot = three_channel_output_sorted_rot(rowstokeep,:);

rowstokeep = three_channel_output_sorted_rot(:,2)<x_max;
three_channel_output_sorted_rot = three_channel_output_sorted_rot(rowstokeep,:);

rowstokeep = three_channel_output_sorted_rot(:,2)>x_min;
three_channel_output_sorted_rot = three_channel_output_sorted_rot(rowstokeep,:);

rowstokeep = three_channel_output_sorted_rot(:,3)<y_max;
three_channel_output_sorted_rot = three_channel_output_sorted_rot(rowstokeep,:);

rowstokeep = three_channel_output_sorted_rot(:,3)>y_min;
three_channel_output_sorted_rot = three_channel_output_sorted_rot(rowstokeep,:);

size(three_channel_output_sorted_rot)

bin_y = 200e-6;
bin_t = 50e-6;
num_t_bins = ceil((t_max-t_min)/bin_t);
num_y_bins = ceil((y_max-y_min)/bin_y);

figure(321)
scatter3(three_channel_output_sorted_rot(:,2),three_channel_output_sorted_rot(:,3),three_channel_output_sorted_rot(:,1),'.')

% bin_y = (y_max-y_min)/100;
% bin_t = (t_max-t_min)/100;
% num_t_bins = 100;
% num_y_bins = 100;

image_x = zeros(num_y_bins,num_t_bins);



for hit_num = 1:length(three_channel_output_sorted_rot)
%     (three_channel_output_sorted_rot(hit_num,1)-t_min)
%     (three_channel_output_sorted_rot(hit_num,3)-y_min)
    which_t_bin = ceil((three_channel_output_sorted_rot(hit_num,1)-t_min)/bin_t);
    which_y_bin = ceil((three_channel_output_sorted_rot(hit_num,3)-y_min)/bin_y);
    if (which_t_bin > 0 && which_t_bin <= num_t_bins && which_y_bin > 0 && which_y_bin <= num_y_bins)
        image_x(which_y_bin,which_t_bin) = image_x(which_y_bin,which_t_bin)+1;
    end
end

figure(123)
pcolor(image_x)
colormap gray
colorbar

end