clear all
format long e

% %variables which will probably go in function call
filename = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\d';
% normalise_time_flag = 0;    %if 1 subtracts hardware trigger time from all event times
% reconst_4_corners_nomcp_flag = 1;   %if 1 reconstructs events with 4 corners but no matching MCP pulse
% reconst_3_corners_flag = 1; %if 1 reconstructs events with only 3 corners (both with and without MCP pulse)
% %%%%%%%%%%%%%%%%

% filename = filename_call; %must be a string, without the .txt if the file has this extension
% normalise_time_flag = normalise_time_flag_call; %boolean variable: if 1 subtracts hardware trigger time from all event times
% reconst_4_corners_nomcp_flag = reconst_4_corners_nomcp_flag_call;   %boolean variable: if 1 reconstructs events with 4 corners but no matching MCP pulse
% reconst_3_corners_flag = reconst_3_corners_flag_call; %boolean variable: if 1 reconstructs events with only 2/3 corners (both with and without MCP pulse)

max_group_time = 3400; %maximum time between first and last event in group
dead_time = 400; %time after 1 group to wait before looking for next group
tsum = 3200;%3200;
tolerance = 200;    %tolerance in bins

bin_time = 25e-12;                                          %% DLD bin size of 25 ps
v_perp_x = 5.2632e+005;
v_perp_y = 5.2632e+005;

dld_output_raw = [];
index = 29;%3368:3374;
for c_ind = index
    
filename_input = [filename,num2str(c_ind),'.txt'];
%BRYCE CHANGE
%saves ~10%
dld_output_raw =[dld_output_raw;tdc_importer(filename_input)];
end

%dld_output_raw = dlmread(filename_input, ',');

%dld_output_raw = dld_output_raw(1:1000,:);   %only used during debugging to reduce file size.  Delete at other times.

%bryce change single step
% number_triggers_matrix = size(dld_output_raw);
% num_trigs = number_triggers_matrix(1);  %number of clicks on TDC
%no time improvement
num_trigs = size(dld_output_raw,1);


[~, NewRowNumber] = sort(dld_output_raw(:,2));

dld_output_sorted = [];

dld_output_sorted = dld_output_raw(NewRowNumber,:);

% dld_output_sorted = dld_output_raw;

t_zero = dld_output_raw(1,2); %time of triggerpulse

    dld_output_sorted(:,2) = dld_output_sorted(:,2)-t_zero;

%dld_output_sorted;

%tsum reconstructing values.  tsumy correcting factor = a_y*x^3+b_y*x^2+c_y*x +d_y 
% Then tsumy = tsum + tsumy correcting factor
a_x = 0.000897;%2.0434e-18;
b_x = 0.021122;%3.6569e-12;
c_x = -2.8726;%-3.7798e-5;
d_x = 54.091;

a_y = 0.0012839;%2.9248e-18;
b_y = 0.014443;%2.5006e-12;
c_y = -2.6612;%-3.60163-5;
d_y = -104.81;


  dld_output_sorted_temp = dld_output_sorted;
  mcp_vec = dld_output_sorted_temp(dld_output_sorted_temp(:,1)==4,2);
%     x1_vec = dld_output_sorted_temp(dld_output_sorted_temp(:,1)==0,2);
%     x1_index = find(dld_output_sorted_temp(:,1)==0);
%     [~,I_tol,IC_tol] = uniquetol(x1_vec,2e-6/bin_time,'DataScale',1);
%     lgc = true(length(x1_vec),1);
%     lgc(I_tol) = false;
%     dld_output_sorted(x1_index(lgc)) = nan;
% %     
%     x2_vec = dld_output_sorted_temp(dld_output_sorted_temp(:,1)==1,2);
%     x2_index = find(dld_output_sorted_temp(:,1)==1);
%     [~,I_tol,IC_tol] = uniquetol(x2_vec,5e-6/bin_time,'DataScale',1);
%     lgc = true(length(x2_vec),1);
%     lgc(I_tol) = false;
%     dld_output_sorted(x2_index(lgc)) = nan;
    
%     y1_vec = dld_output_sorted_temp(dld_output_sorted_temp(:,1)==2,2);
%     y1_index = find(dld_output_sorted_temp(:,1)==2);
%     [~,I_tol,IC_tol] = uniquetol(y1_vec,5e-6/bin_time,'DataScale',1);
%     lgc = true(length(y1_vec),1);
%     lgc(I_tol) = false;
%     dld_output_sorted(y1_index(lgc)) = nan;
%     
%     y2_vec = dld_output_sorted_temp(dld_output_sorted_temp(:,1)==3,2);
%     y2_index = find(dld_output_sorted_temp(:,1)==3);
%     [~,I_tol,IC_tol] = uniquetol(y2_vec,3e-6/bin_time,'DataScale',1);
%     lgc = true(length(y2_vec),1);
%     lgc(I_tol) = false;
%     dld_output_sorted(y2_index(lgc)) = nan;
    
    number_detections_matrix = size(dld_output_sorted);                 %% Will equal [5*n + 1 2] for n detections ideally, the +1 is a master trigger to throw out
    number_detections = number_detections_matrix(1);                    %% Possibly overestimates size since errors will reduce this below what it should be
    number_successes = 0;                                               %% Tally successful hits
    which_row = 0;                                                      %% Index of matrix row to write to
    T_sum = tsum;%+1000;%-2000;    %Already defined                                                   %% Is precisely the time taken for signal to travel 8cm, speed 1e6m/s, bins 25e-12
    tolerance_throw = 200;%+10000;                                              %% Tolerance in what data to throw away
    tolerance_keep = 200;%200;%+10000;                                               %% Tolerance in what data to keep
    search_no = 9*4;%36;                                                     %% Seach over 9 (=36/4) complete hits
    T_spread = 0;                                                       %% Spreads in times reconstructed
    T_sum_spread = 0;
    T_sum_x = 0;
    T_sum_y = 0;

    %%%%%%Loop Constants%%%%%
    T_sum_tol_throw = T_sum + tolerance_throw;
    T_sum_tol_keep = T_sum + tolerance_keep;
    two_tolerance_keep = 2*tolerance_keep;

    five_channel_output_4corners = NaN(floor(number_detections/4),5);

    count = int32(0);
    dummy = int8(0);
    count2 = int8(0);
    dld_pulses_used = [];
    for count = 2:number_detections                                     %% For each detection event. First hit is a trigger to ignore
        if dld_output_sorted(count,1) == 0                              %% Pick out values of x1 (channel = 0) to match everything to
            x1_val = dld_output_sorted(count,2);                        %% This is x1 to match everything else to
            dummy = 0;                                                  %% Dummy variable to keep search local
            how_many_x = 0;
            %out_of_range_x = [0,0];

            for count2 = 1:search_no                                    %% Now looks at other detections near x1
                dummy = dummy+count2*(-1)^(count2);                     %% Go back 1 row, forward 2, back 3, to keep search efficient

                if (count+dummy) <= 0 ;                                  %% Quit if we go out of range
                    continue
                end
                if (count+dummy) > number_detections
                    continue
                end
                if abs(dld_output_sorted(count+dummy,2) - x1_val) > T_sum;                                  %% Quit if we go out of range
                    %                 out_of_range_x(mod(count2,2)+1) = out_of_range_x(mod(count2,2)+1) + 1;
                    %                 out_of_range_check_x = out_of_range_x(1)*out_of_range_x(2);
                    %                 if out_of_range_check_x ~= 0
                    %                     break
                    %                 end
                    continue
                end

                if dld_output_sorted(count+dummy,1) == 1                %% Pick out values of x2 (channel = 1)
                    x2_val = dld_output_sorted(count+dummy,2);

                    if abs(x1_val - x2_val) < (T_sum_tol_throw) %% See if x2 is close enough to check
                        how_many_x = how_many_x + 1;
                        if how_many_x > 1                               %% Breaks if multihit
                            break
                        else
                            x1 = x1_val;
                            x2 = x2_val;
                            which_row = which_row + 1;
                            tx = 0.5*(x1+x2-T_sum);                     %% Stores the value if everything goes correctly
                            five_channel_output_4corners(which_row,2) = x1 - tx;
                             five_channel_output_4corners(which_row,3) = x2 - tx;

                            %x1_index = count;
                            x2_index = count + dummy;   %if not frozen now, then risk dummy moving on

                            %five_channel_output_4corners(which_row,1) = tx;    %can lead to problems if no y's to match the x. Plus is unnecessary
                        end
                    end
                end
            end                                                         %% At this point, each X is now defined

            if how_many_x == ~1                                         %% If cant find X or is not unique, go to next x1
                continue
            end

            dummy2 = 0;
            how_many_y = 0;

            for count3 = 1:search_no                                    %% For each x1 now match the y1 and y2
                dummy2 = dummy2+count3*(-1)^(count3);
                %out_of_range_y = [0,0];

                if (count+dummy2) <= 0
                    continue
                end
                if (count+dummy2) > number_detections
                    continue
                end
                if abs(dld_output_sorted(count+dummy2,2) - x1_val) > T_sum;                                  %% Quit if we go out of range
                    %                 out_of_range_y(mod(count3,2)+1) = out_of_range_y(mod(count3,2)+1) + 1;
                    %                 out_of_range_check_y = out_of_range_y(1)*out_of_range_y(2);
                    %                 if out_of_range_check_y ~= 0
                    %                     break
                    %                 end
                    continue
                end

                if dld_output_sorted(count+dummy2,1) == 2               %% Pick out values of y1 (channel = 2)
                    y1_val = dld_output_sorted(count+dummy2,2);

                    dummy3 = 0;

                    for count4 = 1:search_no
                        dummy3 = dummy3+count4*(-1)^(count4);

                        if (count+dummy3) <= 0
                            continue
                        end
                        if (count+dummy3) > number_detections
                            continue
                        end
                        if abs(dld_output_sorted(count+dummy3,2) - x1_val) > T_sum;                                  %% Quit if we go out of range
                            continue
                        end

                        if dld_output_sorted(count+dummy3,1) == 3       %% Pick out values of y2 (channel = 3)
                            y2_val = dld_output_sorted(count+dummy3,2);
                            if and(abs(y1_val-y2_val) < (T_sum_tol_keep), abs(x1+x2-y1_val-y2_val) < two_tolerance_keep)%-86.6
                                how_many_y = how_many_y + 1;
                                if how_many_y > 1
                                    break
                                else
                                    y1 = y1_val;
                                    y2 = y2_val;
                                    ty = 0.5*(y1+y2-T_sum);
                                    t = (tx + ty)/2;
                                    five_channel_output_4corners(which_row,4) = y1 - ty;
                                    five_channel_output_4corners(which_row,5) = y2 - ty;
                                    five_channel_output_4corners(which_row,1) = t;

                                    dld_output_sorted(count,1) = NaN;   %now that we've read the data, NaN the entries in dld_output_sorted we've read for later removal
                                    dld_output_sorted(x2_index,1) = NaN;
                                    dld_output_sorted(count+dummy2,1) = NaN;
                                    dld_output_sorted(count+dummy3,1) = NaN;

                                    number_successes = number_successes + 1;
                                    
                                    dld_pulses_used(which_row,:) = [count,x2_index,count+dummy2,count+dummy3];

                                    %                                T_spread = T_spread + abs(tx-ty)^2;     %this is just used for statistics on tsum.  Comment out normally
                                    %                                T_sum_x = T_sum_x +(x1+x2-2*tx);
                                    %                                T_sum_y = T_sum_y +(y1+y2-2*ty);
                                    %                                T_sum_spread = T_sum_spread + abs((x1+x2-2*t) - T_sum)^2 + abs((y1+y2-2*t) - T_sum)^2;
                                end
                            end
                        end
                    end
                end
            end                                                         %% At this point, each Y is now defined
            if how_many_y == 0
                which_row = which_row - 1;
            end
        end
    end

    good_five_channel_output_4corners =~ isnan(five_channel_output_4corners(:,1));
    five_channel_output_4corners = five_channel_output_4corners(good_five_channel_output_4corners,:);

    dld_output_sorted_rowstokeep =~ isnan(dld_output_sorted(:,1));  %clear NaN rows from dld_output_sorted
    dld_output_sorted = dld_output_sorted(dld_output_sorted_rowstokeep,:);

    howmanyrows_4corners=size(five_channel_output_4corners);               %% Number of successful counts
    three_channel_output_4corners=zeros(howmanyrows_4corners(1),3);               %% Initialise matrix
    Num_4corners_hits = howmanyrows_4corners(1);

    three_channel_output_4corners(:,1)=five_channel_output_4corners(:,1)*bin_time;
    three_channel_output_4corners(:,2)=(five_channel_output_4corners(:,2)-five_channel_output_4corners(:,3))*v_perp_x*bin_time;
    three_channel_output_4corners(:,3)=(five_channel_output_4corners(:,4)-five_channel_output_4corners(:,5))*v_perp_y*bin_time;

    % testing pulses
%     mask=(five_channel_output_4corners(:,3)-five_channel_output_4corners(:,2))<-5.0e-9/bin_time;%three_channel_output_4corners(:,3)>0e-3 & three_channel_output_4corners(:,2)>0e-3;
%     three_channel_output_4corners = three_channel_output_4corners(mask,:);
    
%     three_channel_output_4corners(:,2)=(five_channel_output_4corners(mask,3)-five_channel_output_4corners(mask,2))*v_perp_x*bin_time;
%     three_channel_output_4corners(:,3)=(five_channel_output_4corners(mask,5)-five_channel_output_4corners(mask,4))*v_perp_y*bin_time;



temp1 = dld_output_sorted_temp(dld_pulses_used(:,1),:);
temp2 = dld_output_sorted_temp(dld_pulses_used(:,2),:);
temp3 = dld_output_sorted_temp(dld_pulses_used(:,3),:);
temp4 = dld_output_sorted_temp(dld_pulses_used(:,4),:);
% temp5 = dld_output_sorted_temp(dld_pulses_used(:,5),:);

mask=(five_channel_output_4corners(:,3)-five_channel_output_4corners(:,2))<-5.0e-9/bin_time;%
%%
out = temp1(mask,2)-temp1(~mask,2).';
out = out*bin_time;
x1diff = -out(out<-500e-9&out>-10e-6);%out(out<10e-6&out>-10e-6);
figure
subplot(4,1,1)
h = histogram(x1diff.*1e6,300);
xlabel('Time (mus)')

out = temp2(mask,2)-temp2(~mask,2).';
out = out*bin_time;
x2diff = out(out<10e-6&out>-10e-6);
subplot(4,1,2)
histogram(x2diff.*1e6,300)
xlabel('Time (mus)')
out = temp3(mask,2)-temp3(~mask,2).';
out = out*bin_time;
y1diff = out(out<10e-6&out>-10e-6);
subplot(4,1,3)
histogram(y1diff.*1e6,300)
xlabel('Time (mus)')
out = temp4(mask,2)-temp4(~mask,2).';
out = out*bin_time;
y2diff = out(out<10e-6&out>-10e-6);
subplot(4,1,4)
histogram(y2diff.*1e6,300)
xlabel('Time (mus)')

%%
out = temp1(mask,2)-temp1(mask,2).';
out = out*bin_time;
x1diff = out(out<10e-6&out>1e-9);
figure
subplot(4,1,1)
histogram(x1diff.*1e6,100)
xlabel('Time (mus)')

out = temp2(mask,2)-temp2(mask,2).';
out = out*bin_time;
x2diff = out(out<10e-6&out>1e-9);
subplot(4,1,2)
histogram(x2diff.*1e6,100)
xlabel('Time (mus)')
out = temp3(mask,2)-temp3(mask,2).';
out = out*bin_time;
y1diff = out(out<10e-6&out>1e-9);
subplot(4,1,3)
histogram(y1diff.*1e6,100)
xlabel('Time (mus)')
out = temp4(mask,2)-temp4(mask,2).';
out = out*bin_time;
y2diff = out(out<10e-6&out>1e-9);
subplot(4,1,4)
histogram(y2diff.*1e6,100)
xlabel('Time (mus)')

%%
out = temp1(~mask,2)-temp1(~mask,2).';
out = out*bin_time;
x1diff = out(out<10e-6&out>1e-9);
figure
subplot(4,1,1)
histogram(x1diff,100)
xlabel('Time')
ylabel('Correlation')

out = temp2(~mask,2)-temp2(~mask,2).';
out = out*bin_time;
x2diff = out(out<10e-6&out>1e-9);
subplot(4,1,2)
histogram(x2diff,100)
out = temp3(~mask,2)-temp3(~mask,2).';
out = out*bin_time;
y1diff = out(out<10e-6&out>1e-9);
subplot(4,1,3)
histogram(y1diff,100)
out = temp4(~mask,2)-temp4(~mask,2).';
out = out*bin_time;
y2diff = out(out<10e-6&out>1e-9);
subplot(4,1,4)
histogram(y2diff,100)

%%

%TDC_importer
%a super fast import function for the TDC data
%this does not handle the wrong format gracefully
%fairly sure the format is u(one decimal),u(64)
%its a bit messy to give the output as a 64bit matrix but previouse
%code is set up to theis type, also matlab wont speed up from u8 vs u64
function  data=tdc_importer(filepath)
fileID = fopen(filepath);
data = textscan(fileID,'%f64%f64',...%f32
    'Delimiter', ',', ...
    'MultipleDelimsAsOne',0,...
    'ReturnOnError', true,...
    'NumCharactersToSkip',0,... 
    'CollectOutput', true);
    %'EndOfLine','\n')
fclose(fileID);
data=data{1};
end