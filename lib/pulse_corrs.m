x1_vec = dld_output_sorted_temp(dld_output_sorted_temp(:,1)==0,2);
x2_vec = dld_output_sorted_temp(dld_output_sorted_temp(:,1)==1,2);
y1_vec = dld_output_sorted_temp(dld_output_sorted_temp(:,1)==2,2);
y2_vec = dld_output_sorted_temp(dld_output_sorted_temp(:,1)==3,2);

mcp_vec = dld_output_sorted_temp(dld_output_sorted_temp(:,1)==4,2);

% x1_vec = uniquetol(x1_vec,1000e-9/bin_time,'DataScale',1);
% x2_vec = uniquetol(x2_vec,1000e-9/bin_time,'DataScale',1);
% y1_vec = uniquetol(y1_vec,1000e-9/bin_time,'DataScale',1);
% y2_vec = uniquetol(y2_vec,1000e-9/bin_time,'DataScale',1);


out = y2_vec(:) - y1_vec(:)';
out = out*bin_time;
ydiff = out(-0.3e-6<out&out<0.3e-6);

out = x2_vec(:) - x1_vec(:)';
out = out*bin_time;
xdiff = out(-0.3e-6<out&out<0.3e-6);
% out = out(out<0.3e-6);

out = x1_vec(:) - x1_vec(:)';
out = out*bin_time;
x1diff = out(out<5.5e-6&out>0.1e-9);

out = x2_vec(:) - x2_vec(:)';
out = out*bin_time;
x2diff = out(out<5.5e-6&out>0.1e-9);

out = y1_vec(:) - y1_vec(:)';
out = out*bin_time;
y1diff = out(out<5.5e-6&out>0.1e-9);

out = y2_vec(:) - y2_vec(:)';
out = out*bin_time;
y2diff = out(out<5.5e-6&out>0.11e-9);

figure(400011)
clf
subplot(2,1,1)
hist(xdiff.*1e9,200)
xlabel('X2-X1 (ns)')
subplot(2,1,2)
hist(ydiff.*1e9,200)
xlabel('Y2-Y1 (ns)')
figure(3121)
clf
subplot(4,1,1)
hist(x1diff.*1e9,100)
subplot(4,1,2)
hist(x2diff.*1e9,100)
subplot(4,1,3)
hist(y1diff.*1e9,100)
subplot(4,1,4)
hist(y2diff.*1e9,100)

%%

temp1 = dld_output_sorted_temp(dld_pulses_used(:,1),:);
temp2 = dld_output_sorted_temp(dld_pulses_used(:,2),:);
temp3 = dld_output_sorted_temp(dld_pulses_used(:,3),:);
temp4 = dld_output_sorted_temp(dld_pulses_used(:,4),:);
temp5 = dld_output_sorted_temp(dld_pulses_used(:,5),:);

mask=(five_channel_output_4corners(:,3)-five_channel_output_4corners(:,2))<30.0e-9/bin_time;%-5.0e-9/bin_time;%

out = temp1(mask,2)-temp1(~mask,2).';
out = out*bin_time;
x1diff = out(out<10e-6&out>-10e-6);
figure
subplot(4,1,1)
hist(x1diff,100)

out = temp2(mask,2)-temp2(~mask,2).';
out = out*bin_time;
x2diff = out(out<10e-6&out>-10e-6);
subplot(4,1,2)
hist(x2diff,100)
out = temp3(mask,2)-temp3(~mask,2).';
out = out*bin_time;
y1diff = out(out<10e-6&out>-10e-6);
subplot(4,1,3)
hist(y1diff,100)
out = temp4(mask,2)-temp4(~mask,2).';
out = out*bin_time;
y2diff = out(out<10e-6&out>-10e-6);
subplot(4,1,4)
hist(y2diff,100)

%%
out = temp1(mask,2)-temp1(mask,2).';
out = out*bin_time;
x1diff = out(out<10e-6&out>1e-9);
figure
subplot(4,1,1)
hist(x1diff,100)

out = temp2(mask,2)-temp2(mask,2).';
out = out*bin_time;
x2diff = out(out<10e-6&out>1e-9);
subplot(4,1,2)
hist(x2diff,100)
out = temp3(mask,2)-temp3(mask,2).';
out = out*bin_time;
y1diff = out(out<10e-6&out>1e-9);
subplot(4,1,3)
hist(y1diff,100)
out = temp4(mask,2)-temp4(mask,2).';
out = out*bin_time;
y2diff = out(out<10e-6&out>1e-9);
subplot(4,1,4)
hist(y2diff,100)

%%
out = temp1(~mask,2)-temp1(~mask,2).';
out = out*bin_time;
x1diff = out(out<10e-6&out>1e-9);
figure
subplot(4,1,1)
hist(x1diff,100)
xlabel('Time')
ylabel('Correlation')

out = temp2(~mask,2)-temp2(~mask,2).';
out = out*bin_time;
x2diff = out(out<10e-6&out>1e-9);
subplot(4,1,2)
hist(x2diff,100)
out = temp3(~mask,2)-temp3(~mask,2).';
out = out*bin_time;
y1diff = out(out<10e-6&out>1e-9);
subplot(4,1,3)
hist(y1diff,100)
out = temp4(~mask,2)-temp4(~mask,2).';
out = out*bin_time;
y2diff = out(out<10e-6&out>1e-9);
subplot(4,1,4)
hist(y2diff,100)

%%
xsum = temp1(mask,2)+temp2(mask,2)-temp5(mask,2).*2;
ysum = temp3(mask,2)+temp4(mask,2)-temp5(mask,2).*2;

xsum_g = temp1(~mask,2)+temp2(~mask,2)-temp5(~mask,2).*2;
ysum_g = temp3(~mask,2)+temp4(~mask,2)-temp5(~mask,2).*2;

figure
histogram(xsum+ysum,300,'Normalization','probability')
hold on
histogram(xsum_g+ysum_g,300,'Normalization','probability')

figure
histogram(xsum,300,'Normalization','probability')
hold on
histogram(xsum_g,300,'Normalization','probability')

figure
histogram(ysum,300,'Normalization','probability')
hold on
histogram(ysum_g,300,'Normalization','probability')