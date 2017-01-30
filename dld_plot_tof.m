function dld_plot_tof(hObject,handles) %integrated to 1d plots
%this function handles ploting the 1d information
%it should be upgraded to do bimodal fits and then use 2.58, 7.144 from pethick
%to calc T/Tc
 
t_window_min=str2double(get(handles.t_window_min_h,'String'));
t_window_max=str2double(get(handles.t_window_max_h,'String'));
T_bin_width=str2double(get(handles.time_binsize,'String'));

ymin = str2double(get(handles.ymin_h,'String'))/1000;
ymax = str2double(get(handles.ymax_h,'String'))/1000;
xmin = str2double(get(handles.xmin_h,'String'))/1000;
xmax = str2double(get(handles.xmax_h,'String'))/1000;
xybins = str2double(get(handles.oned_spatial_bins_h,'String'));

totdisplays=get(handles.oned_time_checkbox,'Value')+...
    get(handles.oned_X_checkbox,'Value')+get(handles.oned_Y_checkbox,'Value');
dispcount=1;
if get(handles.FFT_checkbox,'Value')
    width=2;
else
    width=1;
end
sfigure(150);



if get(handles.oned_time_checkbox,'Value')
    subplot(totdisplays,width,dispcount)
    dispcount=dispcount+width;
    %specify the bin edges
    T_bin_num=round((t_window_max-t_window_min)/T_bin_width);
    T_bin_num=2*floor( T_bin_num/2)+1; %round to an odd number
    T_bin_width=(t_window_max-t_window_min)/T_bin_num; %change the value to the rounded one
    %set(handles.status_handle ,'String',['No Counts in 2d Window']);
    %set(handles.time_binsize,'String',num2str(T_bin_width)); %upate the
    %used odd bins time window
    T_bin=linspace(t_window_min,t_window_max,T_bin_num);
    [TOF_counts,edges]=histcounts(handles.txy_data_windowed(:,1),T_bin);
    t_centers=mean([edges(1:end-1);edges(2:end)]);
    ydat=10^-3*TOF_counts/T_bin_width;
    plot(t_centers,ydat,'k')
    xlabel('t, time(s)');
    ylabel('Count Rate(kHz)');
    set(gcf,'Color',[1 1 1]);
    
    if get(handles.oned_fit_checkbox,'Value')
        if get(handles.fit_bimod_checkbox,'Value') 
            bmparam(1,:,:)=fit_cond(hObject,handles,t_centers,ydat,1);
        else  
            fit_therm(hObject,handles,t_centers,ydat,1,0);
        end
    end

end

if get(handles.oned_X_checkbox,'Value')
    subplot(totdisplays,width,dispcount)
    dispcount=dispcount+width;
    %specify the bin edges
    X_bin=linspace(xmin,xmax,xybins);
    X_bin_width=(xmax-xmin)/xybins;
    [X1d_counts,X1d_edges]=histcounts(handles.txy_data_windowed(:,2),X_bin);
    X1d_centers=mean([X1d_edges(1:end-1);X1d_edges(2:end)]);

    plot(X1d_centers*1000,X1d_counts/X_bin_width,'k')
    xlabel('x(mm)');
    ylabel('Linear Count Density (m^{-1})');
    set(gcf,'Color',[1 1 1]);
    
    if get(handles.oned_fit_checkbox,'Value')
        if get(handles.fit_bimod_checkbox,'Value') 
            bmparam(2,:,:)=fit_cond(hObject,handles,X1d_centers*1000,X1d_counts/X_bin_width,0);
        else
            fit_therm(hObject,handles,X1d_centers*1000,X1d_counts/X_bin_width,0,0);
        end
    end
    
end

if get(handles.oned_Y_checkbox,'Value')
    subplot(totdisplays,width,dispcount)
    dispcount=dispcount+width;
    %specify the bin edges
    Y_bin=linspace(ymin,ymax,xybins);
    Y_bin_width=(ymax-ymin)/xybins;
    [Y1d_counts,Y1d_edges]=histcounts(handles.txy_data_windowed(:,3),Y_bin);
    Y1d_centers=mean([Y1d_edges(1:end-1);Y1d_edges(2:end)]);
  
    plot(Y1d_centers*1000,Y1d_counts/Y_bin_width,'k')
    xlabel('y(mm)');
    ylabel('Linear Count Density (m^{-1})');
    set(gcf,'Color',[1 1 1]);
    
    if get(handles.oned_fit_checkbox,'Value')
        if get(handles.fit_bimod_checkbox,'Value') 
            bmparam(3,:,:)=fit_cond(hObject,handles,Y1d_centers*1000,Y1d_counts/Y_bin_width,0);
        else
            fit_therm(hObject,handles,Y1d_centers*1000,Y1d_counts/Y_bin_width,0,0);
        end
    end
end

%params in order 
    %1 parabola height (unnorm)
    %2 center
    %3 TF rad
    %4 gauss peak height
    %5 gauss width
if get(handles.fit_bimod_checkbox,'Value') && get(handles.oned_Y_checkbox,'Value') && get(handles.oned_time_checkbox,'Value')
   
end 


if get(handles.FFT_checkbox,'Value')
    %use the matlab example code https://au.mathworks.com/help/matlab/ref/fft.html
    use_window=1; %turn on the hamming window function
    start_clip=2; %clip the first n points from the fft to hide DC spike
    dispcount=2;
    if get(handles.oned_time_checkbox,'Value')
        subplot(totdisplays,width,dispcount)
        dispcount=dispcount+width;
        L=length(TOF_counts);
        Fs=1/T_bin_width;
        if use_window
            fftdat = fft(hamming(L)'.*TOF_counts/T_bin_width);    
        else
            fftdat = fft(TOF_counts/T_bin_width);   
        end
        P2 = abs(fftdat/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        f=f(start_clip:end); %remove part of the DC spike
        P1=P1(start_clip:end);
        if get(handles.FFT_log_checkbox,'Value')
            semilogy(f,P1,'k')
        else
            plot(f,P1,'k')
        end
        if use_window
        title('Single-Sided Amplitude Spectrum (hamming)')
        else
        title('Single-Sided Amplitude Spectrum (unwindowed)')
        end
        set(gcf,'Color',[1 1 1]);
        xlabel('f (Hz)')
        ylabel('Flux amplitude modulation in Hz |P1(f)|')
    end
     if get(handles.oned_X_checkbox,'Value')
        subplot(totdisplays,width,dispcount)
        dispcount=dispcount+width;
        L=length(X1d_counts);
        Fs=1/T_bin_width;
        if use_window
            fftdat = fft(hamming(L)'.*X1d_counts/X_bin_width);    
        else
            fftdat = fft(X1d_counts/X_bin_width);   
        end
        P2 = abs(fftdat/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        f=f(start_clip:end); %remove part of the DC spike
        P1=P1(start_clip:end);
        if get(handles.FFT_log_checkbox,'Value')
            semilogy(f,P1,'k')
        else
            plot(f,P1,'k')
        end
        if use_window
        title('Single-Sided Amplitude Spectrum (hamming)')
        else
        title('Single-Sided Amplitude Spectrum (unwindowed)')
        end
        set(gcf,'Color',[1 1 1]);
        xlabel('f (Hz)')
        ylabel('Flux amplitude modulation in m^{-1} |P1(f)|')
     end
     if get(handles.oned_Y_checkbox,'Value')
        subplot(totdisplays,width,dispcount)
        dispcount=dispcount+width;
        L=length(Y1d_counts);
        Fs=1/T_bin_width;
        if use_window
            fftdat = fft(hamming(L)'.*Y1d_counts/Y_bin_width);    
        else
            fftdat = fft(Y1d_counts/Y_bin_width);   
        end
        P2 = abs(fftdat/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        f=f(start_clip:end); %remove part of the DC spike
        P1=P1(start_clip:end);
        if get(handles.FFT_log_checkbox,'Value')
            semilogy(f,P1,'k')
        else
            plot(f,P1,'k')
        end
        if use_window
        title('Single-Sided Amplitude Spectrum (hamming)')
        else
        title('Single-Sided Amplitude Spectrum (unwindowed)')
        end
        set(gcf,'Color',[1 1 1]);
        xlabel('f (Hz)')
        ylabel('Flux amplitude modulation in m^{-1} |P1(f)|')
    end
    
    
    
    
  
    
    
end


% fit_temp_flag = handles.temp_fit;
% 
% if handles.tof_distribution == 1;
%     
%     z = .848; %fall distance of 848mm
%     tau = .416; %fall time of 416ms
%     TF_radius_s = get(handles.TF_radius_guesss_h,'String');%2e-4; %thomas fermi radius guess in seconds.  Add to front panel
%     TF_radius = str2double(TF_radius_s);
%     fr = 2000;  %(radial trapping freq)
%     fa = 17;  %(Axial trapping freq)
%     nmax = 200;
%     T_guess_s = get(handles.T_guess_input_h,'String');%Temp guess for thermal cloud
%     T_guess = 1e-6*str2double(T_guess_s);
%     
%     t = t_window_min:bin_size:t_window_max-bin_size;
%     
%     binned_time = zeros(floor((t_window_max-t_window_min)/bin_size),1);
%     
%     three_channel_output_sorted = handles.txy_data;
%     
%     count6a = size(three_channel_output_sorted);
%     count6 = count6a(1);
%     
%     max_bin = floor(t_window_max/bin_size);
%     min_bin = ceil(t_window_min/bin_size);
%     
%     for count7=1:count6
%         bin_for_I_plot = ceil((three_channel_output_sorted(count7,1))/bin_size)-min_bin;
%         if bin_for_I_plot > 0 && bin_for_I_plot < max_bin - min_bin
%             binned_time(bin_for_I_plot) = binned_time(bin_for_I_plot) + 1;
%         end
%         
%     end
%     
%     d = binned_time';
%     hits_in_window = sum(binned_time);
%     hits_in_window_s = num2str(hits_in_window);
%     set(handles.num_hits_tof_h ,'String',hits_in_window_s);
%     
%     min_index = min(length(t),length(d));
%     t = t(1:min_index);
%     d = d(1:min_index);
%     
%     
%     if fit_temp_flag == 1;
%         plot_i_data = 0;
%         plot_f_res = 1;
%         [t,rw,y,thermID,condID,T,Tv0,TFr,thermfrac,condfrac,mu,N01,N02,Ntherm] =TC_and_below_b(z,tau,TF_radius,fr,fa,nmax,t,d,plot_i_data,plot_f_res,T_guess) ;
%         TFr_str = num2str(TFr);
%         set(handles.text48 ,'String',TFr_str);
%         T_s = num2str(T*1e6);
%         set(handles.fitted_temp_handle ,'String',T_s);
%         thermfrac_s = num2str(thermfrac);
%         set(handles.thermfrac_h ,'String',thermfrac_s);
%         condfrac_s = num2str(condfrac);
%         set(handles.condfrac_h ,'String',condfrac_s);
%     else
%         plot_i_data = 1;
%         plot_f_res = 0;
%         TC_and_below_b(z,tau,TF_radius,fr,fa,nmax,t,d,plot_i_data,plot_f_res,T_guess) ;
% 
%     end
%     
% else %otherwise just use thermal (gauss) dist
%     
%     
%     
%     [~,T,hits_in_window] = dld_tof_plotter_window_a(handles.txy_data_windowed,t_window_min,t_window_max,150,bin_size,fit_temp_flag,spatial_window);
%     %[binned_time,T] = dld_tof_plotter(handles.txy_data,t_window_min_input,t_window_max_input,figure_number,bin_size,fit_temp_flag)
%     hits_in_window_s = num2str(hits_in_window);
%     set(handles.num_hits_tof_h ,'String',hits_in_window_s);
%     
%     if fit_temp_flag == 1;
%         T_s = num2str(T*1e6);
%         set(handles.fitted_temp_handle ,'String',T_s);
% 	end
%     
% 	if get(handles.FFT_checkbox,'Value')
% 		%pass to FFT processing:	
% 		fprintf('FFT...')
% 		tof__=handles.txy_data(:, 1);
% 		tof__windowed= tof__(tof__<=t_window_max & tof__>t_window_min);
% 		n_bins=max(tof__windowed)/bin_size;	
% 		[amplitude_fft, t_fft]=hist(tof__windowed, n_bins);			
% 		procFFT(amplitude_fft, t_fft);
% 		clear amplitude_fft t_fft tof__windowed
% 	end
% end

%change the mouse cursor to an arrow
%set(handles.figure1,'Pointer','arrow');
set(handles.status_handle ,'String','Idle');
pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update

% guidata(handles);

end

function fit_params=fit_therm(hObject,handles,xdata,ydata,istime,quiet)
%idealy this would be done with a constrained fit but ticky to implement in
%matlab

%dim = [0.2 0.6 0.3 0.3];
%str = {'Straight Line Plot','from 1 to 10'};
%annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');
amp_guess=max(ydata);
a=xlim; %cant seem to do [a,b]=xlim
xmin=a(1);
xmax=a(2);
mu_guess=sum(xdata.*ydata)/sum(ydata); %compute the weighted mean
%mu_guess=(xmin+xmax)/2;
sig_guess=sum((xdata-mu_guess).^2.*ydata)/sum(ydata); %compute the mean square weighted deviation
%sig_guess=(xmax-xmin)/1; %seems better to overestimate the width by a lot
fo = statset('TolFun',10^-6,...
    'TolX',10^-10,...
    'MaxIter',10^10,...
    'UseParallel',1);
% 'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',...
fitobject=fitnlm(xdata,ydata,...
    'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))',...
   [amp_guess,mu_guess,sig_guess],...
    'CoefficientNames',{'amp','mu','sig'},'Options',fo);
xvalues=linspace(min(xdata),max(xdata),300);
fit_params=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];

if ~quiet
    hold on
    plot(xvalues,feval(fitobject,xvalues),'r')
    %legend('data','fit')
    hold off
    %other fit params
    % num2str(fit_params(1,1),3)
    % num2str(fit_params(2,1),3)
    % num2str(fit_params(4,1),3)

    %to convert this to temp we use the expression 2.41 from pethick
    %momentum width=(mkT)^{1/2}
    %T=sqrt(pwidth)/mk
    %to convert to what we have
    %pwidth=xwidth*m/tfall
    %for the time axis we must convert to spatail using velocity
    %xwidth=twidth*velocity=vwidth*1/2 g t^2
    %total expression then T=(twidth^2*velocity^2)/(k*tfall^2*tfall^2)
    %handles.falldist = .848; %fall distance of 848mm
    %handles.falltime = .416; %fall time of 416ms

    if istime
        vdet=handles.grav*handles.falltime;
        units='s';
    else
        vdet=1/1000; %to cancel the plot being in mm
        units='mm';
    end
    
    fit_params(3,1)=fit_params(3,1)*vdet;
    fit_params(3,2)=fit_params(3,2)*vdet;

    T=(abs(fit_params(3,1)))^2 *handles.masshe /(handles.boltzconst*handles.falltime^2);

    str=sprintf('GaussFit Width %0.2e±%0.1e%s \nTemp.(no interactions)%0.2ek',abs(fit_params(3,1)),fit_params(3,2),units,T);
    text(0.02,0.9,str,'Units','normalized'); 
end
end


function fit_params=fit_cond(hObject,handles,xdata,ydata,istime)
%idealy this would be done with a constrained fit but ticky to implement in
%matlab
%roughly following https://arxiv.org/pdf/0811.4034.pdf

%dim = [0.2 0.6 0.3 0.3];
%str = {'Straight Line Plot','from 1 to 10'};
%annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');



thermfitparms=fit_therm(hObject,handles,xdata,ydata,istime,1); %{'amp','mu','sig'}


amp_guess=thermfitparms(1,1);
a=xlim; %cant seem to do [a,b]=xlim
xmin=a(1);
xmax=a(2);
mu_guess=thermfitparms(2,1);
sig_guess=thermfitparms(3,1);
%sig_guess=(xmax-xmin)/1; %seems better to overestimate the width by a lot
fo = statset('TolFun',10^-8,...
    'TolX',10^-10,...
    'MaxIter',10^10,...
    'UseParallel',1);
%params in order 
    %1 parabola height (unnorm)
    %2 center
    %3 TF rad
    %4 gauss peak height
    %5 gauss width
modelfun=@bimod;
fitobject=fitnlm(xdata,ydata,modelfun,...
   [amp_guess,mu_guess,sig_guess/10,amp_guess,sig_guess],'Options',fo);


xvalues=linspace(min(xdata),max(xdata),300);
fit_params=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];
hold on
plot(xvalues,justcond(fit_params(:,1),xvalues),'r--')
plot(xvalues,justtherm(fit_params(:,1),xvalues),'r')
plot(xvalues,bimod(fit_params(:,1),xvalues),'b')
%(xvalues,feval(fitobject,xvalues),'r')
%legend('data','fit')
hold off


%to convert this to temp we use the expression 2.41 from pethick
%momentum width=(mkT)^{1/2}
%T=sqrt(pwidth)/mk
%to convert to what we have
%pwidth=xwidth*m/tfall
%for the time axis we must convert to spatail using velocity
%xwidth=twidth*velocity=vwidth*1/2 g t^2
%total expression then T=(twidth^2*velocity^2)/(k*tfall^2*tfall^2)
%handles.falldist = .848; %fall distance of 848mm
%handles.falltime = .416; %fall time of 416ms

if istime
    vdet=handles.grav*handles.falltime;
    units='s';
else
    vdet=1/1000; %to cancel the plot being in mm
    units='mm';
end

%1 parabola height (unnorm)
    %2 center
    %3 TF rad
    %4 gauss peak height
    %5 gauss width
    
fit_params(3,1)=fit_params(3,1)*vdet;
fit_params(3,2)=fit_params(3,2)*vdet;
fit_params(5,1)=fit_params(5,1)*vdet;
fit_params(5,2)=fit_params(5,2)*vdet;


 %here i basicaly integrate under the curve and find the ratio
%to give Nc/Ntot beause i have changed the scale of the plots to put things
%into nice units the absolout value of the counts is a bit odd
ampTF=fit_params(1,1);
ampgauss=fit_params(4,1);
countsTF=abs((1/(4*fit_params(3,1)))*ampTF);
countsGauss=abs((1/(sqrt(2*pi)*fit_params(5,1)))*ampgauss);
Ncondfrac=countsTF/(countsTF+countsGauss);
TonTc=(1-Ncondfrac)^(1/3);

T=(abs(fit_params(5,1)))^2 *handles.masshe /(handles.boltzconst*handles.falltime^2);

%should use 2.90 from pethick to correct for our cigar trap
Tc=T/TonTc;

omegabar=(2*pi*2*pi*2*pi*handles.trapfreqrad*handles.trapfreqrad*handles.trapfreqaxial)^(1/3);



%handles.boltzconst*Tc=0.94*handles.hbar*omegabar*N^1/3
Nest=(handles.boltzconst*Tc/(0.94*handles.hbar*omegabar))^3

str=sprintf('GaussFit Radius %0.2e±%0.1e%s \nTemp.(no interactions)%0.2ek\nTF radius %0.2e±%0.1e%s\nCondensate fraction %0.1f%%\nT/Tc %0.1f%%\n Tc %0.2ek\n Est. N %0.2e',...
    abs(fit_params(3,1)),fit_params(3,2),units,T,abs(fit_params(5,1)),abs(fit_params(5,2)),units,Ncondfrac*100,TonTc*100,Tc,Nest);
text(0.02,0.9,str,'Units','normalized','VerticalAlignment','top'); 

end

function out=bimod(b,x)
  %1 parabola height (unnorm)
    %2 center
    %3 TF rad
    %4 gauss peak height
    %5 gauss width   
%therm %'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',..
b(3)=abs(b(3)); %makes the TF radius pos
b(4)=abs(b(4)); % gauss peak pos
b(1)=abs(b(1)); %cond height pos
parabola=(1-((x(:)-b(2))./b(3)).^2).^(3/2);
zerosformax=zeros(length(parabola),1);
therm=b(4)*exp(-1*((x(:)-b(2)).^2)./(2.*b(5).^2));
out=real(b(1).*max(zerosformax,parabola)+therm);
end

function out=justcond(b,x)
  %1 parabola height (unnorm)
    %2 center
    %3 TF rad
    %4 gauss peak height
    %5 gauss width   
%therm %'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',..
b(3)=abs(b(3)); %makes the TF radius pos
b(4)=abs(b(4)); % gauss peak pos
b(1)=abs(b(1)); %cond height pos

% if b(5)<b(3)*1.5
%     b(1)=0; %this is a very hacky way to set a limit
% end

parabola=(1-((x(:)-b(2))./b(3)).^2).^(3/2);
zerosformax=zeros(length(parabola),1);
out=real(b(1).*max(zerosformax,parabola));
end

function out=justtherm(b,x)
  %1 parabola height (unnorm)
    %2 center
    %3 TF rad
    %4 gauss peak height
    %5 gauss width   
%therm %'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',..
b(3)=abs(b(3)); %makes the TF radius pos
b(4)=abs(b(4)); % gauss peak pos
b(1)=abs(b(1)); %cond height pos

out=b(4)*exp(-1*((x(:)-b(2)).^2)./(2.*b(5).^2));
end
 