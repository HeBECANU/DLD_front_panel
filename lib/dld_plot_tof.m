function dld_plot_tof(hObject,handles) 
%integrated to 1d plots
%this function handles ploting the 1d information
%it should be upgraded to do bimodal fits and then use 2.58, 7.144 from pethick
%to calc T/Tc
 

num_files = handles.files_imported;
t_window_min=str2double(get(handles.t_window_min_h,'String'));
t_window_max=str2double(get(handles.t_window_max_h,'String'));
T_bin_width=str2double(get(handles.time_binsize,'String'))*1e-3;

ymin = str2double(get(handles.ymin_h,'String'))*1e-3;
ymax = str2double(get(handles.ymax_h,'String'))*1e-3;
xmin = str2double(get(handles.xmin_h,'String'))*1e-3;
xmax = str2double(get(handles.xmax_h,'String'))*1e-3;
xybinw = str2double(get(handles.oned_spatial_binw_h,'String'))*1e-3;

totdisplays=get(handles.oned_time_checkbox,'Value')+...
    get(handles.oned_X_checkbox,'Value')+get(handles.oned_Y_checkbox,'Value');
dispcount=1;
if get(handles.FFT_checkbox,'Value')
    width=2;
else
    width=1;
end

fig = stfig('DLD Front Panel: 1d count rate pofiles');



if get(handles.oned_time_checkbox,'Value')
    subplot(totdisplays,width,dispcount)
    dispcount=dispcount+width;
    shist_out_t=smooth_hist(handles.txy_data_windowed(:,1),'lims',[t_window_min,t_window_max],'sigma',T_bin_width);
    shist_out_t.count_rate.smooth=shist_out_t.count_rate.smooth*1e-3/num_files;
    plot(shist_out_t.bin.centers,shist_out_t.count_rate.smooth,'k')
    xlabel('t, time(s)');
    ylabel('Count Rate(kHz/File)');
    set(gcf,'Color',[1 1 1]);
    
    if get(handles.oned_fit_checkbox,'Value')
        if get(handles.fit_bimod_checkbox,'Value') 
            bmparam(1,:,:)=fit_condensate(hObject,handles,shist_out_t.bin.centers,shist_out_t.count_rate.smooth,1);
        else  
            fit_thermal(hObject,handles,shist_out_t.bin.centers,shist_out_t.count_rate.smooth,1,0);
        end
    end
end

if get(handles.oned_X_checkbox,'Value')
    subplot(totdisplays,width,dispcount)
    dispcount=dispcount+width;
    shist_out_x=smooth_hist(handles.txy_data_windowed(:,2),'lims',[xmin,xmax],'sigma',xybinw);
    plot(shist_out_x.bin.centers*1e3,shist_out_x.count_rate.smooth,'k')
    xlabel('x(mm)');
    ylabel('Linear Count Density (m^{-1}/File)');
    set(gcf,'Color',[1 1 1]);
    
    if get(handles.oned_fit_checkbox,'Value')
        if get(handles.fit_bimod_checkbox,'Value') 
            bmparam(2,:,:)=fit_condensate(hObject,handles,shist_out_x.bin.centers,shist_out_x.count_rate.smooth,0);
        else
            fit_thermal(hObject,handles,shist_out_x.bin.centers,shist_out_x.count_rate.smooth,0,0);
        end
    end
    
end

if get(handles.oned_Y_checkbox,'Value')
    subplot(totdisplays,width,dispcount)
    dispcount=dispcount+width;
    
    shist_out=smooth_hist(handles.txy_data_windowed(:,3),'lims',[ymin,ymax],'sigma',xybinw);
    plot(shist_out.bin.centers*1e3,shist_out.count_rate.smooth,'k')
    xlabel('y(mm)');
    ylabel('Linear Count Density (m^{-1}/File)');
    set(gcf,'Color',[1 1 1]);
    
    if get(handles.oned_fit_checkbox,'Value')
        if get(handles.fit_bimod_checkbox,'Value') 
            bmparam(3,:,:)=fit_condensate(hObject,handles,shist_out.bin.centers,shist_out.count_rate.smooth,0);
        else
            fit_thermal(hObject,handles,shist_out.bin.centers,shist_out.count_rate.smooth,0,0);
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
%     fig.xlim
    if get(handles.oned_time_checkbox,'Value')
        subplot(totdisplays,width,dispcount)
        fft_ax = gca;
        f_window = fft_ax.XLim;
        h_window = fft_ax.YLim;
        dispcount=dispcount+width;
        fft_out = fft_tx(shist_out_t.bin.centers,shist_out_t.count_rate.raw,'padding',10,'window','hamming');
        fft_out=fft_out(:,start_clip:end);%remove part of the DC spike
        plot_time_fft=plot(fft_out(1,:),abs(fft_out(2,:)),'k');
        if get(handles.FFT_log_checkbox,'Value')
            set(gca,'yscale','log')
        end

        if get(handles.zoom_checkbox, 'Value')
            fft_ax.YLim = h_window;
            fft_ax.XLim = f_window;
        end
        if use_window
        title('Single-Sided Amplitude Spectrum (hamming)')
        else
        title('Single-Sided Amplitude Spectrum (unwindowed)')
        end
        set(gcf,'Color',[1 1 1]);
        xlabel('f (Hz)')
        ylabel('Flux amplitude modulation in KHz/File |P1(f)|')
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
        P = abs(fftdat/L);
        P(2:end-1) = 2*P(2:end-1); 
        f = (0:L-1)*Fs/L;
        % if mod(size(x,1),2)==0
        %      error('odd Len')
        % end
        if mod(L,2)~=0 %odd
            f = f(1:(L+1)/2);
            P=P(1:(L+1)/2);
        else %even
            f = f(1:L/2);  % sample fs/2
            P=P(1:L/2);
        end
        
        f=f(start_clip:end); %remove part of the DC spike
        P=P(start_clip:end);
        if get(handles.FFT_log_checkbox,'Value')
            semilogy(f,P,'k')
        else
            plot(f,P,'k')
        end
        if use_window
        title('Single-Sided Amplitude Spectrum (hamming)')
        else
        title('Single-Sided Amplitude Spectrum (unwindowed)')
        end
        set(gcf,'Color',[1 1 1]);
        xlabel('f (Hz)')
        ylabel('Flux amplitude modulation in m^{-1}/File |P1(f)|')
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
        P = abs(fftdat/L);
        P(2:end-1) = 2*P(2:end-1); 
        f = (0:L-1)*Fs/L;
        % if mod(size(x,1),2)==0
        %      error('odd Len')
        % end
        if mod(L,2)~=0 %odd
            f = f(1:(L+1)/2);
            P=P(1:(L+1)/2);
        else %even
            f = f(1:L/2);  % sample fs/2
            P=P(1:L/2);
        end
        
        f=f(start_clip:end); %remove part of the DC spike
        P=P(start_clip:end);
        
        if get(handles.FFT_log_checkbox,'Value')
            semilogy(f,P,'k')
        else
            plot(f,P,'k')
        end
        if use_window
        title('Single-Sided Amplitude Spectrum (hamming)')
        else
        title('Single-Sided Amplitude Spectrum (unwindowed)')
        end
        set(gcf,'Color',[1 1 1]);
        xlabel('f (Hz)')
        ylabel('Flux amplitude modulation in m^{-1}/File |P1(f)|')
    end
    
end


set(handles.status_handle ,'String','Idle');
pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update

% guidata(handles);

end

function fit_params=fit_thermal(hObject,handles,xdata,ydata,is_time_axis,quiet)
%idealy this would be done with a constrained fit but ticky to implement in
%matlab

global const

%dim = [0.2 0.6 0.3 0.3];
%str = {'Straight Line Plot','from 1 to 10'};
%annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');
amp_guess=max(ydata);
a=xlim; %cant seem to do [a,b]=xlim
xmin=a(1);
xmax=a(2);
mu_guess=wmean(xdata,ydata); %compute the weighted mean
%mu_guess=(xmin+xmax)/2;
sig_guess=sqrt(sum((xdata-mu_guess).^2.*ydata)/sum(ydata)); %compute the mean square weighted deviation
fo = statset('TolFun',10^-6,...
    'TolX',1e-10,...
    'MaxIter',1e4,...
    'UseParallel',1);
% 'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))+off',...
fitobject=fitnlm(xdata,ydata,...
    'y~amp*exp(-1*((x1-mu)^2)/(2*sig^2))',...
   [amp_guess,mu_guess,sig_guess],...
    'CoefficientNames',{'amp','mu','sig'},'Options',fo);

fit_params=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];

if ~quiet
    x_sample_fit=col_vec(linspace(min(xdata),max(xdata),1e3));
    if is_time_axis
        xscaling=1;
    else
        xscaling=1e3;
    end
    [ysamp_val,ysamp_ci]=predict(fitobject,x_sample_fit,'Alpha',1-erf(1/sqrt(2)));
    hold on
    plot(xscaling*x_sample_fit,ysamp_ci,'color',[1,1,1].*0.5)
    plot(xscaling*x_sample_fit,ysamp_val,'r')
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

    
    %to convert this to temp we use the expression 2.41 from pethick
    %pwidth=(mkT)^{1/2}
    %vwidth=sqrt(kb T/ m )
    %T=vwidth^2 m/(kb)
    %for the time axis we must convert to spatial width using velocity
    %xwidth=twidth*velocity
    %then
    % v_width=x_width/time_fall
    %T=vwidth^2 m/(kb)
    %T=(x_width/time_fall)^2 m/(kb)
    %T=(twidth*velocity/time_fall)^2 m/(kb)
    
    %total expression then (twidth*velocity/time_fall)^2 * m/kb
    
    %handles.falldist = .848; %fall distance of 848mm
    %handles.falltime = .416; %fall time of 416ms


    if is_time_axis
        vdet=handles.grav*handles.falltime;
        width_units='s';
    else
        vdet=1; %to cancel the plot being in mm
        width_units='mm';
    end
    
    fit_params(3,1)=fit_params(3,1)*vdet;
    fit_params(3,2)=fit_params(3,2)*vdet;

    temperature_val=(abs(fit_params(3,1))/handles.falltime)^2 *const.mhe/const.kb;
    temperature_unc=temperature_val*2*fit_params(3,2)/abs(fit_params(3,1));
    temperature_str=string_value_with_unc(1e6*temperature_val,1e6*temperature_unc,'type','b','separator',0);
    width_str=string_value_with_unc(abs(fit_params(3,1)),fit_params(3,2),'type','b','separator',0);
    cen_str=string_value_with_unc(abs(fit_params(2,1)),fit_params(2,2),'type','b','separator',0);
    str=sprintf('Gauss fit \nCen %s %s \nWidth %s %s \nTemp.(no interactions)%s uk',...
        cen_str,width_units,width_str,width_units,temperature_str);
    text(0.01,0.9,str,'Units','normalized','VerticalAlignment','top','FontSize',14); 
end
end


function fit_params=fit_condensate(hObject,handles,xdata,ydata,istime)
%idealy this would be done with a constrained fit but ticky to implement in
%matlab
%roughly following https://arxiv.org/pdf/0811.4034.pdf

%dim = [0.2 0.6 0.3 0.3];
%str = {'Straight Line Plot','from 1 to 10'};
%annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');

    if istime
        xscaling=1;
    else
        xscaling=1e3;
    end

thermfitparms=fit_thermal(hObject,handles,xdata,ydata,istime,1); %{'amp','mu','sig'}


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
plot(xvalues.*xscaling,justcond(fit_params(:,1),xvalues),'r--')
plot(xvalues.*xscaling,justtherm(fit_params(:,1),xvalues),'r')
plot(xvalues.*xscaling,bimod(fit_params(:,1),xvalues),'b')
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
    vdet=1000; %to cancel the plot being in mm
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

T=(abs(fit_params(5,1)/vdet))^2 *handles.masshe /(handles.boltzconst*handles.falltime^2);
dT=(abs(fit_params(5,2)/vdet))^2 *handles.masshe /(handles.boltzconst*handles.falltime^2);

%should use 2.90 from pethick to correct for our cigar trap
Tc=T/TonTc;

omegabar=(2*pi*2*pi*2*pi*handles.trapfreqrad*handles.trapfreqrad*handles.trapfreqaxial)^(1/3);



%handles.boltzconst*Tc=0.94*handles.hbar*omegabar*N^1/3
Nest=(handles.boltzconst*Tc/(0.94*handles.hbar*omegabar))^3;

str=sprintf('Bimodal Fit \nGaussFit Radius %s %s \nTemp.(no interactions)%s nK\nTF radius %s %s\nCondensate fraction %0.1f%%\nT/Tc %0.1f%%\nTc %0.2ek\nEst. N %0.2e',...
    string_value_with_unc(abs(fit_params(3,1)),fit_params(3,2),'type','b','separator',0),units,...
    string_value_with_unc(T*1e9,dT*1e9,'type','b','separator',0),...
    string_value_with_unc(abs(fit_params(5,1)),abs(fit_params(5,2)),'type','b','separator',0),units,...
    Ncondfrac*100,TonTc*100,Tc,Nest);
text(0.02,0.9,str,'Units','normalized','VerticalAlignment','top','FontSize',14); 

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
 