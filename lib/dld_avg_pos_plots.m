function dld_avg_pos_plots(handles,figure_number)

use_window=1; %use hamming window for fft ?
start_clip=2; %clip the first n points of the fft to hide DC

tmin = str2double(get(handles.t_min_2d_handle,'String'));
tmax = str2double(get(handles.t_max_2d_handle,'String'));
ymin = str2double(get(handles.ymin_h,'String'))/1000;
ymax = str2double(get(handles.ymax_h,'String'))/1000;
xmin = str2double(get(handles.xmin_h,'String'))/1000;
xmax =  str2double(get(handles.xmax_h,'String'))/1000;
t_window_min=str2double(get(handles.t_window_min_h,'String'));
t_window_max=str2double(get(handles.t_window_max_h,'String'));
T_bin_width=str2double(get(handles.time_binsize,'String'));
log_fft_bol=get(handles.avg_pos_fft_log_checkbox,'Value');

T_bin_num=round((t_window_max-t_window_min)/T_bin_width);
T_bin_num=round(T_bin_num/2)+1;
T_bin_width=(t_window_max-t_window_min)/T_bin_num; %change the value to the rounded one


if ~length(handles.txy_data_windowed)==0
    T_edges=linspace(t_window_min,t_window_max,T_bin_num);
    t_centers=mean([T_edges(1:end-1);T_edges(2:end)]);
    avgsd=zeros(1,length(T_edges)-1);
    %X pos/sd
    for n=1:(length(T_edges)-1)
        mask=handles.txy_data_windowed(:,1)>T_edges(n) & handles.txy_data_windowed(:,1)<T_edges(n+1);
        xavgsd(n,:)=[mean(handles.txy_data_windowed(mask,2));std(handles.txy_data_windowed(mask,2))]*1000;
        yavgsd(n,:)=[mean(handles.txy_data_windowed(mask,3));std(handles.txy_data_windowed(mask,3))]*1000;
    end
    subplot_counter=1;
    if get(handles.avg_pos_fft_checkbox,'Value');
        subplot_width=2;
    else
         subplot_width=1;
    end
    figure(figure_number)
    set(gcf,'Color',[1 1 1]);
    subplot(4,subplot_width,subplot_counter)
    subplot_counter=subplot_counter+subplot_width;
    plot(t_centers,xavgsd(:,1),'k')
    xlabel('t(s)')
    ylabel('X avg(mm)')
    subplot(4,subplot_width,subplot_counter)
    subplot_counter=subplot_counter+subplot_width;
    plot(t_centers,xavgsd(:,2),'k')
    xlabel('t(s)')
    ylabel('X sd(mm)')
    subplot(4,subplot_width,subplot_counter)
    subplot_counter=subplot_counter+subplot_width;
    plot(t_centers,yavgsd(:,1),'k')
    xlabel('t(s)')
    ylabel('Y avg(mm)')
    subplot(4,subplot_width,subplot_counter)
    plot(t_centers,yavgsd(:,2),'k')
    xlabel('t(s)')
    ylabel('Y sd(mm)')
    
    if get(handles.avg_pos_fft_checkbox,'Value');
        subplot_counter=2;
        subplot(4,subplot_width,subplot_counter)
        subplot_counter=subplot_counter+subplot_width;
        plot_fft_with_options(t_centers,T_bin_width,use_window,xavgsd(:,1),'Yavg',start_clip,log_fft_bol)
        
        subplot(4,subplot_width,subplot_counter)
        subplot_counter=subplot_counter+subplot_width;
        plot_fft_with_options(t_centers,T_bin_width,use_window,xavgsd(:,2),'Xsd',start_clip,log_fft_bol)
        
        subplot(4,subplot_width,subplot_counter)
        subplot_counter=subplot_counter+subplot_width;
        plot_fft_with_options(t_centers,T_bin_width,use_window,yavgsd(:,1),'Yavg',start_clip,log_fft_bol)
        
        subplot(4,subplot_width,subplot_counter)
        subplot_counter=subplot_counter+subplot_width;
        plot_fft_with_options(t_centers,T_bin_width,use_window,yavgsd(:,2),'Ysd',start_clip,log_fft_bol)
        
    end
    
    
    
    
end

end


function plot_fft_with_options(t_centers,T_bin_width,use_window,dat,axislabel,start_clip,log_fft_bol)
        L=length(t_centers);
        Fs=1/T_bin_width;
        if use_window
            fftdat = fft(hamming(L).*dat);    
        else
            fftdat = fft(dat);   
        end
        P2 = abs(fftdat/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        f=f(start_clip:end); %remove part of the DC spike
        P1=P1(start_clip:end);
        if log_fft_bol
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
        ylabel([axislabel,'Postion modulation (mm) |P1(f)|'])
end