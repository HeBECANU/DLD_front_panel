function[xwidth_pix,ywidth_pix] = dld_2dplotter_window_a(handles)

    %%%%%%%%%% 2d_data_plotter_window %%%%%%%%%%
    %   Takes a handle object with
        %handles.txy_data_windowed
        %handles.spatial_bins
        %handles.ymin_h
        %handles.ymax_h
        %handles.xmin_h
        %handles.xmax_h
        %handles.x_fit_guess
        %handles.y_fit_guess
        %handles.radius_fit_guess
    %   Takes data three_channel_input of the form (t,x,y) with t in s, and x,y in m. 
    %   The displayed figure has number figure_number, while zoom_input and
    %   bins_per_pixel control the zoom and resolution of the image. 
    
%factor for compressing the dynamic range
dyn_range_pow=0.5; %power between 0 and 1
    
tmin = str2double(get(handles.t_min_2d_handle,'String'));
tmax = str2double(get(handles.t_max_2d_handle,'String'));
ymin = str2double(get(handles.ymin_h,'String'))/1000;
ymax = str2double(get(handles.ymax_h,'String'))/1000;
xmin = str2double(get(handles.xmin_h,'String'))/1000;
xmax =  str2double(get(handles.xmax_h,'String'))/1000;
XY_bool=get(handles.XY_checkbox,'Value');
XT_bool=get(handles.XT_checkbox,'Value');
YT_bool=get(handles.YT_checkbox,'Value');
range_comp_bool=get(handles.LogScale2d_checkbox,'Value');
spatial_blur=str2double(get(handles.spatial_blur,'String'));
fit_flag=handles.spatial_fit;
bins=str2double(get(handles.spatial_bins_h,'String'));
if spatial_blur<0
    spatial_blur=0;
    set(handles.spatial_blur,'String','0')
end
if spatial_blur>bins/2
    spatial_blur=0;
    set(handles.spatial_blur,'String','0')
end

%set up the colormap
cmap=viridis();
if range_comp_bool
    cmap=nonlinear_colormap(cmap,'power',[dyn_range_pow]);
end

fig = stfig('DLD Front Panel: 2d count rate pofiles');

if mod(bins,2) == 0 %if the number of bins are even make them odd
    bins=bins+1;
    set(handles.spatial_bins_h,'String',int2str(bins));
end

if ~length(handles.txy_data_windowed)==0
    XEdges=linspace(xmin,xmax,bins);
    YEdges=linspace(ymin,ymax,bins);
    TEdges=linspace(tmin,tmax,bins);

    panes=sum([XY_bool,XT_bool,YT_bool]);
    pane_counter=1;

    stfig(fig);
    set(gcf,'Units','normal')
    set(gca,'Position',[0 0 1 1])
    if XY_bool
        subplot(1,panes,pane_counter);
        pane_counter=pane_counter+1;
        bin_area=((xmax-xmin)/bins)*((ymax-ymin)/bins);
        [counts,centers]=hist3(handles.txy_data_windowed(:,2:3),'edges',{XEdges,YEdges});
        counts=counts/bin_area;
        if  ~spatial_blur==0
            counts=imgaussfilt(counts,spatial_blur);
        end
        imagesc(10^3*centers{1},10^3*centers{2},transpose(counts))
        colormap(gca,cmap)
        set(gca,'Ydir','normal')
        set(gcf,'Color',[1 1 1]);
        title('Spatial Dist. TOP')
        xlabel('X(mm)')
        ylabel('Y(mm)')
        h=colorbar;
        xlabel(h,'Count Density (m^{-2})')

    end
    if XT_bool
        subplot(1,panes,pane_counter);
        pane_counter=pane_counter+1;
        bin_area=((xmax-xmin)/bins)*((tmax-tmin)/bins);
        [counts,centers]=hist3(handles.txy_data_windowed(:,1:2),'edges',{TEdges,XEdges});
        counts=counts/bin_area; 
        if  ~spatial_blur==0
            counts=imgaussfilt(counts,spatial_blur);
        end
        imagesc(10^3*centers{2},centers{1},counts)
        colormap(gca,cmap)
        set(gca,'Ydir','normal')
        set(gcf,'Color',[1 1 1]);
        title('Spatial Dist. XT')
        xlabel('X(mm)')
        ylabel('T(s)')
        h=colorbar;
        if range_comp_bool
            xlabel(h,sprintf('Count Density^{%.2f} (m^{-1}s^{-1})',dyn_range_pow))
        else
            xlabel(h,'Count Density (m^{-1}s^{-1})')
        end
    end
    if YT_bool
        subplot(1,panes,pane_counter);
        pane_counter=pane_counter+1;
        bin_area=((ymax-ymin)/bins)*((tmax-tmin)/bins);
        [counts,centers]=hist3(handles.txy_data_windowed(:,[1,3]),'edges',{TEdges,YEdges});
        counts=counts/bin_area;
        if  ~spatial_blur==0
            counts=imgaussfilt(counts,spatial_blur);
        end
        imagesc(10^3*centers{2},centers{1},counts)
        colormap(gca,cmap)
        set(gca,'Ydir','normal')
        set(gcf,'Color',[1 1 1]);
        title('Spatial Dist. YT')
        xlabel('Y(mm)')
        ylabel('T(s)')
        h=colorbar;
        if range_comp_bool
            xlabel(h,sprintf('Count Density^{%.2f} (m^{-1}s^{-1})',dyn_range_pow))
        else
            xlabel(h,'Count Density (m^{-1}s^{-1})')
        end
    end

     %colormap(hot);

    % 
    % binned_output = binned_output';
    % 
    % %%%%%%%%%% Plotting atoms on DLD %%%%%%%%%%
    % size(binned_output);
    % 
    % pcolor(binned_output)
    % colormap gray
    % colorbar
    % 
    %     xwidth_pix = 'Not fitted';
    %     ywidth_pix = 'Not fitted';
    % 
    % %binned_output = binned_output';

    if fit_flag == 1
        guess = [1, radius_fit_guess_num, radius_fit_guess_num, 0, y_fit_guess_num, x_fit_guess_num];

        if isnan(radius_fit_guess_num)
            guess(2) = 5*20/bins_per_pixel;
            guess(3) = 5*20/bins_per_pixel;
        end

        if isnan(x_fit_guess_num)
            guess(5) = 125*20/bins_per_pixel;
        end    

        if isnan(y_fit_guess_num)
            guess(6) = 100*20/bins_per_pixel;
        end  

        size(binned_output)

        xminbin = max(ceil(guess(5)-10*guess(2)),1);
        xmaxbin = min(ceil(guess(5)+10*guess(2)),min(size(binned_output)));
        yminbin = max(ceil(guess(6)-10*guess(3)),1);
        ymaxbin = min(ceil(guess(6)+10*guess(3)),min(size(binned_output)));

        binned_output = binned_output(xminbin:xmaxbin,yminbin:ymaxbin);

        sfigure(figure_number+3e6)
        pcolor(binned_output)
        colormap gray
        colorbar

        guess(5) = ceil(max(size(binned_output))/2);
        guess(6) = ceil(max(size(binned_output))/2);
        brightest_pixel = max(max(binned_output));
        binned_output = binned_output./brightest_pixel;

        [Param_out,l2error] = fminsearch(@(Param)fit_2d_gaussian(Param,binned_output),guess);
        fitted_profile = zeros(max(size(binned_output)),max(size(binned_output)));
        for xpos = 1:max(size(binned_output))
            for ypos = 1:max(size(binned_output))
                fitted_profile(xpos,ypos) = Param_out(1)*exp(-0.5*((xpos-Param_out(5))^2/(Param_out(2))^2))*exp(-0.5*((ypos-Param_out(6))^2/(Param_out(3))^2)) + Param_out(4);
            end
        end
        figure(figure_number+1e6)
        pcolor(fitted_profile)
        colormap gray
        colorbar
        xwidth_pix = Param_out(2);
        ywidth_pix = Param_out(3);

    end

else
    set(handles.status_handle ,'String',['No Counts in 2d Window']);
end