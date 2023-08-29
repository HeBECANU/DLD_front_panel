function varargout = dld_front_panel(varargin)
% DLD_FRONT_PANEL M-file for dld_front_panel.fig
%      DLD_FRONT_PANEL, by itself, creates a new DLD_FRONT_PANEL or raises the existing
%      singleton*.
% 
%      H = DLD_FRONT_PANEL returns the handle to a new DLD_FRONT_PANEL or the handle to
%      the existing singleton*.
% 
%      DLD_FRONT_PANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DLD_FRONT_PANEL.M with the given input arguments.
% 
%      DLD_FRONT_PANEL('Property','Value',...) creates a new DLD_FRONT_PANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dld_front_panel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dld_front_panel_OpeningFcn via
%      varargin.
%      
%      Files which dld_front_panel needs to run are: fileExists, dld_read_5channels_reconst_multi_a,
%      dld_tof_plotter_window_a, guassfun_a, tof_gauss_dld_a and dld_2dplotter_window
% 
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

% 
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dld_front_panel

% Last Modified by GUIDE v2.5 29-Jun-2018 15:42:39



%Improvement Log--------------------------------------------------------------------------
%v9 done
%fix bug if no files imported
%error message if not enough counts in tof to plot

%V9 plans
%import all with cashing
%update import script
%monitor indiciaror
%robust,faster monitoring
%constants file


%V8
%added in hold functionality for FFT
%correct normalization of fft with number of files

%V7:
%Added perceptually uniform colormaps 
%Fixed logscale blurring

%Known BUGS/ possible Improvements--------------------------------------------------------
%ocasional runtime bug with monitor
%fix TOF winowing for no counts case/specify bins for histogram
%idiot proof min max values
%faster TF fit
%the monitor indicator ... is broken
%catch zero counts to histogram
%import all option or incremental import
%clean up code to use more struct variables



%while drawnow update prevents focus stealing you cant use the buttons! for
%example to stop the import





% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @dld_front_panel_OpeningFcn, ...
    'gui_OutputFcn',  @dld_front_panel_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before dld_front_panel is made visible.
function dld_front_panel_OpeningFcn(hObject, eventdata, handles, varargin)

%BEGIN USER VAR-----------------------------------------------------------
%These are some tweaks for the program that are a bit to complex for the
%gui

%add all subfolders to the path
addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path
set_up_project_path('.', {'dev','lib','bin'})


hebec_constants 

handles.files_imported=0;
handles.falldist = .848; %fall distance of 848mm
handles.falltime = 0.416; %fall time of 416ms
handles.trapfreqrad=500;
handles.trapfreqaxial=50;
handles.grav=9.7961; %from https://d28rz98at9flks.cloudfront.net/12293/Rec1969_034.pdf
handles.boltzconst=1.38064852*10^-23; %J/K
handles.hbar=1.0545718*10^(-34);
handles.masshe=6.6464*10^-27; %mass in kg
handles.switchoff=22.954;
handles.dldtrig=20.3;
handles.auto_press_TOF=1;   %runs the Plot TOF button when read data is done
handles.auto_press_2d=1;    %same but for 2d plot 
handles.auto_press_3d=0;    %same but for 3d plot 
handles.auto_press_avgpos=0;   

%END USER VAR-------------------------------------------------------------



% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dld_front_panel (see VARARGIN)
% Choose default command line output for dld_front_panel
handles.output = hObject;
handles.txy_data = [0,0,0];
handles.num_hits = 0;   %global variable, as opposed to handles.num_hits_handle which is the handle to the edit text
handles.temp_fit = 0;
handles.spatial_fit = 0;
handles.free_run_mon_bool = 0;
handles.tof_distribution = 0;   %global variable corresponding to checkbox which determines which tof fit to use
handles.hold_Zoom = 0;

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes dld_front_panel wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function filename_load_Callback(hObject, eventdata, handles)
% hObject    handle to filename_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename_load as text
input = get(hObject,'String');

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
    set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function filename_load_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Outputs from this function are returned to the command line.
function varargout = dld_front_panel_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1 (read data button).
function pushbutton1_Callback(hObject, eventdata, handles)
%disableButtons(handles);
%pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update

if get(handles.pushbutton1,'UserData')
    set(handles.pushbutton1,'String','Stop Import')
    set(handles.pushbutton1,'UserData',0);
    importdata(hObject,eventdata,handles)
else
    disp('import halted');
    set(handles.pushbutton1,'String','Stopping')
    set(handles.pushbutton1,'UserData',1);
    
end


enableButtons(handles);

function importdata(hObject,eventdata,handles)
filename_input = get(handles.filename_load,'String');

% filename_with_ext = [filename_input,'.txt'];
% bool_fileExists = fileExists(filename_with_ext);   %This is a non-standard Matlab function, so must include
% if ~bool_fileExists   %if the read file doesn't exist, reset buttons and return function
%     enableButtons(handles);
%     return
% end

three_channel_output = [];
num_files_s = get(handles.num_files_h,'String');
num_files = str2double(num_files_s);

start_file_s = get(handles.start_file_h,'String');
start_file = str2double(start_file_s);

files_imported=0;
number_hits = 0;
num_hits_mat = [];
three_channel_output = [];
low_files=0;
set(handles.num_low_files,'String',low_files);
refresh(dld_front_panel)
for n = 1:num_files
    current_file=start_file+n-1;
    no_file=0;%flag if there was not a file
    if ~get(handles.pushbutton1,'UserData') %checks if the stop button has been pressed
        current_file_s = num2str(current_file);
        filename_no_ext = [filename_input,current_file_s];
        filename_with_ext = [filename_input,current_file_s,'.txt'];

        set(handles.status_handle ,'String',['Reading File',current_file_s]);
        %drawnow update; %this updates without stealing focus (like pause(1e-5))
        pause(1e-6)%this updates without stealing focus but allows button press unlike drawnow update
        if get(handles.TXY_checkbox,'Value')
            %if the txy does not exist make it
            if ~fileExists([filename_input,'_txy_forc',current_file_s,'.txt']) &&...
                fileExists(filename_with_ext)
                disp('txy does not exist, converting')
                dld_raw_to_txy(filename_input,current_file,current_file);
            %if both exists and the txt is newer than (converted) txy then
            %reconvert
            elseif fileExists([filename_input,'_txy_forc',current_file_s,'.txt']) &&...
                fileExists(filename_with_ext)

                %get the file date for both
                FileInfo=dir(filename_with_ext);
                date_txt=FileInfo.datenum;
                FileInfo=dir([filename_input,'_txy_forc',current_file_s,'.txt']);
                date_txy=FileInfo.datenum;

                %proces data to txt if the data newer than the coresponding txy
                if date_txy<date_txt
                    disp('txt is newer than txy, reconverting')

                    dld_raw_to_txy(filename_input,current_file,current_file);
                end

            end

            %catch if the TXY does not exist
            if ~fileExists([filename_input,'_txy_forc',current_file_s,'.txt'])
                fprintf('file number %i does not exist \n',current_file)
                no_file=1;
                %enableButtons(handles);
                %return
            else
                three_channel_output_single=txy_importer(filename_input,current_file_s);
            end
        else
            if fileExists(filename_with_ext)
                three_channel_output_single=dld_read_5channels_reconst_multi_imp([filename_input,current_file_s],1,0,1,0);
            else
                fprintf('file number %i does not exist \n',current_file)
                no_file=1;
                three_channel_output_single=[];
            end
        end
        if ~no_file %is there a file to process
            number_hits_single=size(three_channel_output_single,1);
            if number_hits_single < 0%1000   %return sub 100 hits in file
                disp('zero counts encountered')
                low_files=low_files+1;
                set(handles.num_low_files,'String',low_files);
            else
                files_imported=files_imported+1;
            end
            number_hits = number_hits + number_hits_single;
            num_hits_mat = cat(1,num_hits_mat,number_hits_single);
            if length(num_hits_mat)>1
                set(handles.num_hits_avg ,'String',mean(num_hits_mat));
                set(handles.num_hits_SD ,'String',std(num_hits_mat));
            end
            rot_angle = str2double(get(handles.rot_angle_h,'String'));
            three_channel_output_sorted_rot = [];
            if rot_angle ~=0
                three_channel_output_raw = three_channel_output_single;
                sin_theta = sin(rot_angle);
                cos_theta = cos(rot_angle);
                three_channel_output_sorted_rot(:,1) = three_channel_output_raw(:,1);
                three_channel_output_sorted_rot(:,2) = three_channel_output_raw(:,2)*cos_theta - three_channel_output_raw(:,3)*sin_theta;
                three_channel_output_sorted_rot(:,3) = three_channel_output_raw(:,2)*sin_theta + three_channel_output_raw(:,3)*cos_theta;
                three_channel_output_single = three_channel_output_sorted_rot;
            end
            three_channel_output{n}=three_channel_output_single;
            set(handles.status_handle ,'String','Idle');
            set(handles.num_hits_handle ,'String',num2str(number_hits));
            pause(1e-5)%update screen
        end%is there a file to process
    end%stop condition
end%looping over files


%% remove hot spots 2019-11-20
% first jam into the standard data stucutre 
% the above code should really be heavily modified to use our standard data structure
if get(handles.hotspot_checkbox,'Value') && files_imported~=0
    data=[]; 
    data.mcp_tdc.counts_txy=three_channel_output;
    data.mcp_tdc.num_counts=cellfun(@(x) size(x,1),data.mcp_tdc.counts_txy);
    data=hotspot_mask(data);
    three_channel_output=data.mcp_tdc.masked.counts_txy;
end

%% cat all the counts together
if files_imported~=0
    %need to handle the fail case better
    three_channel_output=vertcat(three_channel_output{:});
else
    three_channel_output=[];
end
assignin('base','numhits',num_hits_mat); %asigns to the main matlab workspace

handles.files_imported=files_imported; 
%%

handles.txy_data = three_channel_output;
handles.num_hits = number_hits;

guidata(hObject, handles);

if handles.auto_press_TOF   
    tof_button_Callback(hObject, eventdata, guidata(hObject))
end
if handles.auto_press_2d
  plot_2d_button_Callback(hObject, eventdata, guidata(hObject))
end
if handles.auto_press_avgpos
    Avg_Pos_Callback(hObject, eventdata, guidata(hObject))
end
if handles.auto_press_3d %make this one last so that the 3d rotate can work
    plot_3d_button_Callback(hObject, eventdata, guidata(hObject)) 
end

set(handles.pushbutton1,'UserData',1); %gets the button ready to run again
set(handles.pushbutton1,'String','Read Data')

function disableButtons(handles)
%change the mouse cursor to an hourglass
%set(handles.figure1,'Pointer','watch');

%disable all the buttons so they cannot be pressed
set(handles.pushbutton1,'Enable','off');

function enableButtons(handles)
%change the mouse cursor to an arrow
%set(handles.figure1,'Pointer','arrow');

%enable all the buttons so they can be pressed
set(handles.pushbutton1,'Enable','on');
% --- Executes on button press in tof_button.
function tof_button_Callback(hObject, eventdata, handles)
if handles.num_hits <10
    set(handles.status_handle ,'String','Not Enough Counts to Plot TOF');
    pause(1e-5)
    return
end

set(handles.plot_3d_button,'UserData',1)  %stop 3d rotate
set(handles.status_handle ,'String','Plotting 1D profiles');
pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update

windowdata(hObject,handles,1)
handles=guidata(hObject); %update the local copy of the handles to include the new handles.txy_data_windowed
set(handles.num_hits_tof_h ,'String',int2str(size(handles.txy_data_windowed,1)));

dld_plot_tof(hObject,handles)

set(handles.status_handle ,'String','Idle');
pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update

guidata(hObject, handles);

function time_binsize_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of time_binsize as text
%        str2double(get(hObject,'String')) returns contents of time_binsize as a double
bin_size = get(hObject,'String');


%checks to see if input is empty. if so, default bin size to .001
if (isempty(bin_size))
    set(hObject,'String','0.001')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function time_binsize_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function t_window_max_h_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of t_window_max_h as text
%        str2double(get(hObject,'String')) returns contents of t_window_max_h as a double
input = get(hObject,'String');


%checks to see if input is empty. if so, default bin size to .001
if (isempty(input))
    set(hObject,'String','2')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function t_window_max_h_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function t_window_min_h_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of t_window_min_h as text
%        str2double(get(hObject,'String')) returns contents of t_window_min_h as a double

input = get(hObject,'String');


%checks to see if input is empty. if so, default bin size to .001
if (isempty(input))
    set(hObject,'String','0')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function t_window_min_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in temp_fit_checkbox.
function temp_fit_checkbox_Callback(hObject, eventdata, handles)
checkboxStatus = get(handles.temp_fit_checkbox,'Value');
if (checkboxStatus)
    handles.temp_fit = 1;
else
    handles.temp_fit = 0;
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in plot_2d_button.
function plot_2d_button_Callback(hObject, eventdata, handles)
if handles.num_hits <10
    return
end
%guidata(hObject)
plot_2d(hObject,handles)

% bins_per_pix_s = get(handles.spatial_bins,'String');
% bins_per_pix   = str2double(bins_per_pix_s);
% pixel_size_x = bins_per_pix*.013158;
% pixel_size_x_s = num2str(pixel_size_x);
guidata(hObject, handles);

% pixel_size_y = bins_per_pix*.013158;
% pixel_size_x_s = num2str(pixel_size_x);
% set(handles.pixel_size_x_h ,'String',pixel_size_x_s);
function t_min_2d_handle_Callback(hObject, eventdata, handles)
% hObject    handle to t_min_2d_handle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_min_2d_handle as text
%        str2double(get(hObject,'String')) returns contents of t_min_2d_handle as a double
input = get(hObject,'String');


%checks to see if input is empty. if so, default bin size to .001
if (isempty(input))
    set(hObject,'String','0')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function t_min_2d_handle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function t_max_2d_handle_Callback(hObject, eventdata, handles)
input = get(hObject,'String');


%checks to see if input is empty. if so, default bin size to .001
if (isempty(input))
    set(hObject,'String','2')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function t_max_2d_handle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ymin_h_Callback(hObject, eventdata, handles)

input = get(hObject,'String');


%checks to see if input is empty. if so, default bin size to .001
if (isempty(input))
    set(hObject,'String','-0.04')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ymin_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xmin_h_Callback(hObject, eventdata, handles)
input = get(hObject,'String');


%checks to see if input is empty. if so, default bin size to .001
if (isempty(input))
    set(hObject,'String','-0.04')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xmin_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xmax_h_Callback(hObject, eventdata, handles)
input = get(hObject,'String');


%checks to see if input is empty. if so, default bin size to .001
if (isempty(input))
    set(hObject,'String','0.04')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xmax_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function num_files_h_Callback(hObject, eventdata, handles)
input = get(hObject,'String');

%checks to see if input is empty. if so, default bin size to .001
if (isempty(input))
    set(hObject,'String','1')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function num_files_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function start_file_h_Callback(hObject, eventdata, handles)
input = get(hObject,'String');


%checks to see if input is empty. if so, default bin size to .001
if (isempty(input))
    set(hObject,'String','1')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function start_file_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function free_run_shots_Callback(hObject, eventdata, handles)
input = get(hObject,'String');


%checks to see if input is empty. if so, default bin size to .001
if (isempty(input))
    set(hObject,'String','100')
end

guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of free_run_shots as text
%        str2double(get(hObject,'String')) returns contents of free_run_shots as a double

% --- Executes during object creation, after setting all properties.
function free_run_shots_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles);

function monitorrunbutton_ButtonDownFcn(hObject, eventdata, handles)

guidata(hObject, handles);

% --- Executes on button press in monitorrunbutton.
function monitorrunbutton_Callback(hObject, eventdata, handles)



if ~get(handles.monitorrunbutton,'UserData')
    set(handles.monitorrunbutton,'String','Stop mon')
    set(handles.monitorrunbutton,'UserData',1);
    monitor_files(hObject,eventdata,handles)
    
else
    set(handles.monitorrunbutton,'String','Monitor Run')
    set(handles.monitorrunbutton,'UserData',0);
end


pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update
% enableButtons(handles);
guidata(hObject, handles);

function  monitor_files(hObject,eventdata,handles)


filename_input = get(handles.filename_load,'String');
slash_split_filename=strsplit(filename_input,'\');
file_name_nonum=slash_split_filename(end);
monitor_dir=filename_input(1:end-(length(file_name_nonum)+1));
 
%set(handles.monitorrunbutton,'Enable','off');

% disableButtons(handles);
 pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update


mon_txy_flag = get(handles.mon_txy_check_box,'Value');


%need to cut dir name to remove all but the name eg
%'C:\datadir\d'
%to 'C:\datadir\'
% 
% filename_no_ext = [filename_input,current_file_s];
% filename_with_ext = [filename_input,current_file_s,'.txt'];


%Sketch for txt version

%read file names/dates
%cut out those that dont match .txt
%find new/modified
%change path to match
Low_Count_Size=0.0;
wait_for_mod=1;

loop_num=1;
pause on
dir_init_content = dir(monitor_dir);
initial_file_names = {dir_init_content.name};
initial_file_dates={dir_init_content.date};
initial_filenames_dates=[ initial_file_names ; initial_file_dates ]';
%cut . and .. from the listings
initial_filenames_dates=initial_filenames_dates(3:end,:);
initial_file_names=initial_file_names(3:end);
%fprintf('\nMonitoring %s ', watching_dir);
mon_continue=1;

%when only looking at the txy then dont need to wait
if get(handles.mon_txy_check_box,'Value')==1
     wait_for_mod=4;
     Low_Count_Size=0;
end

while mon_continue
    dir_init_content = dir(monitor_dir);
    file_names = {dir_init_content.name};
    file_dates=  {dir_init_content.date};
    filenames_dates=[ file_names ; file_dates ]';
    %cut . and .. from the listings
    file_names=file_names(3:end);
    filenames_dates=filenames_dates(3:end,:);
    new_files = setdiff(file_names,initial_file_names);
    deleted_files=setdiff(initial_file_names,file_names);

    %remove new files from the filenames_dates array
    %dont need to do this for inital as should be the same
    filenames_dates_newcut=filenames_dates(~ismember(file_names,new_files),:);
    
    %cut deleted files from the initial array
    if isempty(file_names)
        error('NO DIRECTORY FOUND AT THIS LOCATION CHECK THE FILENAME BOX')
    end
    
    non_deleted_mask=~ismember(initial_filenames_dates(:,1),deleted_files);
    initial_filenames_dates=initial_filenames_dates(non_deleted_mask,:);

    if ~isequal(size(initial_filenames_dates),size(filenames_dates_newcut))
        fprintf(2,'error file array size does not agree\n')
    end
    if ~isequal(initial_filenames_dates(:,1),filenames_dates_newcut(:,1))
        fprintf(2,'file names do not agree perhaps the order has been mixed up\n')
        disp(initial_filenames_dates)
        disp(filenames_dates_newcut)
    end
    %catch the empty case
    if ~(size(initial_filenames_dates,1)==0 || size(filenames_dates,1)==0)
      modified=datenum(initial_filenames_dates(:,2))<datenum(filenames_dates_newcut(:,2));
    else
        modified=[];
    end
    %add the moded file names to the new_files
    new_files=[new_files, filenames_dates(modified,1)];

    %chop off txy% data and keep txt files
    if get(handles.mon_txy_check_box,'Value')==1
        new_files=new_files(cellfun(@(x) ~isempty(strfind(x,'_txy_forc')),new_files));
        %new_files=new_files(cellfun(@(x) ~isempty(regexp(x,strcat(file_name_nonum,'_txy_forc.*.txt'))),new_files));
    else
        new_files=new_files(cellfun(@(x) isempty(strfind(x,'_txy_forc')),new_files));
        new_files=new_files(cellfun(@(x) ~isempty(regexp(x,[file_name_nonum '.*.txt'])),new_files));
    end
    %new_files=new_files(cellfun(@(x) ~isempty(findstr('.txt',x)),new_files));

    initial_file_names = file_names;
    initial_filenames_dates=filenames_dates;
      
      
    if ~isempty(new_files)
        %dont do a format check if txy, i could make a diffrent condition
        %but this will be left for an improvement
        
            
        if ~(get(handles.mon_txy_check_box,'Value')==1)
            %now we check if the first few lines in the file are the right format
            %prealocate the pass list
            %pass_line_test=zeros(numel(new_files),1);
            pass_line_test=[];
            for k=1:numel(new_files) 
                %wait for the file to be written before checking it
                wait_till_written(monitor_dir,new_files{k},wait_for_mod)
                %new_files{k}
                FID=fopen(fullfile(monitor_dir,new_files{k}),'r');
                FirstLine=fgetl(FID);
                SecondLine=fgetl(FID);
                fclose(FID);
                %check that the line lengths are right with 1 comma
                pass_line_test(k)=FirstLine(1)=='5' && size(SecondLine,2)<=15 ...
                && size(FirstLine,2)<=15 && sum(SecondLine==',')==1 && sum(FirstLine==',')==1
            end

            %if there was one or greater file with the correct format play a
            %sound
%             if sum(pass_line_test)>0
%                 fs = 16000;
%                 t = 0:(1/fs):0.02;
%                 f = 2000;
%                 a = 0.2;
%                 y = a*sin(2*pi*f*t);ewe
%                 sound(y, fs);
%             end

            %this is a bit of a hack as i dont think cells can use the same
            %logical mask approach that matricies use
            %pass_line_test
            new_files=new_files(pass_line_test*1:numel(pass_line_test));
        end
        % deal with the new files
        %meddage with number of new files found
        fprintf('\n(%d) new files detected \n',numel(new_files));
        for k=1:numel(new_files)
            fprintf('Reading %s ...', new_files{k});
            %this must be set to the max timeout time on the TDC front
            %pannel, 5 seconds is rarely exceeded
            %find the creation time of this file to store in the params LOG

            %first i check if it has '_txy_forc' in the file name which
            %indicates that it is a converted file
            
            FileInfo=dir(fullfile(monitor_dir,new_files{k}));
            FileSize=FileInfo.bytes;
            FileSize=FileSize/(2^20);
            %here i check if the modification date is more than 1 second
            %old, this relies on this computers clock being close to the TDC
            %computer

            %first reformat the string to have the path and the file number
            
            %I read the file size, if it is too small then i will not update
            if FileSize>Low_Count_Size
                
                %this just seperates out the number at the end of the file
                %and changes the box on the interface to be equal
                filename=new_files(k);
                filename=fullfile(monitor_dir,new_files{k});%combine C:/dir/d123.txt
                filename=filename(1:end-4); %C:/dir/d123
                numpart=filename(end-5:end);  %/d123 so that can handle up to 10000
                numpart=regexp(numpart,'\d*','Match'); %give number component %/d123
                numpart=numpart{end};
                
                
                %if handles.mon_txy_bool
                
                
                set(handles.filename_load,'string',filename_input);
                set(handles.start_file_h,'string',numpart);
                set(handles.num_files_h,'string','1');
                
                pushbutton1_Callback(hObject, eventdata, guidata(hObject))
                fprintf(' Read \n');
            else
               
                fprintf(2,'\n file is too small(%2.3f MiB) will not update.\n',FileSize)
               
            end
        end
        
        fprintf('Monitoring %s ', monitor_dir);
        loop_num=1;
    else
        
        if mod(loop_num,4)==0
            pause(.2)
            fprintf('\b\b\b')
            loop_num=1;
        else
            pause(.1) %little wait animation
            fprintf('.')
            loop_num=loop_num+1;
        end
        
    end
    mon_continue=get(handles.monitorrunbutton,'UserData');
end

% --- Executes on button press in stop_button.
function stop_button_Callback(hObject, eventdata, handles)
set(handles.monitorrunbutton,'UserData',0);
set(handles.monitorrunbutton,'Enable','on');
guidata(hObject, handles);

% --- Executes on button press in mon_txy_check_box.
function mon_txy_check_box_Callback(hObject, eventdata, handles)
if get(handles.mon_txy_check_box,'Value')
    set(handles.TXY_checkbox,'Value',1);
end

function ymax_h_Callback(hObject, eventdata, handles)

input = get(hObject,'String');


%checks to see if input is empty. if so, default bin size to .001
if (isempty(input))
    set(hObject,'String','-0.04')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ymax_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_2d(hObject,handles)
set(handles.plot_3d_button,'UserData',1)  %stop 3d rotate
set(handles.status_handle ,'String','Plotting 2D profiles');


pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update
windowdata(hObject,handles,0);
handles=guidata(hObject); %update the local copy of the handles to include the new handles.txy_data_windowed
set(handles.num_hits_2d_h ,'String',int2str(size(handles.txy_data_windowed,1)));

pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update
%[xwidth_pix,ywidth_pix] =
dld_2dplotter_window_a(handles);

pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update
%pixel_size_x = bins_per_pix*.013158;
%pixel_size_x_s = num2str(pixel_size_x);
%set(handles.pixel_size_x_h ,'String',pixel_size_x_s);

if handles.spatial_fit == 1
    xwidth = num2str(xwidth_pix*pixel_size_x);
    ywidth = num2str(ywidth_pix*pixel_size_x);
    
    set(handles.text45 ,'String',xwidth);
    set(handles.text46 ,'String',ywidth);

end

set(handles.status_handle ,'String','Idle');
pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update


return

%point cloud
function plot_3d(hObject,handles)

cap_counts=200000;
windowdata(hObject,handles,0);
handles=guidata(hObject); %update the local copy of the handles to include the new handles.txy_data_windowed
set(handles.num_hits_2d_h ,'String',int2str(size(handles.txy_data_windowed,1)));
totcounts=length(handles.txy_data_windowed);
if totcounts>cap_counts
    decimated=1;
    dat=datasample(handles.txy_data_windowed,cap_counts,1,'Replace',false);
else
    decimated=0;
    dat=handles.txy_data_windowed;
end

if ~totcounts==0
    sfigure(260)
    scatter3(dat(:,1),dat(:,2),dat(:,3),1,'k.')
    %set(gcf,'Position',[200 50 800 600])
    set(gcf,'Color',[1 1 1]);
    axis vis3d;
    %axis equal;
    xlabel('time(s)')
    ylabel('X(m)')
    zlabel('y(m)')
    if decimated
        title({'Decimated Counts,',num2str(cap_counts/totcounts,2),' of total'});
    else
        title({'All Counts'});
    end
    
    if get(handles.plot_3d_rot_checkbox,'Value') 
    make_rot_timer(handles)
    else
       set(handles.plot_3d_button,'UserData',1)
       set(handles.plot_3d_button,'String','Point Cloud');
    end
    
    
    %if it gets stuck use delete(timerfind)
else
     set(handles.status_handle ,'String',['No Counts in 2d Window']);
end

return

%define timer
function make_rot_timer(handles)   
ViewZ=[-20,10;-110,10;-190,80;-290,10;-380,10];
frames=200; % number frames
temp_p=(frames-1)/(size(ViewZ,1)-1); % length of each interval
ViewZ_new=zeros(frames,2);
% space view angles, if needed
    for inis=1:(size(ViewZ,1)-1)
        ViewZ_new(round(temp_p*(inis-1)+1):round(temp_p*inis+1),:)=...
            [linspace(ViewZ(inis,1),ViewZ(inis+1,1),...
             round(temp_p*inis)-round(temp_p*(inis-1))+1).',...
             linspace(ViewZ(inis,2),ViewZ(inis+1,2),...
             round(temp_p*inis)-round(temp_p*(inis-1))+1).'];
    end
    
handles.ViewZ=ViewZ_new;
handles.ViewZframes=frames;
RotTimer=timer; 
RotTimer.TimerFcn = @(src,event)timerCallback(src,event,handles);
RotTimer.StartFcn = @initTimer;
RotTimer.Period   = 0.05;
RotTimer.TasksToExecute = Inf;
RotTimer.ExecutionMode  = 'fixedSpacing';
start(RotTimer);
%pause(10);
%delete(RotTimer);
return

function initTimer(src, event)
   src.UserData = 1;
   disp('initialised')
return

function timerCallback(src,event,handles)
    sfigure(260);  %use silent figure https://www.mathworks.com/matlabcentral/fileexchange/8919--smart--silent-figure
    view(handles.ViewZ(src.UserData,:));
    drawnow;
    src.UserData = src.UserData + 1;
    if src.UserData>handles.ViewZframes
        src.UserData=1;
    end
   %disp(src.UserData)
   %disp(get(handles.plot_3d_button,'UserData'))
   if get(handles.plot_3d_button,'UserData')
       stop(src);
        disp('3dplot rot canceled');
       set(handles.plot_3d_button,'String','Point Cloud');
       pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update
   end
return
 
function rot_angle_h_Callback(hObject, eventdata, handles)
input = get(hObject,'String');

%checks to see if input is empty. if so, default rotation angle is 0.93 rad
if (isempty(input))
    set(hObject,'String','0.93')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function rot_angle_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in num_hits_multi_files_h.
function num_hits_multi_files_h_Callback(hObject, eventdata, handles)
set(handles.num_hits_multi_files_h ,'String',hits_in_window_s);

% --- Executes during object creation, after setting all properties.
function num_hits_multi_files_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TF_radius_guesss_h_Callback(hObject, eventdata, handles)
input = get(hObject,'String');

%checks to see if input is empty. if so, default rotation angle is 0.93 rad
if (isempty(input))
    set(hObject,'String','0.0002')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TF_radius_guesss_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function T_guess_input_h_Callback(hObject, eventdata, handles)
input = get(hObject,'String');

%checks to see if input is empty. if so, default rotation angle is 0.93 rad
if (isempty(input))
    set(hObject,'String','1')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function T_guess_input_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in spatial_fit_checkbox.
function spatial_fit_checkbox_Callback(hObject, eventdata, handles)
checkboxStatus = get(handles.spatial_fit_checkbox,'Value');
if (checkboxStatus)
    handles.spatial_fit = 1;
else
    handles.spatial_fit = 0;
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in FFT_checkbox.
function FFT_checkbox_Callback(hObject, eventdata, handles)

% --- Executes on button press in TXY_checkbox.
function TXY_checkbox_Callback(hObject, eventdata, handles)
TXY_checkboxStatus = get(handles.TXY_checkbox,'Value');
if (TXY_checkboxStatus)
    handles.use_TXY = 1;
else
    handles.use_TXY = 0;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on key press with focus on t_window_max_h and none of its controls.
function t_window_max_h_KeyPressFcn(hObject, eventdata, handles)
currChar = get(handles.figure1,'CurrentCharacter');
   if isequal(currChar,char(13)) %char(13) == enter key
       tof_button_Callback(hObject, eventdata, guidata(hObject))
       %call the pushbutton callback
   end
   
guidata(hObject, handles);

%This will replot the TOF on an enter key pressed in the t_min window

% --- Executes on key press with focus on t_window_max_h and none of its controls.
function t_window_min_h_KeyPressFcn(hObject, eventdata, handles)

currChar = get(handles.figure1,'CurrentCharacter');
   if isequal(currChar,char(13)) %char(13) == enter key
       tof_button_Callback(hObject, eventdata, guidata(hObject))
       %call the pushbutton callback
   end
   
guidata(hObject, handles);


% --- Executes on key press with focus on time_binsize and none of its controls.
function time_binsize_KeyPressFcn(hObject, eventdata, handles)

currChar = get(handles.figure1,'CurrentCharacter');
   if isequal(currChar,char(13)) %char(13) == enter key
       tof_button_Callback(hObject, eventdata, guidata(hObject))
       %call the pushbutton callback
   end
guidata(hObject, handles);
   
% --- Executes during object creation, after setting all properties.
function monitorrunbutton_CreateFcn(hObject, eventdata, handles)


% --- Executes on key press with focus on filename_load and none of its controls.
function filename_load_KeyPressFcn(hObject, eventdata, handles)
currChar = get(handles.figure1,'CurrentCharacter');
   if isequal(currChar,char(13)) %char(13) == enter key
       pushbutton1_Callback(hObject, eventdata, guidata(hObject))
       %call the pushbutton callback
   end
   guidata(hObject, handles);


% --- Executes on key press with focus on start_file_h and none of its controls.
function start_file_h_KeyPressFcn(hObject, eventdata, handles)
currChar = get(handles.figure1,'CurrentCharacter');
   if isequal(currChar,char(13)) %char(13) == enter key
       pushbutton1_Callback(hObject, eventdata, guidata(hObject))
       %call the pushbutton callback
   end
guidata(hObject, handles);



function windowdata(hObject,handles,TOF)
% this takes the data and then windows it to prevent duplication
%could be smart and check if the window changes
if TOF
    t_window_min = str2double(get(handles.t_window_min_h,'String'));

    t_window_max = str2double(get(handles.t_window_max_h,'String'));
else
t_window_min = str2double(get(handles.t_min_2d_handle,'String'));
t_window_max = str2double(get(handles.t_max_2d_handle,'String'));

end
ymin = str2double(get(handles.ymin_h,'String'))/1000;
ymax = str2double(get(handles.ymax_h,'String'))/1000;
xmin = str2double(get(handles.xmin_h,'String'))/1000;
xmax = str2double(get(handles.xmax_h,'String'))/1000;

mask=handles.txy_data(:,1)>t_window_min & handles.txy_data(:,1)<t_window_max...
    & handles.txy_data(:,2)>xmin & handles.txy_data(:,2)<xmax ...
    & handles.txy_data(:,3)>ymin & handles.txy_data(:,3)<ymax ;

handles.txy_data_windowed=handles.txy_data(mask,:);
guidata(hObject, handles);
    
function spatial_bins_h_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function spatial_bins_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in XY_checkbox.
function XY_checkbox_Callback(hObject, eventdata, handles)

% --- Executes on button press in XT_checkbox.
function XT_checkbox_Callback(hObject, eventdata, handles)

% --- Executes on button press in YT_checkbox.
function YT_checkbox_Callback(hObject, eventdata, handles)

function velocity_h_Callback(hObject, eventdata, handles)
% hObject    handle to velocity_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of velocity_h as text
%        str2double(get(hObject,'String')) returns contents of velocity_h as a double


% --- Executes during object creation, after setting all properties.
function velocity_h_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Hold_zoom_callback(hObject, eventdata,handles)
    
    hold_checkboxStatus = get(handles.zoom_checkbox,'Value');
    if (hold_checkboxStatus)
        handles.hold_Zoom = 1;
    else
        handles.hold_Zoom = 0;
    end
    % Update handles structure
    guidata(hObject, handles);
    


% --- Executes on button press in plot_3d_button.
function plot_3d_button_Callback(hObject, eventdata, handles)

if get(handles.plot_3d_button,'UserData')
    set(handles.plot_3d_button,'String','Stop Rot')
    set(handles.plot_3d_button,'UserData',0);
    if handles.num_hits <10
        return
    end
    %guidata(hObject)
    plot_3d(hObject,handles)
else
    disp('import halted');
    set(handles.plot_3d_button,'String','Stopping')
    set(handles.plot_3d_button,'UserData',1);   %stop 3d rotate
end


% --- Executes on button press in FFT_log_checkbox.
function FFT_log_checkbox_Callback(hObject, eventdata, handles)



function spatial_blur_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function spatial_blur_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LogScale2d_checkbox.
function LogScale2d_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to LogScale2d_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LogScale2d_checkbox


% --- Executes on button press in oned_time_checkbox.
function oned_time_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to oned_time_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of oned_time_checkbox


% --- Executes on button press in oned_X_checkbox.
function oned_X_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to oned_X_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of oned_X_checkbox


% --- Executes on button press in oned_Y_checkbox.
function oned_Y_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to oned_Y_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of oned_Y_checkbox



function oned_spatial_bins_h_Callback(hObject, eventdata, handles)
% hObject    handle to oned_spatial_bins_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oned_spatial_bins_h as text
%        str2double(get(hObject,'String')) returns contents of oned_spatial_bins_h as a double


% --- Executes during object creation, after setting all properties.
function oned_spatial_bins_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oned_spatial_bins_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Avg_Pos.
function Avg_Pos_Callback(hObject, eventdata, handles)
% hObject    handle to Avg_Pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.status_handle ,'String','Plotting Avg Pos');
pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update

windowdata(hObject,handles,1);
handles=guidata(hObject); %update the local copy of the handles to include the new handles.txy_data_windowed
%set(handles.num_hits_2d_h ,'String',int2str(size(handles.txy_data_windowed,1)));
%pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update
set(handles.plot_3d_button,'UserData',1) %stop 3d rotate
%[xwidth_pix,ywidth_pix] =
dld_avg_pos_plots(handles,300);
pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update
set(handles.status_handle ,'String','Idle');
pause(1e-5)%this updates without stealing focus but allows button press unlike drawnow update
guidata(hObject, handles);


% --- Executes on button press in avg_pos_fft_checkbox.
function avg_pos_fft_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to avg_pos_fft_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of avg_pos_fft_checkbox


% --- Executes on button press in avg_pos_fft_log_checkbox.
function avg_pos_fft_log_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to avg_pos_fft_log_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of avg_pos_fft_log_checkbox


% --- Executes on button press in plot_3d_rot_checkbox.
function plot_3d_rot_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to plot_3d_rot_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_3d_rot_checkbox


% --- Executes on button press in oned_fit_checkbox.
function oned_fit_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to oned_fit_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of oned_fit_checkbox


% --- Executes on button press in fit_bimod_checkbox.
function fit_bimod_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to fit_bimod_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on t_min_2d_handle and none of its controls.
function t_min_2d_handle_KeyPressFcn(hObject, eventdata, handles)
currChar = get(handles.figure1,'CurrentCharacter');
   if isequal(currChar,char(13)) %char(13) == enter key
       plot_2d_button_Callback(hObject, eventdata, handles)
       %call the pushbutton callback
   end
guidata(hObject, handles);
   

% --- Executes on key press with focus on t_max_2d_handle and none of its controls.
function t_max_2d_handle_KeyPressFcn(hObject, eventdata, handles)
currChar = get(handles.figure1,'CurrentCharacter');
   if isequal(currChar,char(13)) %char(13) == enter key
       plot_2d_button_Callback(hObject, eventdata, handles)
       %call the pushbutton callback
   end
guidata(hObject, handles);


% --- Executes on key press with focus on tof_button and none of its controls.
function tof_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to tof_button (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
