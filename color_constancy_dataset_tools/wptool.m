%
% Copyright (c) 2019, Intel Corporation
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the Intel Corporation nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%

function varargout = wptool(varargin)
% WPTOOL MATLAB code for wptool.fig
%      WPTOOL(filename, cameraid) opens the white point annotation tool for
%      the given .plain16 raw Bayer file and cameraid. The .plain16 format
%      refers to two bytes per pixel uncompressed format without any
%      header, in which the pixel values can be read as uint16.
%
%      cameraid = {'SonyIMX135'|'SonyIMX135_BLCCSC'|'Canon5DSR'|'NikonD810'|'Canon5Dmk3'|'Canon5Dmk3Video'}
%
%      The tool will store .wp, .ccm, and .meta files with the same path
%      and filename body as in the given input file.
%
%      To add support for new cameraid:
%      1: Add new frame_format_[new_cameraid].mat with the same
%      structure as in the existing ones, to store the raw frame format
%      save('frame_format_[new_cameraid].mat', '-struct', 'frame_format', '-mat')
%      2: Add new reference_wps_ccms_[new_cameraid].mat with the
%      same structure as in the existing ones, to store the white points
%      and color conversion matrices for the new camera, for the given 10
%      lab light sources (5 Image Engineering LightSTUDIO and 5 x-rite
%      SpectraLight III light sources)
%      3: Add handling of the new cameraid in functions
%      wptool_read_camera_characterization() and wptool_read_raw_frame_properties()
%      Hint: Look for an existing cameraid to find the correct places in
%      the code.
%
%      Note that pedestal subtraction is supported by the tool, but if
%      color shading correction is needed, it has to be applied on the
%      .plain16 raw images prior to calling the tool (in which case also
%      the pedestal subtraction needs to be applied already on the raw
%      input images and zero pedestal level indicated in the frame format).
%
%      .wp file is a tab delimited ascii file and contains two values:
%          RperG BperG
%      To calculate WB gains from the white point:
%          wb_gains = [1/RperG, 1, 1/BperG];
%          mingain = min(wb_gains(:));
%          wb_gains = wb_gains./mingain; % To prevent <1.0 gains
%
%      .ccm file is a tab delimited ascii file and contains 9 values:
%          RinR GinR BinR RinG GinG BinG RinB GinB BinB
%      To convert the stored CCM into the 3x3 matrix:
%          ccm = reshape(ccm, 3, 3)'; % Note also the transpose
%
%      When adding new frame formats, note the bayer_order:
%          0 = [Gr, R, B, Gb]
%          1 = [R, Gr, Gb, B]
%          2 = [B, Gb, Gr, R]
%          3 = [Gb, B, R, Gr]
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wptool

% Last Modified by GUIDE v2.5 01-Nov-2018 12:12:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wptool_OpeningFcn, ...
                   'gui_OutputFcn',  @wptool_OutputFcn, ...
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

% Manual addition: List of the variables that have been added in handles
% handles.filename: the .plain16 file to be processed
% handles.cameraid: {'SonyIMX135'|'SonyIMX135_BLCCSC'|'Canon5DSR'|'NikonD810'|'Canon5Dmk3'|'Canon5Dmk3Video'}
% handles.whitepoint: The final WP
% handles.whitepoint_candidate: The latest WP candidate
% handles.ccm: The selected CCM as a vector (need to reshape to 3x3 before use)
% handles.ref_frame_format: Structure of the raw frame format for this camera
% handles.ref_wps_ccms: Structure of the reference WPs and CCMs for this camera
% handles.other_characterization: Other camera characterization for this camera
% handles.R: Downscaled R component
% handles.G: Downscaled G component
% handles.B: Downscaled B component
% handles.RperG_BperG: List of the [R/G,B/G] pairs that correspond to the current image
% handles.R_filt: The same as R, except noise reduced by 9x9 averaging kernel
% handles.G_filt: The same as G, except noise reduced by 9x9 averaging kernel
% handles.B_filt: The same as B, except noise reduced by 9x9 averaging kernel
% handles.metadata: exposure parameters etc. from EXIF or .i3av4 header stored to .meta file


% --- Executes just before wptool is made visible.
function wptool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wptool (see VARARGIN)

% Choose default command line output for wptool
handles.output = hObject;

% Manual addition: Initialize values
if nargin < 5 || nargin > 6,
    error('Usage: wptool(filename_str, camera_id, white_map_file), where camera_id = {''SonyIMX135''|''SonyIMX135_BLCCSC''|''Canon5DSR''|''NikonD810''|''Canon5Dmk3''|''Canon5Dmk3Video''}, and the file is .plain16 raw Bayer file');
end
handles.filename = varargin{2};
handles.text_filename.String = handles.filename;
handles.cameraid = varargin{1};
if nargin == 6
    fname_whitemap = varargin{3};
    if strcmp(fname_whitemap,'') == 1
        handles.whitemap = [];
    else
        handles.whitemap = load(fname_whitemap);
    end
else
    handles.whitemap = [];
end
handles = wptool_init_metadata(handles);
handles = wptool_read_camera_characterization(handles); % Checks also for .wp file
wptool_update_final_wp_edits(handles);
handles = wptool_read_raw_frame_properties(handles);
wptool_update_raw_frame_properties_table(handles);

% Update metadata file edit box
if exist('wptool_recent_metadata_path.txt', 'file') == 2
    fid = fopen('wptool_recent_metadata_path.txt','r');
    txt = textscan(fid,'%s','delimiter','\n');
    meta_path = txt{1};
    set(handles.edit_metadatapath, 'String', meta_path);
else
   meta_path = '';
   set(handles.edit_metadatapath, 'String', meta_path);
   fid = fopen('wptool_recent_metadata_path.txt','w');
   fprintf(fid, '%s', meta_path);
   fclose(fid);
end

% Read metadata
handles = wptool_read_metadata(handles);

% Other initializations
handles = wptool_init_plots(handles);
handles.gmb_chart_detected = 0;
handles = wptool_detect_gretag_macbeth(handles);
handles = wptool_update_metadata_listbox(handles);
handles = wptool_update_tags_checkboxes(handles);
handles.whitepoint_candidate = [0.0, 0.0];
handles.marked_whitepoints = {};
handles.previous_coordinates = [-1,-1; -1,-1; -1,-1; -1,-1];
handles.previous_coordinates_ind = 1;
wptool_update_plots(handles);

% If Gretag MacBeth chart was detected
if handles.gmb_chart_detected > 0
    reply = questdlg('Gretag MacBeth chart detected. Is it correct?', ...
                     'Gretag MacBeth chart detected', ...
                     'Yes','No','Yes');
    if strcmp(reply,'No') == 1
        % Reset the chart coordinates to invalid values
        handles.metadata.gmb.tl = [-1,-1];
        handles.metadata.gmb.tr = [-1,-1];
        handles.metadata.gmb.br = [-1,-1];
        handles.metadata.gmb.bl = [-1,-1];
        handles.metadata.tag_ground_truth = 0;
        handles.gmb_chart_detected = 0;
        
        handles = wptool_update_metadata_listbox(handles);
        handles = wptool_update_tags_checkboxes(handles);
        wptool_update_plots(handles);
    else
        % Calculate initial white point candidates based on the detected chart
        handles = wptool_add_grey_patches_from_gretag_macbeth(handles);
        handles = wptool_calculate_final_wp(handles);
        wptool_update_plots(handles);
    end
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes wptool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wptool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_wpadd.
function pushbutton_wpadd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_wpadd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: It would be better if the button group would be populated by
% buttons in run-time, based on contents of handles.ref_wps_ccms.wps_ccms.
% Then e.g. the tag of the button would directly indicate the entry in the
% struct, and there would not be any other buttons than the ones actually
% covered by the struct.

h = get(handles.btgrp_wps,'SelectedObject');
wp_data = get(handles.table_wps, 'Data');
sel_wp_name = h.Tag;
new_wp_line = [];
switch sel_wp_name
    case 'rb_sl_d_wp'
        new_wp_line = [wptool_find_wp_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_D'), 1];
    case 'rb_sl_cw_wp'
        new_wp_line = [wptool_find_wp_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_CW'), 1];
    case 'rb_sl_hor_wp'
        new_wp_line = [wptool_find_wp_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_Hor'), 1];
    case 'rb_sl_a_wp'
        new_wp_line = [wptool_find_wp_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_A'), 1];
    case 'rb_sl_tl84u30_wp'
        new_wp_line = [wptool_find_wp_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_TL84U30'), 1];
    case 'rb_ie_d65_wp'
        new_wp_line = [wptool_find_wp_by_name(handles.ref_wps_ccms.wps_ccms, 'IE_D65'), 1];
    case 'rb_ie_d50_wp'
        new_wp_line = [wptool_find_wp_by_name(handles.ref_wps_ccms.wps_ccms, 'IE_D50'), 1];
    case 'rb_ie_f11_wp'
        new_wp_line = [wptool_find_wp_by_name(handles.ref_wps_ccms.wps_ccms, 'IE_F11'), 1];
    case 'rb_ie_f12_wp'
        new_wp_line = [wptool_find_wp_by_name(handles.ref_wps_ccms.wps_ccms, 'IE_F12'), 1];
    case 'rb_ie_a_wp'
        new_wp_line = [wptool_find_wp_by_name(handles.ref_wps_ccms.wps_ccms, 'IE_A'), 1];
end % switch

if numel(wp_data) < 3
    wp_data = new_wp_line;
else
    wp_data = [wp_data; new_wp_line];
end
set(handles.table_wps, 'Data', wp_data);


function edit_wp_rperg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_wp_rperg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_wp_rperg as text
%        str2double(get(hObject,'String')) returns contents of edit_wp_rperg as a double


% --- Executes during object creation, after setting all properties.
function edit_wp_rperg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_wp_rperg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_wp_bperg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_wp_bperg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_wp_bperg as text
%        str2double(get(hObject,'String')) returns contents of edit_wp_bperg as a double


% --- Executes during object creation, after setting all properties.
function edit_wp_bperg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_wp_bperg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_wpcalc: "Calculate final WP"
function pushbutton_wpcalc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_wpcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = wptool_calculate_final_wp(handles);
guidata(hObject, handles);


function [handles_out] = wptool_calculate_final_wp(handles)

% Calculate the new white point as the weighted average of the white points
% in table_wps
data = get(handles.table_wps, 'Data');
datasize = size(data);
if datasize(2) == 3 && datasize(1) > 0
    weightsum = sum(data(:,3));
    data(:,1) = (data(:,1).*data(:,3))./weightsum;
    data(:,2) = (data(:,2).*data(:,3))./weightsum;
    rperg = sum(data(:,1));
    bperg = sum(data(:,2));
    
    % Store the new white point
    handles.whitepoint = [rperg, bperg];
    
    % Show the new white point
    wptool_update_final_wp_edits(handles);
    
    % Plot the [R/G,B/G] -data
    wptool_update_chrplot(handles);
end % if

handles_out = handles;

% --- Executes on button press in pushbutton_writewp.
function pushbutton_writewp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_writewp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

wptool_write_wp_file(handles);

% --- Executes on button press in pushbutton_writeccm.
function pushbutton_writeccm_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_writeccm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

wptool_write_ccm_file(handles);

% --- Executes on button press in pushbutton_updateplots: "Update illustration"
function pushbutton_updateplots_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_updateplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

wptool_update_plots(handles);

% Manual addition: Update the plots
function [handles_out] = wptool_init_plots(handles)

% Read the raw data and crop if needed
fname = handles.filename;
width = handles.ref_frame_format.width;
height = handles.ref_frame_format.height;
raw_data = wptool_plain16_read(fname, width, height); % Read the data from .plain16 file
crop_top = handles.ref_frame_format.crop_top;
crop_left = handles.ref_frame_format.crop_left;
crop_bottom = handles.ref_frame_format.crop_bottom;
crop_right = handles.ref_frame_format.crop_right;
raw_data = raw_data(crop_top+1:end-crop_bottom, crop_left+1:end-crop_right); % Crop out light shielded pixels etc. (if any)
bayer_order = handles.ref_frame_format.bayer_order;
[Gr, R, B, Gb] = wptool_separateChannels(raw_data, bayer_order);
G = (double(Gr)+double(Gb))./2;
R = double(R);
B = double(B);

% Optional rotation
if handles.ref_frame_format.rotate180 > 0
    R = R(end:-1:1, end:-1:1);
    G = G(end:-1:1, end:-1:1);
    B = B(end:-1:1, end:-1:1);
end

% Optional pedestal subtraction, stretch saturation point back if pedestal is subtracted
pedestal = handles.ref_frame_format.pedestal;
if pedestal > 0,
    R = max(0, R - pedestal);
    G = max(0, G - pedestal);
    B = max(0, B - pedestal);
    %bpp = handles.ref_frame_format.bpp;
    %maxval = (2^bpp)-1;
    maxval = handles.ref_frame_format.satpoint;
    extra_gain = maxval/(maxval-pedestal);
    R = R.*extra_gain;
    G = G.*extra_gain;
    B = B.*extra_gain;
end

% Downscale the image, full resolution is not needed (noise is also reduced due to the averaging)
R = imresize(R, [1080 nan], 'bilinear');
G = imresize(G, [1080 nan], 'bilinear');
B = imresize(B, [1080 nan], 'bilinear');

% Store the pre-processed raw image
handles.R = R;
handles.G = G;
handles.B = B;
RperG = round(R./G, 3);
BperG = round(B./G, 3);
handles.RperG_BperG = unique([RperG(:), BperG(:)], 'rows');
handles.R_filt = conv2(R, ones(9)/(9*9), 'same'); % Average inside 9x9 to reduce noise impact
handles.G_filt = conv2(G, ones(9)/(9*9), 'same'); % Average inside 9x9 to reduce noise impact
handles.B_filt = conv2(B, ones(9)/(9*9), 'same'); % Average inside 9x9 to reduce noise impact

handles_out = handles;


% Manual addition: Update the plots
function wptool_update_plots(handles)

% Display the image, using the current WP and CCM
axes(handles.fig_image);
wb_gains = wptool_convert_wp_to_gains(handles.whitepoint);
ccm = reshape(handles.ccm, 3, 3)';
tmp_img = wptool_illustrate_image( handles.R, handles.G, handles.B, wb_gains, ccm, handles.ref_frame_format.satpoint, 0.45 );
if min(structfun(@(x)min(x(:)),handles.metadata.gmb)) > 0 % If Gretag MacBeth chart location is available
    tmp_img = wptool_illustrate_gmb(tmp_img, handles.metadata.gmb);
end
h = imshow(tmp_img);
set(h, 'ButtonDownFcn', @wptool_get_mouse_position); % Register the callback function that is called when the displayed image is clicked

% Plot the [R/G,B/G] -data
wptool_update_chrplot(handles);


% Manual addition: Plot the [R/G,B/G] -data
function wptool_update_chrplot(handles)

axes(handles.fig_chrplot);
h = plot(handles.RperG_BperG(:,1), handles.RperG_BperG(:,2), '.', 'MarkerEdgeColor', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6]);
hold on;
xmin = 100;
xmax = 0;
ymin = 100;
ymax = 0;
for i=1:length(handles.ref_wps_ccms.wps_ccms),
    if strcmp(handles.ref_wps_ccms.wps_ccms(i).name(1:3),'IE_')
        plot(handles.ref_wps_ccms.wps_ccms(i).wp(1), handles.ref_wps_ccms.wps_ccms(i).wp(2), 'o', 'MarkerEdgeColor', [0.75 0 0.75]);
    else
        plot(handles.ref_wps_ccms.wps_ccms(i).wp(1), handles.ref_wps_ccms.wps_ccms(i).wp(2), 'bo');
    end
    xmin = min(handles.ref_wps_ccms.wps_ccms(i).wp(1), xmin);
    xmax = max(handles.ref_wps_ccms.wps_ccms(i).wp(1), xmax);
    ymin = min(handles.ref_wps_ccms.wps_ccms(i).wp(2), ymin);
    ymax = max(handles.ref_wps_ccms.wps_ccms(i).wp(2), ymax);
end
if ~isempty(handles.whitemap)
    plot(handles.whitemap(:,1), handles.whitemap(:,2), 'k-', 'LineWidth', 1);
end
plot(handles.whitepoint(1), handles.whitepoint(2), 'ro', 'LineWidth', 2);
xmin = min(handles.whitepoint(1), xmin);
xmax = max(handles.whitepoint(1), xmax);
ymin = min(handles.whitepoint(2), ymin);
ymax = max(handles.whitepoint(2), ymax);
if handles.whitepoint_candidate(1) > 0 && handles.whitepoint_candidate(2) > 0
    plot(handles.whitepoint_candidate(1), handles.whitepoint_candidate(2), 'r*');
    xmin = min(handles.whitepoint_candidate(1), xmin);
    xmax = max(handles.whitepoint_candidate(1), xmax);
    ymin = min(handles.whitepoint_candidate(2), ymin);
    ymax = max(handles.whitepoint_candidate(2), ymax);
end
if ~isempty(handles.marked_whitepoints)
    for i = 1:length(handles.marked_whitepoints)
        chr = handles.marked_whitepoints{i};
        plot(chr(1), chr(2), 'mx');
        xmin = min(chr(1), xmin);
        xmax = max(chr(1), xmax);
        ymin = min(chr(2), ymin);
        ymax = max(chr(2), ymax);
    end
end
hold off;
xmin = xmin*0.5;
xmax = xmax*1.5;
ymin = ymin*0.5;
ymax = ymax*1.5;
axis([xmin xmax ymin ymax]);
xlabel('R/G');
ylabel('B/G');
set(h, 'ButtonDownFcn', @wptool_get_mouse_position_chrplot); % Register the callback function that is called when the plot is clicked
set(handles.fig_chrplot, 'ButtonDownFcn', @wptool_get_mouse_position_chrplot); % Register the callback function that is called when the plot is clicked


% Manual addition: Update the final WP edit boxes
function wptool_update_final_wp_edits(handles)

% Update Final WP edit boxes
wp = handles.whitepoint;
handles.text_final_rperg.String = num2str(wp(1), 4);
handles.text_final_bperg.String = num2str(wp(2), 4);


% Manual addition: Read the WPs and CCMs and other characterization that
% correspond to the used camera, and set the initial final WP
% (read from .wp if it exists)
function [handles_out] = wptool_read_camera_characterization(handles)

if strcmp('SonyIMX135', handles.cameraid)
    handles.ref_wps_ccms = load('data\reference_wps_ccms_sonyimx135.mat');
    handles.other_characterization = load('data\other_characterization_sonyimx135.mat');
    handles.metadata.aperture = 2.4; % Aperture is constant
elseif strcmp('SonyIMX135_BLCCSC', handles.cameraid)
    handles.ref_wps_ccms = load('data\reference_wps_ccms_sonyimx135.mat');
    handles.other_characterization = load('data\other_characterization_sonyimx135.mat');
    handles.metadata.aperture = 2.4; % Aperture is constant
elseif strcmp('Canon5DSR', handles.cameraid)
    handles.ref_wps_ccms = load('data\reference_wps_ccms_canon5dsr.mat');
    handles.other_characterization = load('data\other_characterization_canon5dsr.mat');
elseif strcmp('NikonD810', handles.cameraid)
    handles.ref_wps_ccms = load('data\reference_wps_ccms_nikond810.mat');
    handles.other_characterization = load('data\other_characterization_nikond810.mat');
elseif strcmp('Canon5Dmk3', handles.cameraid)
    handles.ref_wps_ccms = load('data\reference_wps_ccms_canon5dmk3.mat');
    handles.other_characterization = load('data\other_characterization_canon5dmk3.mat');
elseif strcmp('Canon5Dmk3Video', handles.cameraid)
    handles.ref_wps_ccms = load('data\reference_wps_ccms_canon5dmk3.mat');
    handles.other_characterization = load('data\other_characterization_canon5dmk3.mat');
else
    error('wptool: Incorrect camera ID');
end

% Check if .wp file already exists for this raw file, and use it if available
fname_wp = wptool_replace_file_extension(handles.filename, '.wp');
if exist(fname_wp, 'file') == 2
    wp = load(fname_wp);
    handles.whitepoint = wp;
    set(handles.table_wps, 'Data', [wp, 1]);
else
    set(handles.table_wps, 'Data', []);
    handles.whitepoint = [1.0, 1.0];
end

% Check if .ccm file already exists for this raw file, and use it if available
fname_ccm = wptool_replace_file_extension(handles.filename, '.ccm');
if exist(fname_ccm, 'file') == 2
    handles.ccm = load(fname_ccm);
    set(handles.btgrp_ccms,'SelectedObject',handles.rb_current_ccm_file);
else
    %handles.ccm = [1, 0, 0, 0, 1, 0, 0, 0, 1];
    set(handles.btgrp_ccms,'SelectedObject',handles.rb_sl_d_ccm);
    handles.ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_D');
end
handles.text_ccm.String = num2str(reshape(handles.ccm,3,3)', 4);

handles_out = handles;

% Manual addition: Read the raw frame properties that correspond to the used camera
function [handles_out] = wptool_read_raw_frame_properties(handles)

if strcmp('SonyIMX135', handles.cameraid)
    handles.ref_frame_format = load('data\frame_format_sonyimx135.mat');
elseif strcmp('SonyIMX135_BLCCSC', handles.cameraid)
    handles.ref_frame_format = load('data\frame_format_sonyimx135_blccsc.mat');
elseif strcmp('Canon5DSR', handles.cameraid)
    handles.ref_frame_format = load('data\frame_format_canon5dsr.mat');
elseif strcmp('NikonD810', handles.cameraid)
    handles.ref_frame_format = load('data\frame_format_nikond810.mat');
elseif strcmp('Canon5Dmk3', handles.cameraid)
    handles.ref_frame_format = load('data\frame_format_canon5dmk3.mat');
elseif strcmp('Canon5Dmk3Video', handles.cameraid)
    handles.ref_frame_format = load('data\frame_format_canon5dmk3_video.mat');
else
    error('wptool: Incorrect camera ID');
end

handles.text_cameraid.String = handles.cameraid;
handles_out = handles;


% Manual addition: Update the table that shows the frame format parameters
function wptool_update_raw_frame_properties_table(handles)

names = fieldnames(handles.ref_frame_format);
data = zeros(numel(names),1);
for i = 1:numel(names),
    data(i) = handles.ref_frame_format.(names{i});
end

set(handles.table_raw_frame_properties, 'Data', data);
set(handles.table_raw_frame_properties, 'RowName', names);


% Manual addition: Write .wp file
function wptool_write_wp_file(handles)

wp = handles.whitepoint;
fname_wp = wptool_replace_file_extension(handles.filename, '.wp');

% If the file already exists, confirm from the user whether it is ok to proceed
proceed = 1;
if exist(fname_wp, 'file') == 2
    reply = questdlg('The .wp file already exists. Overwrite?', ...
                     'The .wp file already exists', ...
                     'Yes','No','No');
    if strcmp(reply,'No') == 1
        proceed = 0;
    end
end % if

if proceed == 1
    save(fname_wp, 'wp', '-ascii', '-tabs');
    disp('.wp file written');
end % if

% Save the white point in a common file of 10 most recent white points
if exist('wptool_recent_wps.wp', 'file') == 2
    wplist = load('wptool_recent_wps.wp');
    wplist_sz = size(wplist);
    if wplist_sz(1) > 10,
        wplist = wplist(end-9:end,:);
        save('wptool_recent_wps.wp', 'wplist', '-ascii', '-tabs');
        disp('.wp list file truncated');
    end
end % if
save('wptool_recent_wps.wp', 'wp', '-ascii', '-tabs', '-append');
disp('.wp list file updated');


% Manual addition: Write .ccm file
function wptool_write_ccm_file(handles)

ccm = handles.ccm;
fname_ccm = wptool_replace_file_extension(handles.filename, '.ccm');

% If the file already exists, confirm from the user whether it is ok to proceed
proceed = 1;
if exist(fname_ccm, 'file') == 2
    reply = questdlg('The .ccm file already exists. Overwrite?', ...
                     'The .ccm file already exists', ...
                     'Yes','No','No');
    if strcmp(reply,'No') == 1
        proceed = 0;
    end
end % if

if proceed == 1
    save(fname_ccm, 'ccm', '-ascii', '-tabs');
    disp('.ccm file written');
end % if


% Manual addition: Initialize metadata
function [handles_out] = wptool_init_metadata(handles)

handles.metadata.exposure_time = -1;
handles.metadata.iso = -1;
handles.metadata.analog_gain = -1;
handles.metadata.digital_gain = -1;
handles.metadata.aperture = -1;
handles.metadata.lux = -1;
handles.metadata.lux_estimate = -1;
handles.metadata.gmb.tl = [-1,-1];
handles.metadata.gmb.tr = [-1,-1];
handles.metadata.gmb.br = [-1,-1];
handles.metadata.gmb.bl = [-1,-1];
handles.metadata.tag_ground_truth = 0;
handles.metadata.tag_mixed_illumination = 0;
handles.metadata.normalized_exposure_s = -1;

handles_out = handles;


% Manual addition: Read metadata from .meta file
function [handles_out] = wptool_read_metadata(handles)

fname_meta = wptool_replace_file_extension(handles.filename, '.meta');

if exist(fname_meta, 'file') == 2
    fid = fopen(fname_meta,'r');
    meta_struct = textscan(fid,'%s%f');
    fclose(fid);
    for i = 1:length(meta_struct{1})
        if strcmp(meta_struct{1}(i),'exposure_time')
            handles.metadata.exposure_time = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'iso')
            handles.metadata.iso = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'analog_gain')
            handles.metadata.analog_gain = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'digital_gain')
            handles.metadata.digital_gain = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'aperture')
            handles.metadata.aperture = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'lux')
            handles.metadata.lux = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'lux_estimate')
            handles.metadata.lux_estimate = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'gmb_tl_x')
            handles.metadata.gmb.tl(1) = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'gmb_tl_y')
            handles.metadata.gmb.tl(2) = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'gmb_tr_x')
            handles.metadata.gmb.tr(1) = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'gmb_tr_y')
            handles.metadata.gmb.tr(2) = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'gmb_br_x')
            handles.metadata.gmb.br(1) = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'gmb_br_y')
            handles.metadata.gmb.br(2) = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'gmb_bl_x')
            handles.metadata.gmb.bl(1) = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'gmb_bl_y')
            handles.metadata.gmb.bl(2) = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'tag_ground_truth')
            handles.metadata.tag_ground_truth = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'tag_mixed_illumination')
            handles.metadata.tag_mixed_illumination = meta_struct{2}(i);
        elseif strcmp(meta_struct{1}(i),'normalized_exposure_s')
            handles.metadata.normalized_exposure_s = meta_struct{2}(i);
        end
    end % for i
end % if

handles_out = handles;


% Manual addition: Look for Gretag MacBeth chart in the image
function [handles_out] = wptool_detect_gretag_macbeth(handles)

if handles.metadata.gmb.tl(1) < 0 % No valid chart coordinates was stored in .meta
    % Construct the Bayer GR_BG frame that is supported by the detection function
    raw = zeros(size(handles.R)*2);
    raw(1:2:end,1:2:end) = handles.G;
    raw(1:2:end,2:2:end) = handles.R;
    raw(2:2:end,1:2:end) = handles.B;
    raw(2:2:end,2:2:end) = handles.G;
    raw_size = size(raw);
    raw = raw';
    raw = uint16(raw(:));
    
    % Call the chart detector
	try
		[chart_found, chart_coord] = GMB_Detector(raw_size(2), raw_size(1), raw);
	catch ME
		disp('Error when calling GMB_Detector');
		chart_found = 0;
	end
    
    % Update metadata if the chart was found
    if chart_found > 0
        handles.metadata.gmb.tl = chart_coord(1,:)./[raw_size(2), raw_size(1)];
        handles.metadata.gmb.tr = chart_coord(2,:)./[raw_size(2), raw_size(1)];
        handles.metadata.gmb.br = chart_coord(3,:)./[raw_size(2), raw_size(1)];
        handles.metadata.gmb.bl = chart_coord(4,:)./[raw_size(2), raw_size(1)];
        handles.metadata.tag_ground_truth = 1;
        handles.gmb_chart_detected = 1;
    end
end

handles_out = handles;


% Manual addition: Add white points on the list based on the stored chart
% coordinates
function [handles_out] = wptool_add_grey_patches_from_gretag_macbeth(handles)

if handles.gmb_chart_detected > 0
    width  = size(handles.R_filt, 2);
    height = size(handles.R_filt, 1);
    gmb_coord = handles.metadata.gmb;
    patch_centers = wptool_get_gmb_patch_centers(gmb_coord, width, height);
    maxval = handles.ref_frame_format.satpoint;
    
    % Iterate through the grey patches
    for i = 19:24
        x = round(patch_centers(i, 1));
        y = round(patch_centers(i, 2));
        
        if x <= width && y <= height && x > 0 && y > 0
            % Get the current filtered RGB values
            currR = handles.R_filt(y, x);
            currG = handles.G_filt(y, x);
            currB = handles.B_filt(y, x);
            
            % Proceed if the data is not close to the ends of the dynamic range
            if (currR > 0.05*maxval && currR < 0.9*maxval && ...
                currG > 0.05*maxval && currG < 0.9*maxval && ...
                currB > 0.05*maxval && currB < 0.9*maxval)
                
                % Calculate the [R/G,B/G] of this position
                whitepoint_candidate = [currR./currG, currB./currG];
                
                % Add to the white point list
                new_wp_line = [whitepoint_candidate, 1];
                wp_data = get(handles.table_wps, 'Data');
                if numel(wp_data) < 3
                    wp_data = new_wp_line;
                else
                    wp_data = [wp_data; new_wp_line];
                end
                set(handles.table_wps, 'Data', wp_data);
            end % if data is valid
        end % if coordinate is valid
    end % for i
end

handles_out = handles;


% Manual addition: Write metadata to .meta file
function wptool_write_metadata(handles)

fname_meta = wptool_replace_file_extension(handles.filename, '.meta');

proceed = 1;

if exist(fname_meta, 'file') == 2
    reply = questdlg('The .meta file already exists. Overwrite?', ...
                     'The .meta file already exists', ...
                     'Yes','No','No');
    if strcmp(reply,'No') == 1
        proceed = 0;
    end
end % if

if proceed == 1
    fid = fopen(fname_meta,'w');
    
    fprintf(fid, '%s\t%f\n', 'exposure_time', handles.metadata.exposure_time);
    fprintf(fid, '%s\t%f\n', 'iso', handles.metadata.iso);
    fprintf(fid, '%s\t%f\n', 'analog_gain', handles.metadata.analog_gain);
    fprintf(fid, '%s\t%f\n', 'digital_gain', handles.metadata.digital_gain);
    fprintf(fid, '%s\t%f\n', 'aperture', handles.metadata.aperture);
    fprintf(fid, '%s\t%f\n', 'lux', handles.metadata.lux);
    fprintf(fid, '%s\t%f\n', 'lux_estimate', handles.metadata.lux_estimate);
    fprintf(fid, '%s\t%f\n', 'gmb_tl_x', handles.metadata.gmb.tl(1));
    fprintf(fid, '%s\t%f\n', 'gmb_tl_y', handles.metadata.gmb.tl(2));
    fprintf(fid, '%s\t%f\n', 'gmb_tr_x', handles.metadata.gmb.tr(1));
    fprintf(fid, '%s\t%f\n', 'gmb_tr_y', handles.metadata.gmb.tr(2));
    fprintf(fid, '%s\t%f\n', 'gmb_br_x', handles.metadata.gmb.br(1));
    fprintf(fid, '%s\t%f\n', 'gmb_br_y', handles.metadata.gmb.br(2));
    fprintf(fid, '%s\t%f\n', 'gmb_bl_x', handles.metadata.gmb.bl(1));
    fprintf(fid, '%s\t%f\n', 'gmb_bl_y', handles.metadata.gmb.bl(2));
    fprintf(fid, '%s\t%f\n', 'tag_ground_truth', handles.metadata.tag_ground_truth);
    fprintf(fid, '%s\t%f\n', 'tag_mixed_illumination', handles.metadata.tag_mixed_illumination);
    fprintf(fid, '%s\t%f\n', 'normalized_exposure_s', handles.metadata.normalized_exposure_s);
    
    fclose(fid);
    
    disp('.meta file written');
end % if


% Manual addition: Update the metadata listbox contents
function [handles_out] = wptool_update_metadata_listbox(handles)

handles.listbox_metadata.String = '';

listbox_data{1} = sprintf('exposure_time = %f', handles.metadata.exposure_time);
listbox_data{2} = sprintf('iso = %d', handles.metadata.iso);
listbox_data{3} = sprintf('analog_gain = %.2f', handles.metadata.analog_gain);
listbox_data{4} = sprintf('digital_gain = %.2f', handles.metadata.digital_gain);
listbox_data{5} = sprintf('aperture = %.2f', handles.metadata.aperture);
listbox_data{6} = sprintf('lux = %.2f', handles.metadata.lux);
listbox_data{7} = sprintf('lux_estimate = %.2f', handles.metadata.lux_estimate);
listbox_data{8} = sprintf('GMB = {(%.3f,%.3f), (%.3f,%.3f), (%.3f,%.3f), (%.3f,%.3f)}', ...
                  handles.metadata.gmb.tl(1), handles.metadata.gmb.tl(2), ...
                  handles.metadata.gmb.tr(1), handles.metadata.gmb.tr(2), ...
                  handles.metadata.gmb.br(1), handles.metadata.gmb.br(2), ...
                  handles.metadata.gmb.bl(1), handles.metadata.gmb.bl(2) ...
                  );
if handles.metadata.tag_ground_truth == 1
    tag_ground_truth_str = 'yes';
else
    tag_ground_truth_str = 'no';
end
listbox_data{9} = sprintf('tag_ground_truth = %s', tag_ground_truth_str);
if handles.metadata.tag_mixed_illumination == 1
    tag_mixed_illumination_str = 'yes';
else
    tag_mixed_illumination_str = 'no';
end
listbox_data{10} = sprintf('tag_mixed_illumination = %s', tag_mixed_illumination_str);
listbox_data{11} = sprintf('normalized_exposure_s = %f', handles.metadata.normalized_exposure_s);

set(handles.listbox_metadata, 'string', listbox_data);
set(handles.listbox_metadata, 'Min', 0, 'Max', 2, 'Value', [])
set(handles.listbox_metadata, 'BackgroundColor', [0.9, 0.9, 0.9]);
set(handles.listbox_metadata, 'enable', 'inactive');

handles_out = handles;


% Manual addition: Update the checkboxes in the tags panel
function [handles_out] = wptool_update_tags_checkboxes(handles)

if handles.metadata.tag_ground_truth > 0
    handles.gtcheckbox.Value = 1;
else
    handles.gtcheckbox.Value = 0;
end

if handles.metadata.tag_mixed_illumination > 0
    handles.mixedcheckbox.Value = 1;
else
    handles.mixedcheckbox.Value = 0;
end

handles_out = handles;

    
% Manual addition: Find WP from the struct by the name field value
function [wp] = wptool_find_wp_by_name(wps_ccms, name)

wp = [];
for i = 1:length(wps_ccms)
    if strcmp(wps_ccms(i).name, name)
        wp = wps_ccms(i).wp;
        break;
    end
end % for


% Manual addition: Find CCM from the struct by the name field value
function [ccm] = wptool_find_ccm_by_name(wps_ccms, name)

ccm = diag([1 1 1]);
for i = 1:length(wps_ccms)
    if strcmp(wps_ccms(i).name, name)
        ccm = wps_ccms(i).ccm;
        break;
    end
end % for

% Manual addition: Return an image that is suitable for illustration
function [out_img] = wptool_illustrate_image( R, G, B, wb_gains, ccm, satpoint, gamma )

% Apply WB
img_size = size(R);
tmp = zeros([img_size,3]);
tmp(:,:,1) = double(R).*wb_gains(1);
tmp(:,:,2) = double(G).*wb_gains(2);
tmp(:,:,3) = double(B).*wb_gains(3);
tmp(tmp>satpoint) = satpoint; % Apply clipping to avoid coloration of saturated areas

% Apply CCM
tmp = wptool_apply_ccm(tmp, ccm, satpoint);

% Apply brightness normalization
tmp = wptool_normalize_brightness(tmp, satpoint, 2.0, 0.99);

% Limit to range [0,1] and apply gamma
tmp = tmp./satpoint;
tmp = tmp.^gamma; % Approximation of the sRGB gamma

% Apply sharpening
out_img = imsharpen(tmp, 'Radius', 1.0, 'Amount', 1.3, 'Threshold', 0);


% Manual addition: Illustrate Gretag MacBeth chart patch locations that
% correspond to the given corner co-ordinates
function [img_out] = wptool_illustrate_gmb( img_in, gmb_coord )

% Input image dimensions
width  = size(img_in, 2);
height = size(img_in, 1);

% Get patch centers
patch_centers = wptool_get_gmb_patch_centers(gmb_coord, width, height);

% Draw the illustrations
img_in = insertShape(img_in, 'Circle', patch_centers(1:3*6,:), 'LineWidth', 4.0, 'Color', [1.0, 1.0, 1.0], 'SmoothEdges', true); % Color patches
img_in = insertShape(img_in, 'Circle', patch_centers(3*6+1:4*6,:), 'LineWidth', 6.0, 'Color', [0.6196, 0.7373, 0.2510], 'SmoothEdges', true); % Grey patches

img_out = img_in;


function [gmb_patch_centers] = wptool_get_gmb_patch_centers(gmb_coord, width, height)

tl = gmb_coord.tl; % top-left
tr = gmb_coord.tr; % top-right
br = gmb_coord.br; % bottom-right
bl = gmb_coord.bl; % bottom-left

% Color patch centers
gmb_patch_centers = zeros(4*6,3);
ind = 1;
for j = 1:4
    ls = tl + (bl-tl)*(j*2-1)/8; % chart line start point
    le = tr + (br-tr)*(j*2-1)/8; % chart line end point
    for i = 1:6
        p = ls + (le-ls)*(i*2-1)/12; % patch center point
        gmb_patch_centers(ind,:) = [p(1)*width, p(2)*height, norm(le-ls)/24*height];
        ind = ind + 1;
    end % for i
end % for j


% Manual addition: Normalize brightness of the input image
function [img_out] = wptool_normalize_brightness(img_in, satpoint, max_extra_gain, max_extra_gain_br_ind_def)

% Initializations
br_ind_tgt = 0.9;  % Where the bright end of histogram should be mapped on the dynamic range [0,1]

% Calculate cumulative histogram of all the pixel values
h = hist(img_in(:), 1000);
hcum = h;
for i = 2:1000,
    hcum(i) = hcum(i) + hcum(i-1);
end

% Find the bright index
total_cnt = hcum(1000);
br_ind = 1;
while br_ind < 1000 && hcum(br_ind) <= total_cnt*max_extra_gain_br_ind_def,
   br_ind = br_ind + 1; 
end
br_ind = max(1, br_ind - 1)/1000;

% Calculate and apply extra gain
gain = min(max_extra_gain, max(1.0, br_ind_tgt/br_ind));
img_out = img_in.*gain;
disp(sprintf('Brightness normalization gain %.3f applied', gain));

img_out(img_out>satpoint) = satpoint; % Apply clipping


% Manual addition: Apply 3x3 color conversion matrix
function [out_img] = wptool_apply_ccm(in_img, ccm, satpoint)
    img_size = size(in_img);
    R = in_img(:,:,1);
    G = in_img(:,:,2);
    B = in_img(:,:,3);
    tmp = zeros([3,img_size(1)*img_size(2)]);
    tmp(1,:) = R(:);
    tmp(2,:) = G(:);
    tmp(3,:) = B(:);
    tmp = ccm*tmp;
    tmp(tmp<0) = 0;
    tmp(tmp>satpoint) = satpoint;
    R = tmp(1,:);
    G = tmp(2,:);
    B = tmp(3,:);
    R = reshape(R,img_size(1),img_size(2));
    G = reshape(G,img_size(1),img_size(2));
    B = reshape(B,img_size(1),img_size(2));
    out_img = zeros([img_size(1), img_size(2), 3]);
    out_img(:,:,1) = R;
    out_img(:,:,2) = G;
    out_img(:,:,3) = B;


% --- Executes when selected object is changed in btgrp_ccms.
function btgrp_ccms_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in btgrp_ccms 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: It would be better if the button group would be populated by
% buttons in run-time, based on contents of handles.ref_wps_ccms.wps_ccms.
% Then e.g. the tag of the button would directly indicate the entry in the
% struct, and there would not be any other buttons than the ones actually
% covered by the struct.

ccm = handles.ccm;

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object
    case 'rb_sl_d_ccm'
        ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_D');
    case 'rb_sl_cw_ccm'
        ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_CW');
    case 'rb_sl_hor_ccm'
        ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_Hor');
    case 'rb_sl_a_ccm'
        ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_A');
    case 'rb_sl_tl84u30_ccm'
        ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_TL84U30');
    case 'rb_ie_d65_ccm'
        ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'IE_D65');
    case 'rb_ie_d50_ccm'
        ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'IE_D50');
    case 'rb_ie_f11_ccm'
        ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'IE_F11');
    case 'rb_ie_f12_ccm'
        ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'IE_F12');
    case 'rb_ie_a_ccm'
        ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'IE_A');
    case 'rb_current_ccm_file'
        % Check if .ccm file already exists for this raw file, and use it if available
        fname_ccm = wptool_replace_file_extension(handles.filename, '.ccm');
        if exist(fname_ccm, 'file') == 2
            ccm = load(fname_ccm);
        else
            set(handles.btgrp_ccms,'SelectedObject',handles.rb_sl_d_ccm);
            ccm = wptool_find_ccm_by_name(handles.ref_wps_ccms.wps_ccms, 'SL_D');
            %handles.text_infotext.String = '.ccm file is not available. Switching to SL_D.';
            msgbox('.ccm file is not available. Switching to SL_D.');
        end % if
end % switch

handles.ccm = ccm;

handles.text_ccm.String = num2str(reshape(handles.ccm,3,3)', 4);

% Update handles structure
guidata(hObject, handles);

% Manual addition: Convert [R/G,B/G] to WB gains
function [wb_gains] = wptool_convert_wp_to_gains(wp)

wb_gains = [1/wp(1), 1, 1/wp(2)];
mingain = min(wb_gains(:));
wb_gains = wb_gains./mingain;


% Manual addition: Called when mouse button is pressed over the fig_image
function wptool_get_mouse_position(hObject, ~)

% Get current mouse coordinates
handles = guidata(hObject);
cursorPoint = get(handles.fig_image, 'CurrentPoint');
curX = floor(cursorPoint(1,1));
curY = floor(cursorPoint(1,2));
tmpSize = size(handles.R_filt);
hf = get(handles.fig_image,'parent');
seltype = get(hf,'selectiontype');

if curX <= tmpSize(2) && curY <= tmpSize(1) && curX > 0 && curY > 0
    % Get the current filtered RGB values
    currR = handles.R_filt(curY, curX);
    currG = handles.G_filt(curY, curX);
    currB = handles.B_filt(curY, curX);
    
    % Update the RGB value in the UI
    set(handles.text_linrawrgb, 'String', strcat('Linearized raw RGB = [',num2str(currR,'%.0f'),', ',num2str(currG,'%.0f'),', ',num2str(currB,'%.0f'),']'));
    
    % Calculate the [R/G,B/G] of this position
    handles.whitepoint_candidate = [currR./currG, currB./currG];
    handles.text_cand_rperg.String = num2str(handles.whitepoint_candidate(1),4);
    handles.text_cand_bperg.String = num2str(handles.whitepoint_candidate(2),4);
    
    % Check if the second mouse button was clicked; if yes, then the white
    % point is added in list of marked white points for plotting purposes
    if strcmpi(seltype,'alt') % 'normal' would mean left mouse button, 'alt' means right mouse button
        tmp_len = length(handles.marked_whitepoints);
        handles.marked_whitepoints{tmp_len+1} = handles.whitepoint_candidate;
    end
    
    % Update the list of 4 most recent coordinates (used for the masking functionality)
    handles.previous_coordinates(handles.previous_coordinates_ind,:) = [curX, curY];
    handles.previous_coordinates_ind = mod(handles.previous_coordinates_ind,4)+1;
    
    guidata(hObject, handles);
    
    % Plot the [R/G,B/G] -data
    wptool_update_chrplot(handles);
    
    % Add automatically to the white point list if the check box is checked
    if handles.autoaddcheckbox.Value > 0 && handles.whitepoint_candidate(1) > 0 && handles.whitepoint_candidate(2) > 0
        new_wp_line = [handles.whitepoint_candidate, 1];
        wp_data = get(handles.table_wps, 'Data');
        if numel(wp_data) < 3
            wp_data = new_wp_line;
        else
            wp_data = [wp_data; new_wp_line];
        end
        set(handles.table_wps, 'Data', wp_data);
        
        guidata(hObject, handles);
    end % if
end


% Manual addition: Called when mouse button is pressed over the fig_chrplot
function wptool_get_mouse_position_chrplot(hObject, ~)

% Get current mouse coordinates
handles = guidata(hObject);
cursorPoint = get(handles.fig_chrplot, 'CurrentPoint');
curX = round(cursorPoint(1,1), 4);
curY = round(cursorPoint(1,2), 4);
xlimCur = xlim();
ylimCur = ylim();

if curX < xlimCur(2) && curY < ylimCur(2) && curX > xlimCur(1) && curY > ylimCur(1)
    % Calculate the [R/G,B/G] of this position
    handles.whitepoint_candidate = [curX, curY];
    handles.text_cand_rperg.String = num2str(handles.whitepoint_candidate(1),4);
    handles.text_cand_bperg.String = num2str(handles.whitepoint_candidate(2),4);
    
    guidata(hObject, handles);
    
    % Plot the [R/G,B/G] -data
    wptool_update_chrplot(handles)
end


% --- Executes on mouse press over axes background.
function fig_image_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to fig_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_add_wpcandidate.
function pushbutton_add_wpcandidate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_wpcandidate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.whitepoint_candidate(1) > 0 && handles.whitepoint_candidate(2) > 0
    new_wp_line = [handles.whitepoint_candidate, 1];
    wp_data = get(handles.table_wps, 'Data');
    if numel(wp_data) < 3
        wp_data = new_wp_line;
    else
        wp_data = [wp_data; new_wp_line];
    end
    set(handles.table_wps, 'Data', wp_data);
end % if


% --- Executes on mouse press over axes background.
function fig_chrplot_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to fig_chrplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Manual addition: Read .plain16 raw Bayer file
function [out_img] = wptool_plain16_read(fname, width, height)

f = fopen(fname);
a = dir(fname);
if ~exist(fname,'file') || length(a) ~= 1 || f < 0,
    error(['Cannot open file: ' fname]);
end
if a.bytes ~= width*height*2,
    fclose(f);
    error(['File size is wrong: ' fname])
end

out_img = fread(f, width*height, 'uint16');
fclose(f);
out_img = reshape(out_img, width, height)';


% Manual addition: Write .plain16 raw Bayer file
function [] = wptool_plain16_write(raw_data, out_fname)

f = fopen(out_fname, 'w+');
c = fwrite(f, raw_data', 'uint16'); % written column-wise, hence transpose
fclose(f);

disp(sprintf('%d elements written to file %s', c, out_fname));


% Manual addition: Replace file extension
function [ fname_out ] = wptool_replace_file_extension( fname_in, ext_new )

i = length(fname_in);
if i < 3
    error('Invalid file name');
end
while fname_in(i) ~= '.'
    if fname_in(i) == '\' || fname_in(i) == '/' || i < 3
        error('Invalid file name');
    end
    i = i - 1;
end
fname_out = [fname_in(1:i-1), ext_new];


% --- Executes on button press in pushbutton_recent_wps.
function pushbutton_recent_wps_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_recent_wps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if exist('wptool_recent_wps.wp', 'file') == 2
    wplist = load('wptool_recent_wps.wp');
    wplist_str = num2str(wplist);
    wplist_sz = size(wplist);
    wplist_cell = cell([wplist_sz(1),1]);
    for i=1:wplist_sz(1),
        wplist_cell{i} = wplist_str(i,:);
    end
    
    [ind,ok] = listdlg('PromptString', 'Select a WP:',...
                       'SelectionMode', 'single',...
                       'ListString', wplist_cell);
    
    if ok,
        handles.whitepoint_candidate = wplist(ind,:);
        handles.text_cand_rperg.String = num2str(handles.whitepoint_candidate(1),4);
        handles.text_cand_bperg.String = num2str(handles.whitepoint_candidate(2),4);
        
        guidata(hObject, handles);
        
        % Plot the [R/G,B/G] -data
        wptool_update_chrplot(handles)
    end % if
else
    msgbox('File wptool_recent_wps.wp not found');
end % if


% Manual addition: Separate the color components from raw Bayer image
function [Gr, R, B, Gb] = wptool_separateChannels(image, order)

if size(image, 3) ~= 1
    error('wptool: Input image must be a raw Bayer image');
end % if

switch order
    case 0, [Gr, R, B, Gb] = wptool_convertToChannels(image);
    case 1, [R, Gr, Gb, B] = wptool_convertToChannels(image);
    case 2, [B, Gb, Gr, R] = wptool_convertToChannels(image);
    case 3, [Gb, B, R, Gr] = wptool_convertToChannels(image);
end % switch


% Manual addition: For wptool_separateChannels
function [ch1, ch2, ch3, ch4] = wptool_convertToChannels(image)

ch1 = image(1:2:end, 1:2:end);
ch2 = image(1:2:end, 2:2:end);
ch3 = image(2:2:end, 1:2:end);
ch4 = image(2:2:end, 2:2:end);


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if min(handles.previous_coordinates(:)) > 0 % Check if coordinates have been clicked at least 4 times (the matrix is initialized with values -1)
    prompt = sprintf('Coordinates: [%d, %d; %d, %d; %d, %d; %d, %d]. Mask?', handles.previous_coordinates');
    reply = questdlg({prompt}, 'Masking tool', ...
                     'Yes','No','No');
    if strcmp(reply,'Yes') == 1
        % Check that the coordinates have been clicked clockwise in order such that
        % the vectors enclose the center area
        ps = handles.previous_coordinates;
        c = [mean(ps(:,1)), mean(ps(:,2))]; % center point
        ps = ps - [c;c;c;c]; % move the center point to the origo
        cp1 = ps(1,1)*ps(2,2) - ps(1,2)*ps(2,1);
        cp2 = ps(2,1)*ps(3,2) - ps(2,2)*ps(3,1);
        cp3 = ps(3,1)*ps(4,2) - ps(3,2)*ps(4,1);
        cp4 = ps(4,1)*ps(1,2) - ps(4,2)*ps(1,1);
        if min([cp1, cp2, cp3, cp4]) > 0 % If the points have been clicked clockwise
            ok = 1;
        else
            ok = 0;
            msgbox('The coordinates have not been clicked in clockwise order to enclose the center region');
        end
        ps = ps + [c;c;c;c]; % move the center point back
        
        % Calculate the original raw image coordinates and mask the area that is outside the clicked quadrangle
        if ok == 1
            ds_size = size(handles.R);
            ds_width_cr = ds_size(2);
            ds_height_cr = ds_size(1);
            
            % Rotate the coordinates if the original raw image is rotated
            if handles.ref_frame_format.rotate180 == 1
                for i=1:4
                    ps(i,:) = [ds_width_cr-ps(i,1)+1, ds_height_cr-ps(i,2)+1];
                end % for
            end % if
            
            org_width_cr = handles.ref_frame_format.width - handles.ref_frame_format.crop_left - handles.ref_frame_format.crop_right;
            org_height_cr = handles.ref_frame_format.height - handles.ref_frame_format.crop_top - handles.ref_frame_format.crop_bottom;
            ps_org = ps;
            for i=1:4
                ps_org(i,1) = floor(ps_org(i,1)*org_width_cr/ds_width_cr);
                ps_org(i,2) = floor(ps_org(i,2)*org_height_cr/ds_height_cr);
            end % for
            for i=1:4
                ps_org(i,1) = ps_org(i,1) + handles.ref_frame_format.crop_left;
                ps_org(i,2) = ps_org(i,2) + handles.ref_frame_format.crop_top;
            end % for
            xmin = min(ps_org(:,1)); % For speedup, calculate the outer bounding box [ymin:ymax,xmin:xmax]
            xmin = floor(xmin/2)*2+1;
            xmax = max(ps_org(:,1));
            xmax = floor(xmax/2)*2+1;
            ymin = min(ps_org(:,2));
            ymin = floor(ymin/2)*2+1;
            ymax = max(ps_org(:,2));
            ymax = floor(ymax/2)*2+1;
            tl_ind = 0; % For speedup, calculate the inner box, contained by the pointed quadrangle [ymin2:ymax2,xmin2:xmax2]
            min_sum = org_width_cr+org_height_cr;
            for i=1:4
                if ps_org(i,1)+ps_org(i,2) < min_sum
                    tl_ind = i;
                    min_sum = ps_org(i,1)+ps_org(i,2);
                end
            end % for
            tr_ind = mod(tl_ind,4)+1;
            br_ind = mod(tl_ind+1,4)+1;
            bl_ind = mod(tl_ind+2,4)+1;
            xmin2 = max(ps_org(tl_ind,1), ps_org(bl_ind,1));
            xmin2 = floor(xmin2/2)*2+1;
            xmax2 = min(ps_org(tr_ind,1), ps_org(br_ind,1));
            xmax2 = floor(xmax2/2)*2;
            ymin2 = max(ps_org(tl_ind,2), ps_org(tr_ind,2));
            ymin2 = floor(ymin2/2)*2+1;
            ymax2 = min(ps_org(bl_ind,2), ps_org(br_ind,2));
            ymax2 = floor(ymax2/2)*2;
            
            % Go through the raw image and apply the masking
            fname = handles.filename;
            width = handles.ref_frame_format.width;
            height = handles.ref_frame_format.height;
            raw_data = wptool_plain16_read(fname, width, height); % Read the data from .plain16 file
            if 1 % Rectangle masking
                tmp = zeros(size(raw_data));
                tmp(ymin2:ymax2,xmin2:xmax2) = raw_data(ymin2:ymax2,xmin2:xmax2);
                raw_data = tmp;
            else % Accurate masking
                raw_data(1:end,1:xmin-1) = 0;
                raw_data(1:end,xmax+1:end) = 0;
                raw_data(1:ymin-1,xmin:xmax) = 0;
                raw_data(ymax+1:end,xmin:xmax) = 0;
                for x = xmin:2:xmin2,
                    for y = ymin:2:ymax,
                        in = inpolygon(x,y,ps_org(:,1), ps_org(:,2));
                        if in == 0 % Mask the whole Bayer quad to value 0 if the coordinate is outside the indicated quadrangle
                            raw_data(y,x)     = 0;
                            raw_data(y,x+1)   = 0;
                            raw_data(y+1,x)   = 0;
                            raw_data(y+1,x+1) = 0;
                        end % if
                    end % for y
                end % for x
                for x = xmax2:2:xmax,
                    for y = ymin:2:ymax,
                        in = inpolygon(x,y,ps_org(:,1), ps_org(:,2));
                        if in == 0 % Mask the whole Bayer quad to value 0 if the coordinate is outside the indicated quadrangle
                            raw_data(y,x)     = 0;
                            raw_data(y,x+1)   = 0;
                            raw_data(y+1,x)   = 0;
                            raw_data(y+1,x+1) = 0;
                        end % if
                    end % for y
                end % for x
                for x = xmin2:2:xmax2,
                    for y = ymin:2:ymin2,
                        in = inpolygon(x,y,ps_org(:,1), ps_org(:,2));
                        if in == 0 % Mask the whole Bayer quad to value 0 if the coordinate is outside the indicated quadrangle
                            raw_data(y,x)     = 0;
                            raw_data(y,x+1)   = 0;
                            raw_data(y+1,x)   = 0;
                            raw_data(y+1,x+1) = 0;
                        end % if
                    end % for y
                end % for x
                for x = xmin2:2:xmax2,
                    for y = ymax2:2:ymax,
                        in = inpolygon(x,y,ps_org(:,1), ps_org(:,2));
                        if in == 0 % Mask the whole Bayer quad to value 0 if the coordinate is outside the indicated quadrangle
                            raw_data(y,x)     = 0;
                            raw_data(y,x+1)   = 0;
                            raw_data(y+1,x)   = 0;
                            raw_data(y+1,x+1) = 0;
                        end % if
                    end % for y
                end % for x
            end % if
            
            % Write the result back to file
            fname_bak = wptool_replace_file_extension(handles.filename, '.plain16.nomask');
            if exist(fname_bak, 'file') ~= 2 % Make a backup if not done already
                copyfile(handles.filename, fname_bak);
            end % if
            wptool_plain16_write(raw_data, handles.filename);
            
            % Re-initialize
            handles = wptool_init_plots(handles);
            handles.previous_coordinates = [-1,-1; -1,-1; -1,-1; -1,-1];
            handles.previous_coordinates_ind = 1;
            wptool_update_plots(handles);
            guidata(hObject, handles);
        end % if
    end % if
end % if


% --- Executes on button press in autoaddcheckbox.
function autoaddcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to autoaddcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoaddcheckbox


% --- Executes on button press in gtcheckbox.
function gtcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to gtcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gtcheckbox

gt_tag_enable = get(hObject,'Value');

if gt_tag_enable > 0
    handles.metadata.tag_ground_truth = 1;
else
    handles.metadata.tag_ground_truth = 0;
end

handles = wptool_update_metadata_listbox(handles);

% Update handles structure
guidata(hObject, handles);


function edit_metadatapath_Callback(hObject, eventdata, handles)
% hObject    handle to edit_metadatapath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_metadatapath as text
%        str2double(get(hObject,'String')) returns contents of edit_metadatapath as a double


% --- Executes during object creation, after setting all properties.
function edit_metadatapath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_metadatapath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

wptool_write_metadata(handles);


% --- Executes on button press in mixedcheckbox.
function mixedcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to mixedcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mixedcheckbox

mi_tag_enable = get(hObject,'Value');

if mi_tag_enable > 0
    handles.metadata.tag_mixed_illumination = 1;
else
    handles.metadata.tag_mixed_illumination = 0;
end

handles = wptool_update_metadata_listbox(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if min(handles.previous_coordinates(:)) > 0 % Check if coordinates have been clicked at least 4 times (the matrix is initialized with values -1)
    prompt = sprintf('Coordinates: [%d, %d; %d, %d; %d, %d; %d, %d]. Store Gretag MacBeth chart location?', handles.previous_coordinates');
    reply = questdlg({prompt}, 'Store Gretag MacBeth chart location', ...
                     'Yes','No','Yes');
    if strcmp(reply,'Yes') == 1
        % Check that the coordinates have been clicked clockwise in order such that
        % the vectors enclose the center area
        ps = handles.previous_coordinates;
        c = [mean(ps(:,1)), mean(ps(:,2))]; % center point
        ps = ps - [c;c;c;c]; % move the center point to the origo
        cp1 = ps(1,1)*ps(2,2) - ps(1,2)*ps(2,1);
        cp2 = ps(2,1)*ps(3,2) - ps(2,2)*ps(3,1);
        cp3 = ps(3,1)*ps(4,2) - ps(3,2)*ps(4,1);
        cp4 = ps(4,1)*ps(1,2) - ps(4,2)*ps(1,1);
        if min([cp1, cp2, cp3, cp4]) > 0 % If the points have been clicked clockwise
            ok = 1;
        else
            ok = 0;
            msgbox('The chart corner coordinates have not been clicked in clockwise order. Operation cancelled.');
        end
        ps = ps + [c;c;c;c]; % move the center point back
        
        % Store the Gretag MacBeth chart relative coordinates
        if ok == 1
            ds_size = size(handles.R_filt);
            ds_width = ds_size(2);
            ds_height = ds_size(1);
            
            handles.metadata.gmb.tl = [ps(1,1)/ds_width, ps(1,2)/ds_height];
            handles.metadata.gmb.tr = [ps(2,1)/ds_width, ps(2,2)/ds_height];
            handles.metadata.gmb.br = [ps(3,1)/ds_width, ps(3,2)/ds_height];
            handles.metadata.gmb.bl = [ps(4,1)/ds_width, ps(4,2)/ds_height];
            handles = wptool_update_metadata_listbox(handles);
            wptool_update_plots(handles);
            
            % Re-initialize
            handles.previous_coordinates = [-1,-1; -1,-1; -1,-1; -1,-1];
            handles.previous_coordinates_ind = 1;
            guidata(hObject, handles);
        end % if
    end % if
end % if

% --- Executes on button press in pushbutton12 (lux estimate to metadata)
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get characterization values
base_iso = handles.other_characterization.base_iso;
base_iso_aperture = handles.other_characterization.base_iso_aperture;
unity_gain_iso = handles.other_characterization.unity_gain_iso;

% Get current exposure values
exposure_time_s = handles.metadata.exposure_time;
iso = handles.metadata.iso;
analog_gain = handles.metadata.analog_gain;
digital_gain = max(1.0, handles.metadata.digital_gain);
aperture = handles.metadata.aperture;

% Set aperture gain according to current aperture vs. base ISO aperture
aperture_gain = 1.0;
if aperture > 0 && base_iso_aperture > 0 && aperture ~= base_iso_aperture
    aperture_gain = (base_iso_aperture^2)/(aperture^2);
end

% Calculate sensitivity normalized exposure
ok = 1;
if exposure_time_s > 0 && iso > 0 && base_iso > 0 % The DSLRs have the ISO value in the metadata
    
    norm_tot_exposure_ms = aperture_gain*exposure_time_s*1000*iso/unity_gain_iso*base_iso/100;
    
    disp(sprintf('norm_tot_exposure_ms = %.3f, aperture_gain = x%.3f, exposure_time = %.6fs, iso = %d, base_iso = %d', norm_tot_exposure_ms, aperture_gain, exposure_time_s, iso, base_iso));
    
elseif exposure_time_s > 0 && analog_gain > 0 && digital_gain > 0 && base_iso > 0 % The mobile camera has the gains in the metadata
    
    norm_tot_exposure_ms = aperture_gain*exposure_time_s*1000*analog_gain*digital_gain*base_iso/100;
    
    disp(sprintf('norm_tot_exposure_ms = %.3f, aperture_gain = x%.3f, exposure_time = %.6fs, analog_gain = x%.3f, digital_gain = x%.3f, base_iso = %d', norm_tot_exposure_ms, aperture_gain, exposure_time_s, analog_gain, digital_gain, base_iso));
    
else
    ok = 0;
end

% Check bpp
bpp = handles.ref_frame_format.bpp;
if bpp <= 0
    ok = 0;
end

% Check if the GMB location is stored in the metadata
if min(structfun(@(x)min(x(:)),handles.metadata.gmb)) > 0 && ok > 0

    width = size(handles.R_filt,2);
    height = size(handles.R_filt,1);
    
    % Get patch centers (in terms of handles.X_filt size)
    patch_centers = round(wptool_get_gmb_patch_centers(handles.metadata.gmb, width, height));
    
    % Relative brightnesses of the patches #20..#23 compared to patch #21
    gmb_patch_ratio_20_23_to_21 = [1.582, 1.0, 0.525, 0.266];
    
    % Set lux_mul that is the multiplier to convert the normalized digital
    % count per ms of exposure to lux. The value has been determined
    % empirically from the raw data, cf. the Sensitivities.xlsx.
    lux_mul = 45;
    
    % Iterate through grey patches #20..#23 (omit grey patches #19 and #24)
    lux_estimate = 0;
    for i = 20:1:23
        R = handles.R_filt(patch_centers(i,2), patch_centers(i,1));
        G = handles.G_filt(patch_centers(i,2), patch_centers(i,1));
        B = handles.B_filt(patch_centers(i,2), patch_centers(i,1));
        Y = max([R,G,B]);
        Y = Y/2^(bpp-10);
        
        lux_estimate = lux_estimate + ((Y/gmb_patch_ratio_20_23_to_21(i-19))/norm_tot_exposure_ms)*lux_mul;
    end % for i
    lux_estimate = lux_estimate/4;
    
    % Round the lux estimate
    lux_estimate = round_lux_est(lux_estimate);

    % Update the metadata listbox and handles
    handles.metadata.lux_estimate = lux_estimate;
    handles = wptool_update_metadata_listbox(handles);
    guidata(hObject, handles);
    
    disp(sprintf('lux_estimate = %d', lux_estimate));
    
else
    msgbox('Cannot calculate the lux estimate due to missing data');
end % if


% --- Executes on button press in pushbutton13 (read exposure parameters
% from EXIF)
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current metadata path and write it to wptool_recent_metadata_path.txt
data = get(handles.edit_metadatapath, 'String');
meta_path = data{1};
if ~strcmp('\', meta_path(end))
    meta_path = strcat(meta_path, '\');
end
fid = fopen('wptool_recent_metadata_path.txt','w');
fprintf(fid, '%s', meta_path);
fclose(fid);

% Construct the EXIF JPEG file name from the metadata path and the current
% raw file path and name
raw_path_and_name = handles.filename;
backslash_inds = strfind(raw_path_and_name, '\');
backslash_ind = 0;
if ~isempty(backslash_inds)
    backslash_ind = backslash_inds(end);
end
dot_inds = strfind(raw_path_and_name, '.');
dot_ind = length(raw_path_and_name)+1;
if ~isempty(dot_inds) % Find the index of the next dot from the last backslash
    dot_ind = backslash_ind+1;
    while min(dot_ind ~= dot_inds) > 0
        dot_ind = dot_ind + 1;
    end
end
search_name = raw_path_and_name(backslash_ind+1:dot_ind-1);
search_name = strrep(search_name, '_blc_lsc', '');

% Get the metadata from .i3av4 header in case of IMX135/Semco and .jpg EXIF otherwise
if strcmp('SonyIMX135', handles.cameraid) || strcmp('SonyIMX135_BLCCSC', handles.cameraid)
    selection_wildcard = strcat(search_name,'*.i3av4');
    fprintf('Looking for MKN in %s\n', [meta_path, selection_wildcard]);
    mkn_file_list = dir([meta_path, selection_wildcard]);
    if ~isempty(mkn_file_list)
        mkn_file = [meta_path, mkn_file_list(1).name];
        fprintf('Selected MKN file: %s\n', mkn_file);
        
        try
            mkn = parse_makernote(mkn_file); % parse_makernote needs to be in Matlab path
            exposure_time = double(mkn.E3A_DNID_AE_Result_Exposure_Time.data)/1000000.0;
            analog_gain = double(mkn.E3A_DNID_AE_Result_Analog_Gain.data);
            fprintf('Exposure time = %fs, Analog gain = %f\n', exposure_time, analog_gain);
            
            base_iso = handles.other_characterization.base_iso;
            normalized_exposure_s = calculate_normalized_exposure(exposure_time, analog_gain, -1.0, base_iso, -1.0);
            fprintf('Normalized exposure = %fs\n', normalized_exposure_s);
            
            % Update the metadata in the handles structure
            handles.metadata.exposure_time = exposure_time;
            handles.metadata.analog_gain = analog_gain;
            handles.metadata.normalized_exposure_s = normalized_exposure_s;
            
            % Refresh metadata listbox and update handles structure
            handles = wptool_update_metadata_listbox(handles);
            guidata(hObject, handles);
        catch ME
            msgbox(sprintf('Error reading the MakerNote: %s', ME.identifier));
        end % try-catch
    else
        msgbox('MKN file not found');
    end
else
    selection_wildcard = strcat(search_name,'*.jpg');
    fprintf('Looking for EXIF in %s\n', [meta_path, selection_wildcard]);
    exif_file_list = dir([meta_path, selection_wildcard]);
    if ~isempty(exif_file_list)
        exif_file = [meta_path, exif_file_list(1).name];
        fprintf('Selected EXIF file: %s\n', exif_file);
        
        exif = imfinfo(exif_file);
        exposure_time = exif.DigitalCamera.ExposureTime;
        iso = exif.DigitalCamera.ISOSpeedRatings;
        aperture = exif.DigitalCamera.FNumber;
        fprintf('Exposure time = %fs, ISO = %d, F# = %f\n', exposure_time, iso, aperture);
        
        base_iso = handles.other_characterization.base_iso;
        base_iso_aperture = handles.other_characterization.base_iso_aperture;
        unity_gain_iso = handles.other_characterization.unity_gain_iso;
        normalized_exposure_s = calculate_normalized_exposure(exposure_time, iso/unity_gain_iso, aperture, base_iso, base_iso_aperture);
        fprintf('Normalized exposure = %fs\n', normalized_exposure_s);
        
        % Update the metadata in the handles structure
        handles.metadata.exposure_time = exposure_time;
        handles.metadata.iso = iso;
        handles.metadata.aperture = aperture;
        handles.metadata.normalized_exposure_s = normalized_exposure_s;
        
        % Refresh metadata listbox and update handles structure
        handles = wptool_update_metadata_listbox(handles);
        guidata(hObject, handles);
    else
        msgbox('EXIF file not found');
    end
end


% --- Executes on selection change in listbox_metadata.
function listbox_metadata_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_metadata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_metadata contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_metadata


% --- Executes during object creation, after setting all properties.
function listbox_metadata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_metadata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [normalized_exposure] = calculate_normalized_exposure(exposure_time_s, gain, aperture, base_iso, base_iso_aperture)

% Set aperture gain according to current aperture vs. base ISO aperture
aperture_gain = 1.0;
if aperture > 0 && base_iso_aperture > 0 && aperture ~= base_iso_aperture
    aperture_gain = (base_iso_aperture^2)/(aperture^2);
end

normalized_exposure = exposure_time_s * gain * (base_iso/100) * aperture_gain;

