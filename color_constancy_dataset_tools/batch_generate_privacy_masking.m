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

% INPUTS -->

last_file_to_skip = ''; % Set to empty string if no files are skipped

root_raw_folders = { ...
	'..\example_raw_folder_1\'; ...
	'..\example_raw_folder_2\'  ...
    };

frame_formats = { ...
	'frame_format_xxx.mat'; ... % Set according to camera type
    'frame_format_xxx.mat'; ... % Set according to camera type
    };

selection_wildcard = '*.plain16';

% <-- INPUTS

disp('Generating privacy masking for the raw files in the given folders...');

already_skipped = false;

for i = 1:length(root_raw_folders)
    dir_to_be_processed = root_raw_folders{i};
    fprintf('Processing folder %s\n', dir_to_be_processed);

    filelist = dir([dir_to_be_processed, '**\', selection_wildcard]); % include subfolders
    frame_format = load(frame_formats{i});
    width  = frame_format.width;
    height = frame_format.height;
    satpoint = frame_format.satpoint;
    pedestal = frame_format.pedestal;

    jstart = 1;
    if strcmp(last_file_to_skip,'')==0 && ~already_skipped % Check if some first files should be skipped (e.g. already processed earlier)
        for skip_ind = 1:length(filelist)
            if strcmp(filelist(skip_ind).name,last_file_to_skip)==1
                jstart = skip_ind+1;
                already_skipped = true;
                break;
            end
        end
        fprintf('Skipping the first %d images (last_file_to_skip == %s)\n', jstart-1, last_file_to_skip);
    end

    for j = jstart:length(filelist)
        folder_raw = filelist(j).folder;
        name_raw = filelist(j).name;
        fname_raw = [folder_raw, '\', name_raw];
        name_newraw = replace_file_extension(name_raw, '_privacy.plain16');
        fname_newraw = [folder_raw, '\', name_newraw];

        try
            % Read the raw Bayer data
            raw = plain16_read(fname_raw, width, height, 1);

            % Use the imshow figure as the UI
            fig_i = fig_CreateFigure(raw, fname_newraw, satpoint, pedestal);
            waitfor(fig_i); % Wait for the dialog to close before proceeding to the next file

        catch ME
            fprintf('Error reading raw file %s (%s)\n', fname_raw, ME.identifier);
        end % try
    end % for j

    fprintf('Folder %s processed\n', root_raw_folders{i});

end % for i

% Creation function for the figure
function [fig_i] = fig_CreateFigure(raw, fname_newraw, satpoint, pedestal)

fig_i = imshow(raw, [0, satpoint]);
figdata = guidata(fig_i);
figdata.raw = raw;
figdata.fname_newraw = fname_newraw;
figdata.curInd = 1;
figdata.coords = [-1, -1; -1, -1];
fig_i.ButtonDownFcn = @fig_ButtonDownFcn;

guidata(fig_i, figdata);

fig_img_data = get(fig_i, 'CData');
min_val = pedestal;
max_val = max(fig_img_data(:));
fig_img_data = max(fig_img_data - min_val, 0);
fig_img_data = fig_img_data.*(satpoint/(max_val-min_val));
fig_img_data = round(((fig_img_data./satpoint).^0.45).*satpoint);
set(fig_i, 'CData', fig_img_data);

end % function

% Callback function for mouse button press. Left mouse button marks
% coordinates and right mouse button closes the figure.
function fig_ButtonDownFcn(hObject, ~)

% Get axis handle and mouse button information
fig_axis = hObject.Parent;
fig_figure = fig_axis.Parent;
sel_type = get(fig_figure,'selectiontype');

% Check which mouse button was pressed: 'normal' means left mouse button, 'alt' means right mouse button
if strcmpi(sel_type,'normal')
    % Get current mouse coordinates
    cursorPoint = get(fig_axis, 'CurrentPoint');
    curX = floor(cursorPoint(1,1));
    curY = floor(cursorPoint(1,2));

    % Update the stored coordinate pairs
    figdata = guidata(hObject);
    figdata.coords(figdata.curInd, :) = [curX, curY];
    figdata.curInd = mod(figdata.curInd,2)+1;
    guidata(hObject, figdata);

    % If second coordinate has been clicked, mask the indicated rectangular area
    if figdata.curInd == 1 % TODO
        top    = floor(min(figdata.coords(:,2))/2)*2 + 1; % force odd number
        bottom = round(max(figdata.coords(:,2))/2)*2;     % force even number
        left   = floor(min(figdata.coords(:,1))/2)*2 + 1; % force odd number
        right  = round(max(figdata.coords(:,1))/2)*2;     % force even number
        raw_roi = figdata.raw(top:bottom, left:right);
        masked_roi = apply_mask_on_roi(raw_roi);

        % Apply the mask in the raw data
        figdata.raw(top:bottom, left:right) = masked_roi;
        guidata(hObject, figdata);

        % Apply the mask in the displayed data
        fig_img_data = get(hObject, 'CData');
        fig_img_data(top:bottom, left:right) = masked_roi;
        set(hObject, 'CData', fig_img_data);
    end

    fprintf('Button down [%d,%d] ([%d,%d] and [%d,%d])\n', curX, curY, figdata.coords(1,:), figdata.coords(2,:));
else % mouse button differentiation

    % New filename
    figdata = guidata(hObject);
    fname_newraw = figdata.fname_newraw;

    % Ask if a new file should be created
    reply = questdlg(sprintf('Write new file %s?', fname_newraw), ...
                     'Write new raw file', ...
                     'Yes','No','Yes');
    if strcmp(reply,'Yes') == 1
        try
            plain16_write(figdata.raw, fname_newraw);
        catch ME
            fprintf('Error writing raw file %s (%s)\n', fname_newraw, ME.identifier);
        end % try
    end

    % Close the figure
    fig_f = fig_axis.Parent;
    fig_f.WindowStyle = 'normal';
    close(fig_f);
end % mouse button differentiation

end % function

% Apply mask on the given raw Bayer ROI. Use average of each channel as the
% masking value in order not to introduce new colors outside of the current
% color gamut.
function [out_roi] = apply_mask_on_roi(in_roi)

ch1_avg = round(mean(mean(in_roi(1:2:end-1, 1:2:end-1))));
ch2_avg = round(mean(mean(in_roi(1:2:end-1, 2:2:end))));
ch3_avg = round(mean(mean(in_roi(2:2:end, 1:2:end-1))));
ch4_avg = round(mean(mean(in_roi(2:2:end, 2:2:end))));
fprintf('ROI averages = %d, %d, %d, %d\n', ch1_avg, ch2_avg, ch3_avg, ch4_avg);

out_roi = ones(size(in_roi));
out_roi(1:2:end-1, 1:2:end-1) = ch1_avg;
out_roi(1:2:end-1, 2:2:end)   = ch2_avg;
out_roi(2:2:end, 1:2:end-1)   = ch3_avg;
out_roi(2:2:end, 2:2:end)     = ch4_avg;

end % function

% Write a .plain16 raw Bayer file
function [] = plain16_write(raw_data, fname)

f = fopen(fname, 'w+');
c = fwrite(f, raw_data', 'uint16'); % written column-wise, hence transpose
fclose(f);

fprintf('%d elements written to file %s\n', c, fname);

end % function

