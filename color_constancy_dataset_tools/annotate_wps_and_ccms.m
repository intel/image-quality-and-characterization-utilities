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

function [] = annotate_wps_and_ccms()
%ANNOTATE_WPS_AND_CCMS For calling the annotation tool
%   The main function for calling the white point and CCM annotation tool.
%   Use:
%   1) Enable the correct if-elseif-else branch according to the camera
%   2) Set the folder that contains the raw images (dir_to_be_processed)
%   3) Set which raw images should be processed (selection_wildcard)
%   4) Call annotate_wps_and_ccms()

% INPUTS -->

first_file_to_process = ''; % Give the name of the first file to process if some files from the beginning should be skipped
white_map_file = '';
if 0, % Sony IMX135 Semco CCMs / Without BLC or CSC
    dir_to_be_processed = '..\example_path\';
    selection_wildcard = '*.plain16';
    cameraid = 'SonyIMX135';
elseif 0, % Sony IMX135 Semco / With BLC and CSC
    dir_to_be_processed = '..\example_path\';
    selection_wildcard = '*.plain16';
    cameraid = 'SonyIMX135_BLCCSC';
elseif 0, % Canon EOS 5DSR
    dir_to_be_processed = '..\example_path\';
    selection_wildcard = '*.plain16';
    cameraid = 'Canon5DSR';
elseif 0, % Canon EOS 5Dmk3
    dir_to_be_processed = '..\example_path\';
    selection_wildcard = '*.plain16';
    cameraid = 'Canon5Dmk3';
else % Nikon D810
    dir_to_be_processed = '..\example_path\';
    selection_wildcard = '*.plain16';
    cameraid = 'NikonD810';
end % if

% <-- INPUTS

filelist = dir([dir_to_be_processed, '**\', selection_wildcard]); % include subfolders

istart = 1;
if strcmp(first_file_to_process,'')==0 % Check if some first files should be skipped (e.g. already annotated)
    for i = 1:length(filelist)
        if strcmp(filelist(i).name,first_file_to_process)==1
            istart = i;
            break;
        end
    end
    disp(sprintf('Skipping the first %d images (first_file_to_process == %s)', istart-1, first_file_to_process));
end

for i = istart:length(filelist),
    fname = [filelist(i).folder, '\', filelist(i).name]; % The current file name

    h = wptool(cameraid, fname, white_map_file); % Open the annotation tool

    waitfor(h); % Wait for the dialog to close before proceeding to the next file
end % for

% Duplicate white points and color conversion information for consecutive raw files that were skipped during annotation
duplicate_metadata(dir_to_be_processed, selection_wildcard, '.ccm');
duplicate_metadata(dir_to_be_processed, selection_wildcard, '.wp');

end % function


