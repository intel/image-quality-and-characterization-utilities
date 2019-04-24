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

function [] = duplicate_metadata(dir_to_be_processed, selection_wildcard, duplicated_file_extension)
%DUPLICATE_METADATA Duplicate .wp files or .ccm files if some .plain16 files are
%lacking the corresponding annotation file (to be used if all the files have
%been captured under the same illumination)

filelist = dir([dir_to_be_processed, '**\', selection_wildcard]); % include subfolders
fname_previous_wp = '';
copycount = 0;

for i = 1:length(filelist),
    fname = [filelist(i).folder, '\', filelist(i).name]; % The current file name
    fname_wp = dm_replace_file_extension(fname, duplicated_file_extension);
    
    if exist(fname_wp, 'file') == 2
        fname_previous_wp = fname_wp; % Mark this file as the most recent white point file
    elseif exist(fname_previous_wp, 'file') == 2
        copyfile(fname_previous_wp, fname_wp); % Copy the .wp from previous annotated white point file
        copycount = copycount + 1;
    end % if
end % for

fprintf('%d %s files copied\n', copycount, duplicated_file_extension);

end % function


function [ fname_out ] = dm_replace_file_extension( fname_in, ext_new )

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

end % function
