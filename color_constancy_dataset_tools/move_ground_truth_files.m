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

function [] = move_ground_truth_files()
%MOVE_GROUND_TRUTH_FILES Move groud truth files to target folder
%   - Move all the files whose .meta file contains line
%     "tag_ground_truth	1.000000", i.e. the ground truth tag has value 1.
%   - The script can also be used to move any set of files that have e.g.
%     .tag or other similar file with the same name.

% INPUTS -->

dir_to_be_processed = '..\example_path\';
gt_tag_wildcard = '*.meta';
target_folder = '..\example_target_path\';

% <-- INPUTS

filelist = dir([dir_to_be_processed, gt_tag_wildcard]);

gt_inds = [1:length(filelist)];

% Remove those .meta files from filelist that do not contain
% tag_ground_truth with value 1
if strcmp(gt_tag_wildcard,'*.meta')==0
    for i = 1:length(filelist)
        fname_meta = [dir_to_be_processed, filelist(i).name]; % The current file name
        is_ok = false;
        if exist(fname_meta, 'file') == 2
            fid = fopen(fname_meta,'r');
            meta_struct = textscan(fid,'%s%f');
            fclose(fid);
            
            % Check if one of the lines in the .meta file contains the ground truth tag
            for j = 1:length(meta_struct{1})
                if strcmp(meta_struct{1}(j),'tag_ground_truth') && meta_struct{2}(j)>0
                    is_ok = true;
                end
            end % for j
        end % if
        
        % Remove file index from gt_inds if suitable tag was not found
        if ~is_ok
            gt_inds(gt_inds==i) = [];
        end
    end % for i
end % if

fprintf('Moving files to folder %s, %d tag files found\n', target_folder, length(gt_inds));
for i = gt_inds
    fname = [dir_to_be_processed, filelist(i).name]; % The current file name
    fname_base = replace_file_extension(fname, '');
    files_to_move = [fname_base, '*'];
    fprintf('Moving files %s\n', files_to_move);
    movefile(files_to_move, target_folder);
end % for

end % function

