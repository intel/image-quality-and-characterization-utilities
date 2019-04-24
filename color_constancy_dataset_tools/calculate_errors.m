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

function [] = calculate_errors()
%CALCULATE_ERRORS Calculate angular errors based on white point files

% INPUTS -->

dir_to_be_processed = '..\example_path\';

process_subdirs = 0;

ground_truth_file_extension = '.wp';

evaluated_file_extension = '.wp_example_awb';

max_file_cnt = 10000;

wpformat = 'WP_RperG_BperG'; % 'WP_RperG_BperG' or 'WP_RGB'

% <-- INPUTS

file_ind = 0;
name_list = cell(max_file_cnt,1);
err_list = zeros(max_file_cnt,1);

[name_list, err_list, file_ind] = process_directory(dir_to_be_processed, ground_truth_file_extension, evaluated_file_extension, process_subdirs, name_list, err_list, file_ind, wpformat);

file_id = fopen('angular_errors.txt','w');
fprintf(file_id,'%s\t%s\n','File name','Angular error');
for i=1:file_ind,
    fprintf(file_id,'%s\t%.4f\n',name_list{i},err_list(i));
end % for i
fclose(file_id);

end % function calculate_errors

function [name_list, err_list, file_ind] = process_directory(input_dir, ground_truth_file_extension, evaluated_file_extension, process_subdirs, name_list, err_list, file_ind, wpformat)

filelist = dir([input_dir, ['*',ground_truth_file_extension]]);

for i = 1:length(filelist),
    % Load white points
    fname = [input_dir, filelist(i).name];
    gt_wp = load(fname); % Read the ground truth white point from file
    fname_eval = replace_file_extension_nondot(fname, ground_truth_file_extension, evaluated_file_extension);
    %fname_eval = replace_file_extension(fname, evaluated_file_extension);
    eval_wp = load(fname_eval); % Read the evaluated white point from file
    % TODO: Handle the error cases in which the evaluated white point file does not exist

    % Calculate angular error
    switch wpformat
        case 'WP_RperG_BperG'
            gt_wpv = [gt_wp(1), 1, gt_wp(2)]; % Convert [R/G, B/G] to RGB vector
            eval_wpv = [eval_wp(1), 1, eval_wp(2)]; % Convert [R/G, B/G] to RGB vector
        case 'WP_RGB'
            gt_wpv = gt_wp;
            eval_wpv = eval_wp;
        otherwise
            error('Invalid wpformat selection');
    end % swtich
    ang_err = acosd(dot(gt_wpv,eval_wpv)/(norm(gt_wpv)*norm(eval_wpv))); % Angular error between ground truth and evaluted white point vectors

    % Update the results list
    file_ind = file_ind + 1;
    fname_tmp = replace_file_extension_nondot(fname, ground_truth_file_extension, '');
    name_list{file_ind} = fname_tmp;
    err_list(file_ind) = ang_err;

end % for i

if process_subdirs
    filelist = dir([input_dir, '*']);
    input_dir_org = input_dir;

    % Recursively call process_directory() for each subdirectory
    for i = 1:length(filelist),
        if filelist(i).isdir && ~strcmp(filelist(i).name,'.') && ~strcmp(filelist(i).name,'..')
            input_dir = [input_dir_org, filelist(i).name, '\'];
            [name_list, err_list, file_ind] = process_directory(input_dir, ground_truth_file_extension, evaluated_file_extension, process_subdirs, name_list, err_list, file_ind);
        end % if isdir
    end % for i
end % if process_subdirs

end % function process_directory
