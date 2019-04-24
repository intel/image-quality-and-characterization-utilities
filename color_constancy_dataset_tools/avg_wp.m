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

function [ RperG, BperG ] = avg_wp( input_dir, selection_wildcard )
%AVG_WP Calculate the average white point within the given folder and files

filelist = dir([input_dir, '**\', selection_wildcard]); % include subfolders

filecnt = length(filelist);
wplist = zeros(filecnt,2);
disp( sprintf('Number of files: %d', filecnt) );

if filecnt < 1
    disp('No files found');
    RperG = 0;
    BperG = 0;
    return;
end

for i = 1:filecnt,
    fname = [filelist(i).folder, '\', filelist(i).name];

    wp = load(fname); % Read the white point from file
    wplist(i,:) = wp;

end % for i

if filecnt > 1
    wpmean = mean(wplist);
else
    wpmean = wplist;
end

fprintf('RperG: mean = %.3f, min = %.3f, max = %.3f\n', wpmean(1), min(wplist(:,1)), max(wplist(:,1)));
fprintf('BperG: mean = %.3f, min = %.3f, max = %.3f\n', wpmean(2), min(wplist(:,2)), max(wplist(:,2)));

RperG = wpmean(1);
BperG = wpmean(2);

end % function

