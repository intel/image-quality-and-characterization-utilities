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

function [] = image_matrix(folder, selection_wildcard)
%IMAGE_MATRIX Create a matrix image out of multiple JPEG images

filelist = dir([folder, selection_wildcard]);
fcnt = length(filelist);

% Determine the matrix size
fname = [folder, filelist(1).name];
img = imread(fname, 'jpeg');
tn_height = 240; % Each thumbnail is 240 lines high, width is set by the aspect ratio
img = imresize(img, [tn_height nan], 'bilinear');
img_size = size(img);
matrix_wperh = (2/3)*3047.0/1269.0; % The aspect ratio of the image matrix
matrix_height = round(sqrt(fcnt/matrix_wperh));
matrix_width = ceil(fcnt/matrix_height);
img_matrix_height = matrix_height*img_size(1);
img_matrix_width  = matrix_width*img_size(2);

img_matrix = zeros(img_matrix_height, img_matrix_width,3);

fprintf('Creating image matrix of size %d x %d, containing %d x %d slots, for %d images in total\n', img_matrix_width, img_matrix_height, matrix_width, matrix_height, fcnt);

for y = 1:matrix_height,
    for x = 1:matrix_width,
        ind = (y-1)*matrix_width + x;
        if ind > fcnt
            break; % The last row might be only partially filled
        end % if
        fname = [folder, filelist(ind).name]; % The current image file name
        disp(sprintf('Processing file %s, index %d [%d,%d]', fname, ind, x, y));
        img = imread(fname, 'jpeg'); % Load the image
        img = double(img)./double(max(img(:)));
        img = imresize(img, img_size(1:2), 'bilinear'); % Downscale the image
        start_y = (y-1)*img_size(1) + 1;
        end_y   = start_y + img_size(1) - 1;
        start_x = (x-1)*img_size(2) + 1;
        end_x   = start_x + img_size(2) - 1;
        img_matrix(start_y:end_y,start_x:end_x,:) = img;
    end % for x
end % for y

fname_out = [folder, 'image_matrix.jpg'];

imwrite(img_matrix, fname_out, 'jpeg', 'Quality', 80);

disp('Done');

end % function

