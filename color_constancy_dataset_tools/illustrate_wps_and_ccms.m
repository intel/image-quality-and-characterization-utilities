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

function [] = illustrate_wps_and_ccms()
%ILLUSTRATE_WPS_AND_CCMS For illustrating annotation results
%   Use:
%   1) Enable the correct if-elseif-else branch according to the camera
%   2) Set the folder that contains the raw images (dir_to_be_processed)
%   3) Set which raw images should be processed (selection_wildcard)
%   4) Call illustrate_wps_and_ccms()

% INPUTS -->

first_file_to_process = '';
default_ccm_file = '';
show_figure = 0;
wp_file_ext = '.wp';
out_file_ext = '.jpg';
if 0, % Sony IMX135 Semco / With BLC and CSC
	dir_to_be_processed = '..\example_path\';
    selection_wildcard = '*.plain16';
    
    pedestal = 0; % BLC and CSC applied already
    width = 3264;
    height = 2448;
    crop_top = 0;
    crop_left = 0;
    bayer_order = 0;
    bpp = 10;
    satpoint = 1023;
    rotate_180 = 0;
    max_extra_gain = 1.75;
    max_extra_gain_br_ind_def = 0.993;
elseif 0, % Sony IMX135 Semco / Without BLC and CSC
	dir_to_be_processed = '..\example_path\';
    selection_wildcard = '*.plain16';
    
    pedestal = 64;
    width = 3264;
    height = 2448;
    crop_top = 0;
    crop_left = 0;
    bayer_order = 0;
    bpp = 10;
    satpoint = 1023;
    rotate_180 = 0;%1;
    max_extra_gain = 1.75;
    max_extra_gain_br_ind_def = 0.993;
elseif 1, % Canon EOS 5DSR
    dir_to_be_processed = '..\example_path\';
    selection_wildcard = '*.plain16';
    
    pedestal = 2047;
    width = 8896;
    height = 5920;
    crop_top = 64;
    crop_left = 160;
    bayer_order = 1;
    bpp = 14;
    satpoint = 15380;
    rotate_180 = 0;
    max_extra_gain = 2.75;
    max_extra_gain_br_ind_def = 0.993;
else % Nikon D810
    dir_to_be_processed = '..\example_path\';
    selection_wildcard = '*.plain16';
    
    pedestal = 601;
    width = 7380;
    height = 4928;
    crop_top = 0;
    crop_left = 0;
    bayer_order = 1;
    bpp = 14;
    satpoint = 16383;
    rotate_180 = 0;
    max_extra_gain = 2.00;
    max_extra_gain_br_ind_def = 0.99;
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
    fprintf('Skipping the first %d images (first_file_to_process == %s)\n', istart-1, first_file_to_process);
end

iend = length(filelist);
for i = istart:iend
    fname = [filelist(i).folder, '\', filelist(i).name];
    fprintf('Image %d/%d: %s\n', i, iend, fname);
    data = plain16_read(fname, width, height);
    data = data(crop_top+1:end, crop_left+1:end); % Crop out light shielded pixels etc. (if any)
    [Gr, R, B, Gb] = iwc_separateChannels(data, bayer_order);
    G = (double(Gr)+double(Gb))./2;
    R = double(R);
    B = double(B);
    
    % Optional rotation
    if rotate_180 > 0
        R = R(end:-1:1, end:-1:1);
        G = G(end:-1:1, end:-1:1);
        B = B(end:-1:1, end:-1:1);
    end
    
    % Pedestal subtraction
    if pedestal > 0,
        R = max(0, R - pedestal);
        G = max(0, G - pedestal);
        B = max(0, B - pedestal);
        extra_gain = satpoint/(satpoint-pedestal);
        R = R.*extra_gain;
        G = G.*extra_gain;
        B = B.*extra_gain;
    end
    
    fname_wp  = replace_file_extension( fname, wp_file_ext );
    fname_wp = strrep(fname_wp, '_privacy', '');
    wp = load(fname_wp);
    wb_gains = [1/wp(1), 1, 1/wp(2)];
    mingain = min(wb_gains(:));
    wb_gains = wb_gains./mingain;
    try
        fname_ccm = replace_file_extension( fname, '.ccm' );
        fname_ccm = strrep(fname_ccm, '_privacy', '');
        ccm = load(fname_ccm);
    catch
        disp(sprintf('CCM file not found, trying default CCM instead: %s', default_ccm_file));
        ccm = load(default_ccm_file);
    end
    ccm = reshape(ccm, [3,3])';
    gamma = 0.45;
    fname_out = replace_file_extension( fname, out_file_ext );
    
    iwc_illustrate_image( R, G, B, wb_gains, [1080 nan], ccm, bpp, gamma, filelist(i).name, fname_out, max_extra_gain, max_extra_gain_br_ind_def, satpoint, show_figure );
end % for

end % function

function [ ] = iwc_illustrate_image( R, G, B, wb_gains, new_size, ccm, bpp, gamma, title_str, filename, max_extra_gain, max_extra_gain_br_ind_def, satpoint, show_figure )

% Resize if new size has been given
if sum(isnan(new_size)) < 2
    R = imresize(R, new_size, 'bilinear');
    G = imresize(G, new_size, 'bilinear');
    B = imresize(B, new_size, 'bilinear');
end

% Set saturation point according to bpp, if not given explicitly
if isnan(satpoint)
    satpoint = (2^bpp)-1;
end

% Apply WB
img_size = size(R);
tmp = zeros([img_size,3]);
tmp(:,:,1) = double(R).*wb_gains(1);
tmp(:,:,2) = double(G).*wb_gains(2);
tmp(:,:,3) = double(B).*wb_gains(3);
tmp(tmp>satpoint) = satpoint; % Apply clipping to avoid coloration of saturated areas

% Apply CCM if given
if ( (sum(sum(isnan(ccm)))==0) && (~isnan(bpp)) )
    tmp = iwc_apply_ccm(tmp, ccm, bpp, satpoint);
end

% Apply brightness normalization
tmp = iwc_ii_normalize_brightness(tmp, satpoint, max_extra_gain, max_extra_gain_br_ind_def);

% Limit to range [0,1] and apply gamma
tmp = tmp./satpoint;
if gamma < 1.0
    %tmp = tmp.^gamma;
    tmp = srgbgamma(tmp, 'inverse'); % From OptProp; more accurate sRGB gamma in the dark part than the 0.45 exponent approximation
end

% Apply sharpening
tmp = imsharpen(tmp, 'Radius', 1.0, 'Amount', 1.3, 'Threshold', 0);

% Show the image
if show_figure > 0
    figure; imagesc(tmp);
    title(title_str);
end

% Write JPEG if filename has been given
if sum(isnan(filename)) == 0
    imwrite(tmp, filename, 'jpeg', 'Quality', 90);
end

end % function

function [img_out] = iwc_ii_normalize_brightness(img_in, satpoint, max_extra_gain, max_extra_gain_br_ind_def)

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

end % function

function [out_img] = iwc_apply_ccm(in_img, ccm, bpp, satpoint)
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
    if isnan(bpp)
        tmp(tmp>1.0) = 1.0;
    else
        tmp(tmp>satpoint) = satpoint;
    end
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
end % function

function [Gr, R, B, Gb] = iwc_separateChannels(image, order)

if size(image, 3) ~= 1
    error('Input image must be a raw Bayer image');
end % if

switch order
    case 0, [Gr, R, B, Gb] = iwc_convertToChannels(image);
    case 1, [R, Gr, Gb, B] = iwc_convertToChannels(image);
    case 2, [B, Gb, Gr, R] = iwc_convertToChannels(image);
    case 3, [Gb, B, R, Gr] = iwc_convertToChannels(image);
end % switch

end % function

% Manual addition: For wptool_separateChannels
function [ch1, ch2, ch3, ch4] = iwc_convertToChannels(image)

ch1 = image(1:2:end, 1:2:end);
ch2 = image(1:2:end, 2:2:end);
ch3 = image(2:2:end, 1:2:end);
ch4 = image(2:2:end, 2:2:end);

end % function
