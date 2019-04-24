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

function [figh] = plot_wps()
%PLOT_WPS Plot white point coordinates from .wp files

% INPUT --> 

cameras = cell(4,1);

data.camera_name = 'Sony IMX135 Semco w. BLC and CSC';
data.dirs_to_be_processed = {'..\Sony_IMX135_BLCCSC\Sony_IMX135_BLCCSC_field\',...
                             '..\Sony_IMX135_BLCCSC\Sony_IMX135_BLCCSC_field_part2\',...
                             '..\Sony_IMX135_BLCCSC\Sony_IMX135_BLCCSC_lab_printouts\',...
                             '..\Sony_IMX135_BLCCSC\Sony_IMX135_BLCCSC_lab_realscene\'};
data.plotstr = 'b^';
data.plotcolor = [0 0 0.5];
cameras{1} = data;

data.camera_name = 'Canon5DSR';
data.dirs_to_be_processed = {'..\Canon_5DSR\Canon_5DSR_field\',...
                             '..\Canon_5DSR\Canon_5DSR_field_part2\',...
                             '..\Canon_5DSR\Canon_5DSR_lab_printouts\',...
                             '..\Canon_5DSR\Canon_5DSR_lab_realscene\'};
data.plotstr = 'ro';
data.plotcolor = [0.5 0 0];
cameras{2} = data;

data.camera_name = 'NikonD810';
data.dirs_to_be_processed = {'..\Nikon_D810\Nikon_D810_field\',...
                             '..\Nikon_D810\Nikon_D810_field_part2\',...
                             '..\Nikon_D810\Nikon_D810_lab_printouts\',...
                             '..\Nikon_D810\Nikon_D810_lab_realscene\'};
data.plotstr = 'gs';
data.plotcolor = [0 0.5 0];
cameras{3} = data;

data.camera_name = 'Canon5DSR field2';
data.dirs_to_be_processed = {'..\Canon_5DSR\Canon_5DSR_field2\',...
                             '..\Canon_5DSR\Canon_5DSR_field2_part2\'};
data.plotstr = 'k.';
data.plotcolor = [0 0 0];
cameras{4} = data;

selection_wildcard = '*.wp';
process_subdirs = 1;
max_wp_cnt = 5000;

special_chromaticities_txt = {'SL\_D', 'SL\_CW', 'SL\_Hor', 'SL\_A', 'SL\_TL84', 'IE\_D65', 'IE\_D50', 'IE\_F11', 'IE\_F12', 'IE\_A'};
special_chromaticities_Canon5DSR = [ ...
    0.3971, 0.6626; ... % SL_D
    0.5402, 0.4280; ... % SL_CW
    0.9071, 0.2885; ... % SL_Hor
    0.7261, 0.3539; ... % SL_A
    0.5402, 0.4383; ... % SL_TL84
    0.4092, 0.6283; ... % IE_D65
    0.4533, 0.5799; ... % IE_D50
    0.5314, 0.4759; ... % IE_F11
    0.6704, 0.3755; ... % IE_F12
    0.7264, 0.3506  ... % IE_A
    ];
special_chromaticities_NikonD810 = [ ...
    0.4492, 0.8365; ... % SL_D
    0.5600, 0.5186; ... % SL_CW
    1.0534, 0.3415; ... % SL_Hor
    0.8489, 0.4356; ... % SL_A
    0.6183, 0.5420; ... % SL_TL84
    0.4735, 0.7836; ... % IE_D65
    0.5523, 0.7554; ... % IE_D50
    0.6307, 0.5992; ... % IE_F11
    0.8088, 0.4552; ... % IE_F12
    0.8479, 0.4287  ... % IE_A
    ];
special_chromaticities_SonyIMX135 = [ ...
    0.5226, 0.7237; ... % SL_D
    0.6162, 0.4981; ... % SL_CW
    1.1476, 0.3636; ... % SL_Hor
    0.9343, 0.4129; ... % SL_A
    0.6602, 0.4912; ... % SL_TL84
    0.5201, 0.6553; ... % IE_D65
    0.5882, 0.6066; ... % IE_D50
    0.6693, 0.5120; ... % IE_F11
    0.8258, 0.4138; ... % IE_F12
    0.9171, 0.3998 ... % IE_A
    ];
special_chromaticities_cnt = length(special_chromaticities_Canon5DSR);
% <-- INPUT

if 0
    load('Spectral_light_source_data.mat');
    wl = [362:782];
    special_chromaticities_Canon5DSR_SPD = zeros(10,2);
    [chr1, chr2] = chr_from_spectra(NaN, SL_D, wl, 'sr_canon5dsr.mat', 0); special_chromaticities_Canon5DSR_SPD(1,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_CW, wl, 'sr_canon5dsr.mat', 0); special_chromaticities_Canon5DSR_SPD(2,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_Hor, wl, 'sr_canon5dsr.mat', 0); special_chromaticities_Canon5DSR_SPD(3,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_A, wl, 'sr_canon5dsr.mat', 0); special_chromaticities_Canon5DSR_SPD(4,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_TL84, wl, 'sr_canon5dsr.mat', 0); special_chromaticities_Canon5DSR_SPD(5,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_D65, wl, 'sr_canon5dsr.mat', 0); special_chromaticities_Canon5DSR_SPD(6,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_D50, wl, 'sr_canon5dsr.mat', 0); special_chromaticities_Canon5DSR_SPD(7,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_F11, wl, 'sr_canon5dsr.mat', 0); special_chromaticities_Canon5DSR_SPD(8,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_F12, wl, 'sr_canon5dsr.mat', 0); special_chromaticities_Canon5DSR_SPD(9,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_A, wl, 'sr_canon5dsr.mat', 0); special_chromaticities_Canon5DSR_SPD(10,:) = [chr1, chr2];
    special_chromaticities_NikonD810_SPD = zeros(10,2);
    [chr1, chr2] = chr_from_spectra(NaN, SL_D, wl, 'sr_nikond810.mat', 0); special_chromaticities_NikonD810_SPD(1,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_CW, wl, 'sr_nikond810.mat', 0); special_chromaticities_NikonD810_SPD(2,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_Hor, wl, 'sr_nikond810.mat', 0); special_chromaticities_NikonD810_SPD(3,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_A, wl, 'sr_nikond810.mat', 0); special_chromaticities_NikonD810_SPD(4,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_TL84, wl, 'sr_nikond810.mat', 0); special_chromaticities_NikonD810_SPD(5,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_D65, wl, 'sr_nikond810.mat', 0); special_chromaticities_NikonD810_SPD(6,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_D50, wl, 'sr_nikond810.mat', 0); special_chromaticities_NikonD810_SPD(7,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_F11, wl, 'sr_nikond810.mat', 0); special_chromaticities_NikonD810_SPD(8,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_F12, wl, 'sr_nikond810.mat', 0); special_chromaticities_NikonD810_SPD(9,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_A, wl, 'sr_nikond810.mat', 0); special_chromaticities_NikonD810_SPD(10,:) = [chr1, chr2];
    special_chromaticities_SonyIMX135_SPD = zeros(10,2);
    [chr1, chr2] = chr_from_spectra(NaN, SL_D, wl, 'sr_sonyimx135.mat', 0); special_chromaticities_SonyIMX135_SPD(1,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_CW, wl, 'sr_sonyimx135.mat', 0); special_chromaticities_SonyIMX135_SPD(2,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_Hor, wl, 'sr_sonyimx135.mat', 0); special_chromaticities_SonyIMX135_SPD(3,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_A, wl, 'sr_sonyimx135.mat', 0); special_chromaticities_SonyIMX135_SPD(4,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, SL_TL84, wl, 'sr_sonyimx135.mat', 0); special_chromaticities_SonyIMX135_SPD(5,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_D65, wl, 'sr_sonyimx135.mat', 0); special_chromaticities_SonyIMX135_SPD(6,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_D50, wl, 'sr_sonyimx135.mat', 0); special_chromaticities_SonyIMX135_SPD(7,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_F11, wl, 'sr_sonyimx135.mat', 0); special_chromaticities_SonyIMX135_SPD(8,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_F12, wl, 'sr_sonyimx135.mat', 0); special_chromaticities_SonyIMX135_SPD(9,:) = [chr1, chr2];
    [chr1, chr2] = chr_from_spectra(NaN, IE_A, wl, 'sr_sonyimx135.mat', 0); special_chromaticities_SonyIMX135_SPD(10,:) = [chr1, chr2];
    chr1 = special_chromaticities_Canon5DSR(1:10,:)';
    chr2 = special_chromaticities_Canon5DSR_SPD';
    err = 100*((chr1./chr2)-1);
    disp(sprintf('Canon5DSR [R/G,B/G] errors (RAW vs. SPD -based):'));
    disp(sprintf('%.1f%%, %.1f%%\n', err));
    err = err';
    disp(sprintf('Mean errors: %.1f%%, %.1f%%', mean(err(:,1)), mean(err(:,2))));
    disp(sprintf('Error range: [%.1f%%, %.1f%%]\n', min(err(:)), max(err(:))));
    chr1 = special_chromaticities_NikonD810(1:10,:)';
    chr2 = special_chromaticities_NikonD810_SPD';
    err = 100*((chr1./chr2)-1);
    disp(sprintf('NikonD810 [R/G,B/G] errors (RAW vs. SPD -based):'));
    disp(sprintf('%.1f%%, %.1f%%\n', err));
    err = err';
    disp(sprintf('Mean errors: %.1f%%, %.1f%%', mean(err(:,1)), mean(err(:,2))));
    disp(sprintf('Error range: [%.1f%%, %.1f%%]\n', min(err(:)), max(err(:))));
    chr1 = special_chromaticities_SonyIMX135(1:10,:)';
    chr2 = special_chromaticities_SonyIMX135_SPD';
    err = 100*((chr1./chr2)-1);
    disp(sprintf('SonyIMX135 [R/G,B/G] errors (RAW vs. SPD -based):'));
    disp(sprintf('%.1f%%, %.1f%%\n', err));
    err = err';
    disp(sprintf('Mean errors: %.1f%%, %.1f%%', mean(err(:,1)), mean(err(:,2))));
    disp(sprintf('Error range: [%.1f%%, %.1f%%]\n', min(err(:)), max(err(:))));
end

figh = figure;
hold on;
legend_strs = cell(length(cameras),1);
for j=1:length(cameras)
    data = cameras{j};
    disp(sprintf('Processing %s', data.camera_name));
    wplist = zeros(max_wp_cnt,2);
    wpind = 0;
    for i=1:length(data.dirs_to_be_processed)
        [wplist, wpind] = process_directory(data.dirs_to_be_processed{i}, selection_wildcard, process_subdirs, wplist, wpind);
    end % for i
    plot(wplist(1:wpind,1), wplist(1:wpind,2), data.plotstr, 'MarkerSize', 5, 'MarkerEdgeColor', data.plotcolor);
    legend_strs{j} = data.camera_name;
end % for j
legend(legend_strs, 'Location', 'NorthEast');
xlabel('R/G');
ylabel('B/G');
plot(special_chromaticities_Canon5DSR(:,1), special_chromaticities_Canon5DSR(:,2), 'ro', 'MarkerSize', 7, 'LineWidth', 1, 'MarkerEdgeColor', [0.25,0,0], 'MarkerFaceColor', [0.5,0,0]);
plot(special_chromaticities_NikonD810(:,1), special_chromaticities_NikonD810(:,2), 'gs', 'MarkerSize', 7, 'LineWidth', 1, 'MarkerEdgeColor', [0,0.25,0], 'MarkerFaceColor', [0,0.5,0]);
plot(special_chromaticities_SonyIMX135(:,1), special_chromaticities_SonyIMX135(:,2), 'b^', 'MarkerSize', 7, 'LineWidth', 1, 'MarkerEdgeColor', [0,0,0.25], 'MarkerFaceColor', [0,0,0.5]);

for j=1:special_chromaticities_cnt
    text(special_chromaticities_Canon5DSR(j,1)+0.01, special_chromaticities_Canon5DSR(j,2)+0.01, special_chromaticities_txt{j}, 'Color', [0.5,0,0]);
    text(special_chromaticities_NikonD810(j,1)+0.01, special_chromaticities_NikonD810(j,2)+0.01, special_chromaticities_txt{j}, 'Color', [0,0.5,0]);
    text(special_chromaticities_SonyIMX135(j,1)+0.01, special_chromaticities_SonyIMX135(j,2)+0.01, special_chromaticities_txt{j}, 'Color', [0,0,0.5]);
end % for j

hold off;

end % function

function [wplist, wpind] = process_directory(input_dir, selection_wildcard, process_subdirs, wplist, wpind)

filelist = dir([input_dir, selection_wildcard]);

for i = 1:length(filelist),
    fname = [input_dir, filelist(i).name];
    
    wp = load(fname); % Read the white point from file
    wpind = wpind + 1;
    wplist(wpind,:) = wp;
    
end % for i

if process_subdirs
    filelist = dir([input_dir, '*']);
    input_dir_org = input_dir;
    
    % Recursively call process_directory() for each subdirectory
    for i = 1:length(filelist),
        if filelist(i).isdir && ~strcmp(filelist(i).name,'.') && ~strcmp(filelist(i).name,'..')
            input_dir = [input_dir_org, filelist(i).name, '\'];
            [wplist, wpind] = process_directory(input_dir, selection_wildcard, process_subdirs, wplist, wpind);
        end % if isdir
    end % for i
end % if process_subdirs

end % function
