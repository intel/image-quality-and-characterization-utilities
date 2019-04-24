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

% NOTE: Set path so that "parse_makernote" can be called (for .I3AV4 support)

% The raw files
root_raw_folders = { ...
    '..\example_raw_folder_1\'; ...
    '..\example_raw_folder_2\'  ...
    };

% The files that contain the metadata (e.g. EXIF data)
root_meta_folders = { ...
    '..\example_meta_folder_1\'; ...
    '..\example_meta_folder_2\'  ...
    };

% Type of the metadata (.jpg or .i3av4)
meta_types = { ...
    '.jpg'; ...      % JPEG EXIF
    '.jpg'  ...      % JPEG EXIF
    };

% Camera base ISO etc.
other_characterization = { ...
    'other_characterization_xxx.mat'; ... % Set according to camera type
    'other_characterization_xxx.mat' ...  % Set according to camera type
    };

% Wildcard for selecting the raw files
selection_wildcard = '*.plain16';

% <-- INPUTS

disp('Generating .meta files for the given folders...');

for i = 1:length(root_raw_folders)
    dir_to_be_processed = root_raw_folders{i};
    fprintf('Processing folder %s\n', dir_to_be_processed);

    filelist = dir([dir_to_be_processed, '**\', selection_wildcard]); % include subfolders
    cnt_ok = 0;
    cnt_nok = 0;
    data = load(other_characterization{i});
    base_iso = data.base_iso;
    base_iso_aperture = data.base_iso_aperture;
    unity_gain_iso = data.unity_gain_iso;

    for j = 1:length(filelist)
        folder_raw = filelist(j).folder;
        name_raw = filelist(j).name;
        name_meta = replace_file_extension(name_raw, '.meta');
        folder_metasrc = root_meta_folders{i};

        metasrc_search_name = replace_file_extension(name_raw, '');
        metasrc_search_name = replace_file_extension(metasrc_search_name, ''); % to remove second extension e.g. from .cr2.plain16
        metasrc_search_name = strrep(metasrc_search_name, '_blc_lsc', '');
        metasrc_search_name = strrep(metasrc_search_name, '_SpectraLight_D', '');
        metasrc_search_name = strrep(metasrc_search_name, '_SpectraLight_CW', '');
        metasrc_search_name = strrep(metasrc_search_name, '_SpectraLight_Hor', '');
        metasrc_search_name = strrep(metasrc_search_name, '_SpectraLight_A', '');
        metasrc_search_name = strrep(metasrc_search_name, '_SpectraLight_TL84', '');
        metasrc_search_name = strrep(metasrc_search_name, '_LightSTUDIO_D65', '');
        metasrc_search_name = strrep(metasrc_search_name, '_LightSTUDIO_D50', '');
        metasrc_search_name = strrep(metasrc_search_name, '_LightSTUDIO_F11', '');
        metasrc_search_name = strrep(metasrc_search_name, '_LightSTUDIO_F12', '');
        metasrc_search_name = strrep(metasrc_search_name, '_LightSTUDIO_A', '');
        metasrc_search_name = strcat(metasrc_search_name, '*', meta_types{i});

        metasrc_file_list = dir([folder_metasrc, '\**\', metasrc_search_name]);

        if ~isempty(metasrc_file_list)
            for file_index = 1:length(metasrc_file_list)
                folder_metasrc = metasrc_file_list(file_index).folder;
                name_metasrc = metasrc_file_list(file_index).name;
                metasrc_filename = [folder_metasrc, '\', name_metasrc];
                raw_filename = [folder_raw, '\', name_raw];
                meta_type = meta_types{i};

                % Read metadata from metasrc_filename
                [exposure_time_s, analog_gain, iso, aperture, normalized_exposure_s] = read_metadata(metasrc_filename, base_iso, base_iso_aperture, unity_gain_iso, raw_filename, meta_type);

                % Step out from the for loop if valid parameters were extracted
                if sum([~isnan(exposure_time_s), ~isnan(analog_gain), ~isnan(iso), ~isnan(aperture), ~isnan(normalized_exposure_s)]) > 0
                    break;
                end
            end % for file_index

            % Update .meta file if there are valid parameters to update
            if sum([~isnan(exposure_time_s), ~isnan(analog_gain), ~isnan(iso), ~isnan(aperture), ~isnan(normalized_exposure_s)]) > 0
                update_meta_file(exposure_time_s, analog_gain, iso, aperture, normalized_exposure_s, [folder_raw, '\', name_meta]);
                cnt_ok = cnt_ok + 1;
            else
                cnt_nok = cnt_nok + 1;
            end % if there are valid values to be written
        else
            % Metadata file not found
            cnt_nok = cnt_nok + 1;
            fprintf('Metadata not found for %s (searching for %s)\n', [folder_raw, '\', name_raw], [folder_metasrc, '\**\', metasrc_search_name]);
        end
    end % for j

    fprintf('Folder %s: %d valid .meta files processed, %d failed\n', root_raw_folders{i}, cnt_ok, cnt_nok);

end % for i

function [normalized_exposure] = calculate_normalized_exposure(exposure_time_s, gain, aperture, base_iso, base_iso_aperture)

% Set aperture gain according to current aperture vs. base ISO aperture
aperture_gain = 1.0;
if aperture > 0 && base_iso_aperture > 0 && aperture ~= base_iso_aperture
    aperture_gain = (base_iso_aperture^2)/(aperture^2);
end

normalized_exposure = exposure_time_s * gain * (base_iso/100) * aperture_gain;

end % function

function [] = update_meta_file(exposure_time_s, analog_gain, iso, aperture, normalized_exposure_s, meta_filename)

% Other contents of the .meta file must be preserved if the
% file already exists. Otherwise just write the valid values.
if exist(meta_filename, 'file') == 2
    fid = fopen(meta_filename,'r');
    meta_struct = textscan(fid,'%s%f');
    fclose(fid);
    exposure_time_processed = false;
    analog_gain_processed = false;
    iso_processed = false;
    aperture_processed = false;
    normalized_exposure_s_processed = false;
    for k = 1:length(meta_struct{1})
        if strcmp(meta_struct{1}(k),'exposure_time')
            if ~isnan(exposure_time_s)
                meta_struct{2}(k) = exposure_time_s;
            end
            exposure_time_processed = true;
        elseif strcmp(meta_struct{1}(k),'analog_gain')
            if ~isnan(analog_gain)
                meta_struct{2}(k) = analog_gain;
            end
            analog_gain_processed = true;
        elseif strcmp(meta_struct{1}(k),'iso')
            if ~isnan(iso)
                meta_struct{2}(k) = iso;
            end
            iso_processed = true;
        elseif strcmp(meta_struct{1}(k),'aperture')
            if ~isnan(aperture)
                meta_struct{2}(k) = aperture;
            end
            aperture_processed = true;
        elseif strcmp(meta_struct{1}(k),'normalized_exposure_s')
            if ~isnan(normalized_exposure_s)
                meta_struct{2}(k) = normalized_exposure_s;
            end
            normalized_exposure_s_processed = true;
        end
    end % for k

    fid = fopen(meta_filename,'w');
    for k = 1:length(meta_struct{1})
        fprintf(fid, '%s\t%f\n', meta_struct{1}{k}, meta_struct{2}(k));
    end % for k
    if ~exposure_time_processed
        fprintf(fid, '%s\t%f\n', 'exposure_time', exposure_time_s);
    end
    if ~iso_processed
        fprintf(fid, '%s\t%f\n', 'iso', iso);
    end
    if ~analog_gain_processed
        fprintf(fid, '%s\t%f\n', 'analog_gain', analog_gain);
    end
    if ~aperture_processed
        fprintf(fid, '%s\t%f\n', 'aperture', aperture);
    end
    if ~normalized_exposure_s_processed
        fprintf(fid, '%s\t%f\n', 'normalized_exposure_s', normalized_exposure_s);
    end
    fclose(fid);
else
    fid = fopen(meta_filename,'w');
    if ~isnan(exposure_time_s)
        fprintf(fid, '%s\t%f\n', 'exposure_time', exposure_time_s);
    end
    if ~isnan(iso)
        fprintf(fid, '%s\t%f\n', 'iso', iso);
    end
    if ~isnan(analog_gain)
        fprintf(fid, '%s\t%f\n', 'analog_gain', analog_gain);
    end
    if ~isnan(aperture)
        fprintf(fid, '%s\t%f\n', 'aperture', aperture);
    end
    if ~isnan(normalized_exposure_s)
        fprintf(fid, '%s\t%f\n', 'normalized_exposure_s', normalized_exposure_s);
    end
    fclose(fid);
end % if file exists or not

end % function

function [exposure_time_s, analog_gain, iso, aperture, normalized_exposure_s] = read_metadata(metasrc_filename, base_iso, base_iso_aperture, unity_gain_iso, raw_filename, meta_type)

% .I3AV4 MakerNote data
if strcmp(meta_type, '.i3av4')
    try
        mkn = parse_makernote(metasrc_filename); % parse_makernote needs to be in Matlab path
        exposure_time_s = double(mkn.E3A_DNID_AE_Result_Exposure_Time.data)/1000000.0;
        analog_gain = double(mkn.E3A_DNID_AE_Result_Analog_Gain.data);
        iso = NaN;
        aperture = base_iso_aperture;

        normalized_exposure_s = calculate_normalized_exposure(exposure_time_s, analog_gain, aperture, base_iso, base_iso_aperture);
    catch ME
        fprintf('Error reading the MakerNote for %s (from %s)\n', raw_filename, metasrc_filename);
        fprintf('(%s)\n', ME.identifier);

        exposure_time_s = NaN;
        analog_gain = NaN;
        iso = NaN;
        aperture = NaN;
        normalized_exposure_s = NaN;
    end
% JPEG EXIF data
else
    try
        exif = imfinfo(metasrc_filename);
        exposure_time_s = exif.DigitalCamera.ExposureTime;
        iso = exif.DigitalCamera.ISOSpeedRatings;
        aperture = exif.DigitalCamera.FNumber;
        analog_gain = NaN;

        normalized_exposure_s = calculate_normalized_exposure(exposure_time_s, iso/unity_gain_iso, aperture, base_iso, base_iso_aperture);
    catch ME
        fprintf('Error reading the MakerNote for %s (from %s)\n', raw_filename, metasrc_filename);
        fprintf('(%s)\n', ME.identifier);

        exposure_time_s = NaN;
        analog_gain = NaN;
        iso = NaN;
        aperture = NaN;
        normalized_exposure_s = NaN;
    end
end

end % function
