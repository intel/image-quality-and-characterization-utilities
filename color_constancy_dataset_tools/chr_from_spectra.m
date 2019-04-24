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

function [ RperG, BperG ] = chr_from_spectra( illuminant_txt, illuminant_spd, illuminant_wl, camera_spectral_response_mat, verbose )
%CHR_FROM_SPECTRA Calculate [R/G, B/G] from spectra
%  camera_spectral_response.mat data format:
%    Column 1: Wavelengths
%    Column 2: Gr
%    Column 3: R
%    Column 4: B
%    Column 5: Gb
%
%  illuminant_txt is expected to follow OptProp notation, e.g. 'd65'
%
%  illuminant_spd and illuminant_wl can be used instead of illuminant_txt
%  to specify the illuminant SPD, in which case illuminant_txt should be
%  set to NaN

spectral_data = load(camera_spectral_response_mat);
spectral_data = spectral_data.data;
camera_wl = spectral_data(:,1);
g_sr = (spectral_data(:,2) + spectral_data(:,5))./2;

if sum(isnan(illuminant_txt)) == 0
    ill_spd = illuminant(illuminant_txt, camera_wl);
elseif sum(isnan(illuminant_spd)) == 0 && sum(isnan(illuminant_wl)) == 0
    ill_spd = interp1(illuminant_wl, illuminant_spd, camera_wl, 'linear', 'extrap');
else
    error('Illuminant SPD must be given either as a OptProp text string or numerically along with the wavelengths');
end

if verbose > 0
    figure;
    plot(spectral_data(:,1), ill_spd, 'k-');
    hold on;
    plot(spectral_data(:,1), spectral_data(:,3), 'r-');
    plot(spectral_data(:,1), g_sr, 'g-');
    plot(spectral_data(:,1), spectral_data(:,4), 'b-');
    hold off;
end

R = sum(spectral_data(:,3).*ill_spd);
G = sum(g_sr.*ill_spd);
B = sum(spectral_data(:,4).*ill_spd);

if G == 0
    error('G is zero');
end

RperG = R/G;
BperG = B/G;

end % function

