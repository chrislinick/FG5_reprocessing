function [gamma,htransfer,lam_laser,fmod,ns,ne,pscale,nsigma,gfinal,sfinal] = loadprojectfile(fname)
%LOADPROJECTFILE Loads FG5 configuration from g9 project file.
% Parses a g9 project file to retrieve configuration parameters used for
% re-processing raw drop data into gravity estimates.
%
% INPUT:
% fname = name of project file (include path to file if running from
% different directory)
%
% OUTPUTS:
% gamma = measured gravitational gradient
% htransfer = transfer height
% lam_laser = center frequencies (WEO100 laser only, in order ID/IE/IF/IG/IH/II/IJ)
% fmod = laser modulation frequency (WEO100 laser only)
% ns = number of starting fringe in observation "window" used for estimation
% ne = number of ending fringe in oservation "window" used for estimation
% pscale = data acquisition pre-scale factor (multiplex * scale)
% nsigma = sigma reject
% gfinal = final estimate of g
% sfinal = standard dev of final estimate of g
%
% Author: Christopher Linick
% Date created: 2023-07-18
% Comments:
%
%

% Test whether input project file exists
assert(exist(fname,"file"),sprintf("%s %s", "Project file not found. Check filename (and path, if needed): ",fname))

fid = fopen(fname);
tline = fgetl(fid);
while ischar(tline)
    
    t2 = append(tline,'                    '); % pad with spaces to make the parsing script work

    % Skip empty lines (only spaces)
    if isempty(strtrim(t2))

    % Measured gradient    
    elseif (strcmp("Gradient", t2(1:8)) && ~exist("gamma","var"))
        gamma = str2double(t2(11:16));

    % Transfer height
    elseif (strcmp("Transfer Height", t2(1:15)) && ~exist("htransfer","var"))
        htransfer = str2double(tline(17:end-2));

    % WEO100 laser center frequencies
    elseif strcmp("ID",t2(1:2))
        lam_laser(1) = str2double(t2(5:16));
    elseif strcmp("IE",t2(1:2))
        lam_laser(2) = str2double(t2(5:16));
    elseif strcmp("IF",t2(1:2))
        lam_laser(3) = str2double(t2(5:16));
    elseif strcmp("IG",t2(1:2))
        lam_laser(4) = str2double(t2(5:16));
    elseif strcmp("IH",t2(1:2))
        lam_laser(5) = str2double(t2(5:16));
    elseif strcmp("II",t2(1:2))
        lam_laser(6) = str2double(t2(5:16));
    elseif strcmp("IJ",t2(1:2))
        lam_laser(7) = str2double(t2(5:16));

    % WE100 laser modulation frequencies
    elseif strcmp("Modulation",t2(1:10))
        fmod = str2double(tline(22:end-2));

    % Start fringe number
    elseif strcmp("Fringe Start",t2(1:12))
        ns = str2double(t2(14:end));

    % End fringe number
    elseif strcmp("Processed Fringes",t2(1:17))
        ne = str2double(t2(19:end))+ns-1;

    % Pre-scale factor settings
    elseif strcmp("GuideCard Multiplex",t2(1:19))
        gmx = str2double(t2(21:end));
    elseif strcmp("GuideCard Scale",t2(1:15))
        gsf = str2double(t2(24:end));

    % Sigma reject (for detecting outliers before averaging)
    elseif strcmp("Sigma Reject",t2(1:12))
        nsigma = str2double(t2(14:end));

    % Final g value
    elseif strcmp("Gravity:",t2(1:8))
        gfinal = str2double(tline(9:end-4));

    % Set scatter
    elseif strcmp("Set Scatter:",t2(1:12))
        sfinal = str2double(tline(13:end-5));

    end

    tline = fgetl(fid);
end

% Compute prescale factor from the acq card settings
if exist("gmx","var") && exist("gsf","var")
    pscale = gmx*gsf;
else
    error("Unable to determine prescale factor from project file")
end


fclose(fid);
end