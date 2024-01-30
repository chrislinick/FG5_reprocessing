% reproFG5.m
% Re-process FG5 observations 
% - Parameterization can be chosen by user:
%   1) conventional 5 terms (z0,v0,g0 + 2 laser demodulation terms)
%   2) conventional 5 terms + Kren et al. (2021) additional 2 laser demodulation terms
% - Drop residuals are generated for all fringe crossings (g9 shows 
%   residuals only for fringe crossings used in estimation of g).
%
% Inputs: exported raw, drop, set, and project files from g9 program
% Outputs: drop g estimates, set g estimates, final g estimates, together
% with associated statistics
%
% Further documentation is in progress.
%
% Author: Chris Linick, The University of Texas at Austin
% Date Created: Feb 2020

%% USER SPECIFIES PATHS TO DATA FILES, AND ESTIMATION SETTINGS IN THIS SECTION

% Retrieve the fringe crossing times from raw file
gtimes = loadrawfile("Merge.raw.txt");    

% Retrieve g9 drop g estimates and environmental corrections from the drop file
gdrop = loaddropfile("Merge.drop.txt");

% Retrieve g9's set g estimates from the set file
gset = loadsetfile("Merge.set.txt"); 

% Retrieve FG5 hardware parameters and g9's estimation settings from project file
[gam,htransfer,lam_laser_vec,fmod,ns,ne,pscale,nsigma_max,gfinal_g9,sfinal_g9] = loadprojectfile("Merge.project.txt");

% Settings for alternate parameterizations
param_method = 2;   % 1 = conventional (5 params); 2 = additional 2 laser demod terms (Kren et al., 2021)

%% USER CAN CHANGE HARDWARE AND ESTIMATE DEFAULTS LOADED FROM THE ABOVE FILES WITHIN THIS SECTION

% Gradient (determined by survey)
gam = +3.079; % ugal/cm

% FG5 hardware configuration constants
lam_laser = 632.99119473; % laser wavelength (nm) (iodine-E)
fmod = 1171.880;             % WEO100 laser mod frequency [Hz]

% FG5 data acquisition 
pscale = 1000;   % every 1000th fringe crossing time recorded

% Settings for g estimation for each drop
ns = 19;    % starting fringe of estimation
ne = 629;   % end fringe of estimation

% Averaging: outlier detection threshold
nsigma_max = 3;     % sigma threshold for determining g outliers

%% Other physical and hardware constants

c = 2.99792458E8;            % speed of light [m/s]
w = 2*pi*fmod;               % angular freq of laser modulation
nfringe = size(gtimes,2);    % number of fringes recorded for each drop
ndrop = length(gdrop.t);  % total number of drops
nset = gdrop.Set(end);    % total number of sets

%% Initializations
cresid = zeros(ndrop,nfringe);    % for saving residuals
cresidt = zeros(ndrop,nfringe);   % for saving adjusted times of each residual

cdrop = table;                  % table for saving g estimates, stdev of resid, and standard error of parameters for each drop
cdrop.t = NaT(ndrop,1);
cdrop.t = gdrop.t;
cdrop.graw = zeros(ndrop,1);
cdrop.gcorr = zeros(ndrop,1);
cdrop.stdresid = zeros(ndrop,1);
cdrop.params = zeros(ndrop,7);
cdrop.serror = zeros(ndrop,7);
cdrop.accept = true(ndrop,1);

cset = table;                   % table for saving set data
cset.mean = zeros(nset,1);
cset.sigma = zeros(nset,1);
cset.naccept = zeros(nset,1);

%% Loop over all drops, estimating parameters and associated statistics
for dnum = 1:ndrop

    % Determine laser wavelength used for drop
    switch(gdrop.LaserLock(dnum))
        case "D"
            lam_laser = lam_laser_vec(1)*1E-9;
        case "E"
            lam_laser = lam_laser_vec(2)*1E-9;
        case "F"
            lam_laser = lam_laser_vec(3)*1E-9;
        case "G"
            lam_laser = lam_laser_vec(4)*1E-9;
        case "H"
            lam_laser = lam_laser_vec(5)*1E-9;
        case "I"
            lam_laser = lam_laser_vec(6)*1E-9;
        case "J"
            lam_laser = lam_laser_vec(7)*1E-9;
        otherwise
            error("Laser wavelength is not an expected value in g9 drop file. WEO100 is only valid type for this software currently")
    end

    d = (1:nfringe)'*lam_laser*(pscale/2);  % distances fallen (z's) at each fringe crossing

    t0 = gtimes(dnum,:)';  % fringe times for this drop

    % Initial state vector
    if (param_method == 1)
        m0 = [ 0; 0; 9.8; 0; 0];
    elseif (param_method == 2)
        m0 = [ 0; 0; 9.8; 0; 0; 0; 0];
    else
        disp("param_method invalid; must be 1 or 2.")
    end

    % Shift times to adjusted time (to correct for TOF)
    z0 = m0(1);
    ta = t0 + (d-z0)/c;  % adjusted times

    % Calculating time powers to simplify G matrix
    ta2 = ta.^2;
    ta3 = ta.^3;
    ta4 = ta.^4;

    % Form partials matrix G (model that maps params to obs)
    if (param_method == 1)
        % Conventional parameterization
        G = [ 1+0.5*gam*ta2 ta+(1/6)*gam*ta3 0.5*(ta2+(1/12)*gam*ta4) ...
            cos(w*ta) sin(w*ta)];
    else
        %  Kren et al. 2021 augmented terms and 2 additional terms for
        %  laser demodulation
        G = [ 1+0.5*gam*ta2 ta+(1/6)*gam*ta3 0.5*(ta2+(1/12)*gam*ta4) ...
            sin(w*ta).*cos(w*2*d/c) cos(w*ta).*cos(w*2*d/c)...
            sin(w*ta).*sin(w*2*d/c) cos(w*ta).*sin(w*2*d/c)];
    end

    % Calculate LS estimate of parameters m
    mhat = G(ns:ne,:)\d(ns:ne); % use only subset of G and d in the chosen window

    % Calculate residuals (O-C)
    r = d - G*mhat;     % full window

    % Save residuals and parameters in table
    cresid(dnum,:) = r'*1E9;            % residual distance (nm)
    cresidt(dnum,:) = ta';               % time (sec)
    cdrop.graw(dnum) = mhat(3)*1E8;     % drop gravity estimate (ugal)
    cdrop.params(dnum,:) = mhat';          % estimated parameters

    % Calculate covariance matrix and standard errors
    nobs = length(ns:ne);           % num observations
    npar = size(G,2);               % num parameters
    vn = (r'*r)/(nobs-npar);        % noise variance estimate
    Cmat = pinv(G'*G)*vn;              % covariance matrix (GtG)^-1 * variance of residual (may need to check assumptions here for LS)
    cdrop.serror(dnum,:) = sqrt(diag(Cmat))';        % standard errors are sqrt of diagonal elements (meter is distance unit)

    % Stdev residual (of window region only); this is the "Error" column in g9
    % drop files
    cdrop.stdresid(dnum) = std(r(ns:ne)*1E9);    % units of nm

end

%% Correct raw g estimates for environment
dg = gdrop.Tide+gdrop.Load+gdrop.Baro+...
    gdrop.Polar+gdrop.Transfer+gdrop.Refxo+gdrop.Tilt;

cdrop.gcorr = cdrop.graw + dg;

% Now form set averages by detecting outliers in each set greater than
% nsigma_max from the mean (removing outliers and then repeating until
% there are no more outliers)
for j = 1:nset

    gcj = cdrop.gcorr((j-1)*100+1:j*100);
    jkeep = true(100,1);
    noutlier_old = -1;
    noutlier_new = 0;

    while (noutlier_new > noutlier_old)

        mj = mean(gcj(jkeep));
        sj = std(gcj(jkeep));
        nsigma = abs((gcj-mj)/sj);

        jkeep = nsigma < nsigma_max;
        cdrop.accept((j-1)*100+1:j*100) = jkeep;
        noutlier_old = noutlier_new;
        noutlier_new = 100 - nnz(jkeep);

    end

    cset.mean(j) = mj;
    cset.sigma(j) = sj;
    cset.naccept(j) = nnz(jkeep);

end

% Form final gravity estimate and std (both are weighted)
wfinal = (1./cset.sigma).^2;
gfinal = (cset.mean'*wfinal)/sum(wfinal);
sfinal = std(cset.mean,1./(cset.sigma.^2));

change_g =  gfinal - gfinal_g9;
change_std = sfinal - sfinal_g9;