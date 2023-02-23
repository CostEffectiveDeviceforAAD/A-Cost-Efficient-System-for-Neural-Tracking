% 
% Curry example script
%
% This m-file's name must not contain spaces.
%
% The input and output arrays are context-dependent.
% input arrays start with the letters 'in'
% input arrays are put in the workspace before the m-file is called
% output arrays start with 'out'
% output arrays need not necessarily be created by the m_file, 
%        but if they are created, precision and dimensions must match expectations
%
% units are uV, fT, uAmm, mm, respectively
%
% The context is indicated by the first three letters of the m-file's name.
% ACQ1 acquisition data stream
% FD1 functional data filtering (pre)
% FD2 functional data filtering (post)
%     indat:          type single, input data, nChannels x nSamples
%     inevents:       type int32, array of input events within current datablock, [sample_wthin_block ; event_type] (only used for contiguous data)
%     inlabels:       type char, list of channel labels (separated by linebreak)
%     insampleratehz: type double, sampling rate of input data in Hz
%     instartsample:  type uint32, absolute startsample of the block
%     inepochtype:    type uint32, event type of epoch (only used for epoched data)
%     outdat:         type single, processed input data, nChannels x nSamples
%     outevents:      type int32, array of output events within current datablock, [sample_wthin_block ; event_type] (only used for contiguous data)
% FD3 functional data event detection and channel selection
% ID1 image data processing
%     inimg:          type uint8,  input data, 256 x 256 x 256
%     inloc:          type double, locations (cursor, reference, NAS, PAL, PAR, AC, PC), 3 x 7
%     outloc:         type double, changed locations (cursor, reference), 3 x 2
% SR1 volume conductor model
%     inloc:          type double, source locations, 3 x nDipoles
%     innor:          type double, source orientations, 3 x nDipoles (column norm must be 1)
%     insen:          type double, sensor locations, 3 x nChannels
%     outlfm:         type double, lead field matrix, nChannels x nDipoles
%     outerr:         type int32, number of locations outside of source space, 1 x 1
% SR2 current density analysis
%     indat:          type double, input data, nChannels x nSamples
%     inlfm:          type double, lead field matrix, nChannels x nDipoles
%     inlam:          type double, lambda vector, nSamples x 1
%     inpar:          type double, UI parameter vector (misfit, pmodel, pdata, qmax), 4 x 1
%     outdat:         type double, forward calculated data matrix, nChannels x nSamples
%     outcdr:         type double, computed currents matrix, nDipoles x nSamples
%     outlam:         type double, used lambda vector, nSamples x 1
% SR3 CDR dipoles
%     inloc:          type double, CDR locations, 3 x nDipoles
%     incdr:          type double, CDR currents matrix, nDipoles x nSamples
%     outloc:         type double, CDR dipole locations, 3 x nLocations
% SR4 Source Coherence
% RE1 Statistics (receive data only)
%

% Get basic info
[numChannels,numSamples] = size(indat);
splittedLabels = strsplit(inlabels, '\n');
chanLabels = cellstr(splittedLabels(1:numChannels));

% Create time vector
startTime = instartsample/insampleratehz;
endTime = double(instartsample+numSamples)/insampleratehz;
t = linspace(double(startTime),double(endTime),numSamples);

% Reorder data so it's plotted as in CURRY
data = flip(indat)';
chanLabels = flip(chanLabels);

% Calculate shift for channels in plot
minmax = [min(min(data)) max(max(data))];
shift = repmat(mean(abs(minmax)), numChannels, 1);
shift = cumsum([0;abs(shift(1:end-1))]);
shift = repmat(shift,1,numSamples)';

% Stacked plot
plot(t, data + shift)
xlabel('Time [s]'), grid on, axis tight
set(gca,'ytick',shift(1,:),'yticklabel',chanLabels, 'fontsize', 8)
ylim([min(min(shift + data)) max(max(shift + data))])
