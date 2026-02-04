%% quantify_response_peaks.m
% Quantify interval-specific responses from calcium imaging mean traces using OASIS,
% with plotting of traces + interval shading (GitHub-shareable).
%
% Requirements (on MATLAB path):
%   - OASIS / CaImAn MATLAB functions: deconvolveCa, GetSn, estimate_time_constant
%
% Expected .mat contents:
%   - mean_trace (vector)
%   - OPTIONAL: all_traces (nROIs x T) if you want multi-trace plotting
%
% Output:
%   - max_peaks: [nFiles x nIntervals]
%   - integrals: [nFiles x nIntervals]

clearvars; close all; clc;

%% ------------------- User settings -------------------
DATA_DIR     = 'D:\Dropbox\halotag\halotag_shared\CaImaging_gfi\revision\dots\';
FILES        = {'max_normalized_data.mat'};   % list of .mat files to process

SCALE_FACTOR = 2;       % 1 = no upsampling, 2 = 2x, etc.
PLOT_SUMMARY = true;    % bar plot of max_peaks per file
PLOT_FITS    = true;    % OASIS fit diagnostic plots
PLOT_TRACES  = true;    % plot mean_trace (and all_traces if available) + interval shading

% "Revision dots" intervals (MANUAL ADJUSTMENT) in original (non-upsampled) timebase
intervals = [ ...
    144, 157; ... % first dark background
    158, 178; ... % dot1
    195, 221; ... % dot2
    246, 258; ... % dark back - appearance grating
    260, 275; ... % grating1
    292, 305; ... % grating2
    313, 333; ... % grating3
    342, 355; ... % grating4
    356, 370  ... % dark background
];

% Shift and scale (as in your original script)
intervals = (intervals - 75) * SCALE_FACTOR;
%% -----------------------------------------------------

nFiles = numel(FILES);
nInt   = size(intervals, 1);

integrals = nan(nFiles, nInt);
max_peaks = nan(nFiles, nInt);

for ff = 1:nFiles
    inFile = fullfile(DATA_DIR, FILES{ff});
    S = load(inFile);

    if ~isfield(S, 'mean_trace')
        error('File "%s" does not contain variable "mean_trace".', inFile);
    end

    % Mean trace: column + baseline correction
    y0 = double(S.mean_trace(:));
    y0 = y0 - prctile(y0, 5);

    % Optional upsampling
    if SCALE_FACTOR ~= 1
        x  = 1:numel(y0);
        xq = linspace(1, numel(y0), numel(y0) * SCALE_FACTOR);
        y  = interp1(x, y0, xq, 'pchip')';
    else
        y = y0;
    end

    % Optional: all_traces handling (only for plotting)
    hasAllTraces = isfield(S, 'all_traces') && ~isempty(S.all_traces);
    if hasAllTraces
        all_traces0 = double(S.all_traces);
        % baseline-correct each trace similarly (5th percentile per trace)
        p5 = prctile(all_traces0, 5, 2);
        all_traces0 = all_traces0 - p5;

        if SCALE_FACTOR ~= 1
            T0 = size(all_traces0, 2);
            x  = 1:T0;
            xq = linspace(1, T0, T0 * SCALE_FACTOR);

            all_traces = zeros(size(all_traces0, 1), numel(xq));
            for i = 1:size(all_traces0, 1)
                all_traces(i, :) = interp1(x, all_traces0(i, :), xq, 'pchip');
            end
        else
            all_traces = all_traces0;
        end
    end

    % Fit + quantify
    [integrals(ff, :), max_peaks(ff, :), fitOut] = oasisFIT(y, intervals, PLOT_FITS);

    % Plot traces with interval shading (mean trace + optional all_traces)
    if PLOT_TRACES
        if hasAllTraces
            plotAllTracesWithIntervals(all_traces, y, intervals, FILES{ff});
        else
            plotMeanTraceWithIntervals(y, intervals, FILES{ff});
        end

        % Also show mean trace with reconstructed components (nice single panel)
        plotReconstructionSummary(y, intervals, fitOut, FILES{ff});
    end

    % Summary bar plot
    if PLOT_SUMMARY
        figure('Name', sprintf('Max peaks: %s', FILES{ff}), 'Color', 'w');
        bar(max_peaks(ff, :));
        xlabel('Interval index');
        ylabel('Max peak (reconstructed c)');
        title(strrep(FILES{ff}, '_', '\_'));
        grid on;
    end
end

%% -------------------- Local functions --------------------
function plotMeanTraceWithIntervals(y, intervals, labelStr)
figure('Name', ['Mean trace + intervals: ' labelStr], 'Color', 'w');
plot(y, 'k', 'LineWidth', 1); hold on;
shadeIntervals(gca, intervals, [min(y), max(y)]);
xlabel('Time (samples)');
ylabel('Amplitude');
title(['Mean trace with intervals: ' strrep(labelStr, '_', '\_')]);
grid on;
end

function plotAllTracesWithIntervals(all_traces, mean_trace, intervals, labelStr)
% Similar to your original "PLOT_ALL" section: 3x2 grid up to 6 traces.
nToShow = min(size(all_traces, 1), 6);
nRows = 3; nCols = 2;

figure('Name', ['All traces + intervals: ' labelStr], 'Color', 'w');
for k = 1:nToShow
    subplot(nRows, nCols, k);
    plot(all_traces(k, :), 'k'); hold on;
    yl = [min(mean_trace), max(mean_trace)];
    shadeIntervals(gca, intervals, yl);
    title(sprintf('Trace %d', k));
    xlabel('Time (samples)');
    ylabel('Amplitude');
    grid on;
end
end

function plotReconstructionSummary(y, intervals, fitOut, labelStr)
% One figure showing:
%  - y (black)
%  - c_sumIntervals (red)
%  - c_off (red dashed)
%  - c_total (black dashed)
%  - s_oasis (green, scaled)
c_sumIntervals = fitOut.c_sumIntervals;
c_off          = fitOut.c_off;
c_total        = fitOut.c_total;
s_oasis        = fitOut.s_oasis;

figure('Name', ['Reconstruction summary: ' labelStr], 'Color', 'w');
plot(y, 'k'); hold on;
plot(c_sumIntervals, 'r', 'LineWidth', 2);
plot(c_off, 'r--', 'LineWidth', 1);
plot(c_total, 'k--', 'LineWidth', 2);

% scaled spikes for visualization
sScaled = s_oasis;
if any(sScaled)
    sScaled = sScaled / max(sScaled + eps) * (max(y) - min(y)) * 0.25 + min(y);
    plot(sScaled, 'g', 'LineWidth', 1);
end

shadeIntervals(gca, intervals, [min(y), max(y)]);

legend({'trace (y)', 'c intervals', 'c off', 'c total', 'spikes (scaled)'}, ...
       'Location', 'best');
xlabel('Time (samples)');
ylabel('Amplitude');
title(['OASIS reconstruction: ' strrep(labelStr, '_', '\_')]);
grid on;
end

function shadeIntervals(ax, intervals, yl)
% Light gray shaded rectangles with dashed border (like your original)
axes(ax); %#ok<LAXES>
for i = 1:size(intervals, 1)
    xint = [intervals(i, 1), intervals(i, 2)];
    fill([xint(1) xint(2) xint(2) xint(1)], ...
         [yl(1) yl(1) yl(2) yl(2)], ...
         [0.9 0.9 0.9], ...
         'EdgeColor', [0 0 0], ...
         'LineWidth', 0.1, ...
         'LineStyle', '--', ...
         'FaceAlpha', 0.25);
end
end

function [integrals, max_peaks, out] = oasisFIT(y, intervals, doPlot)
% Fit full trace with OASIS and quantify per-interval reconstructed components.
%
% Returns 'out' for downstream plotting:
%   out.s_oasis, out.c_oasis, out.c_recDIV, out.c_sumIntervals, out.c_off, out.c_total

arguments
    y (:,1) double
    intervals (:,2) double
    doPlot (1,1) logical = false
end

T = numel(y);

% Clamp intervals to valid bounds
intervals = round(intervals);
intervals(:,1) = max(intervals(:,1), 1);
intervals(:,2) = min(intervals(:,2), T);
if any(intervals(:,2) < intervals(:,1))
    error('One or more intervals are invalid after clamping.');
end

% OASIS deconvolution (kernel-window mode as in your original)
window = 150;
[c_oasis, s_oasis, options] = deconvolveCa(y, 'kernel', 'window', window);

% Estimate AR time constants and build kernel h
options.sn = GetSn(y);
g = estimate_time_constant(y, 2, options.sn);

w = window;

if numel(g) == 1
    h = exp(log(g) * (0:(w-1)));
elseif numel(g) == 2
    rr = roots([1, -g(1), -g(2)]);
    d = max(rr);
    r = min(rr);
    h = (exp(log(d)*(1:w)) - exp(log(r)*(1:w))) / (d - r);
else
    h = g(:);
end
h = h(:);

% Convolution matrix K
K = zeros(w);
[u, t_] = meshgrid(1:w, 1:w);
ind = 1 + t_ - u;
K(ind > 0) = h(ind(ind > 0));

% Interval-only spike trains
nInt = size(intervals,1);
s_intDIV = zeros(nInt, T);
for i = 1:nInt
    a = intervals(i,1);
    b = intervals(i,2);
    s_intDIV(i, a:b) = s_oasis(a:b);
end

% Spike train for "inside any interval"
s_int = sum(s_intDIV, 1)';          % column
s_off = s_oasis - s_int;            % column

% Reconstruct interval components
c_recDIV = zeros(nInt, T);
shift = 1;

for i = 1:nInt
    s = s_intDIV(i, :)';
    c = zeros(T, 1);

    t = 1;
    while t <= T - w + 1
        idx = t:(t + w - 1);
        c(idx) = c(idx) + K(:, 1:shift) * s(t:(t + shift - 1));
        t = t + shift;
    end

    tail_w = T - t + 1;
    if tail_w > 0
        c(t:T) = c(t:T) + K((w - tail_w + 1):w, (w - tail_w + 1):w) * s(t:T);
    end

    c_recDIV(i, :) = c;
end

% Reconstruct "off-interval" component
c_off = zeros(T, 1);
t = 1;
while t <= T - w + 1
    idx = t:(t + w - 1);
    c_off(idx) = c_off(idx) + K(:, 1:shift) * s_off(t:(t + shift - 1));
    t = t + shift;
end
tail_w = T - t + 1;
if tail_w > 0
    c_off(t:T) = c_off(t:T) + K((w - tail_w + 1):w, (w - tail_w + 1):w) * s_off(t:T);
end

% Quantify per interval
integrals = sum(c_recDIV, 2)';        % [1 x nInt]
max_peaks = max(c_recDIV, [], 2)';    % [1 x nInt]

% Collect outputs for plotting
c_sumIntervals = sum(c_recDIV, 1)';   % column
c_total = c_sumIntervals + c_off;

out = struct();
out.s_oasis        = s_oasis(:);
out.c_oasis        = c_oasis(:);
out.c_recDIV       = c_recDIV;
out.c_sumIntervals = c_sumIntervals;
out.c_off          = c_off;
out.c_total        = c_total;

% Optional diagnostic figure (close to your original two-panel style)
if doPlot
    figure('Name', 'OASIS diagnostics', 'Color', 'w');

    subplot(1,2,1);
    plot(y, 'k'); hold on;
    plot(c_sumIntervals, 'r', 'LineWidth', 2);
    plot(c_off, 'r', 'LineWidth', 1);
    plot(c_total, 'k--', 'LineWidth', 2);

    sScaled = s_oasis;
    if any(sScaled)
        sScaled = sScaled / max(sScaled + eps) * (max(y) - min(y)) * 0.25 + min(y);
        plot(sScaled, 'g', 'LineWidth', 1);
    end

    shadeIntervals(gca, intervals, [min(y), max(y)]);
    legend({'trace', 'c intervals', 'c off', 'c total', 'spikes (scaled)'}, 'Location', 'best');
    title('Full trace reconstruction');
    grid on;

    subplot(1,2,2);
    plot(y, 'k'); hold on;
    plot(c_recDIV', 'LineWidth', 1.5);
    if any(s_oasis)
        plot(sScaled, 'g', 'LineWidth', 1);
    end
    shadeIntervals(gca, intervals, [min(y), max(y)]);
    title('Per-interval reconstructed components');
    grid on;
end

end
