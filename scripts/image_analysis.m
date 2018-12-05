%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: September 5, 2014.
%
% For: Optical mapping via DNA stretching in conjuction with the Synthetic 
% Biology Platform at the Wyss Institite and Church Lab - Genetics
% Department, Harvard Medical School.
%
% Purpose: This program receives a set of EXCEL files containing statistics
% relevant to DNA, Alexa-dUTP object identification, then generates
% histograms, 2D/3D plots to display various correlations.
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function image_analysis()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           ANALYSIS STARTUP                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
disp('--> Image analysis start');
fprintf('\n');

% Set default number formatting.
format short;

% Conversion factor from pixel to microns, i.e., 1 pixel = 0.13 um.
conv = 0.13;

% Define current working directory.
work_dir = pwd;

% Navigate to raw data directory.
if ~exist('stats_out', 'dir')
  mkdir('stats_out');
end

cd 'stats_out';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         DNA STATISTICS SECTION                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read data values in from DNA statistics file.
dna_length_p = xlsread('StretchQuant_DNA.csv', 'K:K');
dna_length_m = conv * dna_length_p;

% Mean value and standard deviation of DNA length distribution.
mean_dna = mean(dna_length_m);
std_dna = std(dna_length_m);

% Create histogram of DNA length distribution. 
figure(1);
hist(dna_length_m, 100);
title('DNA Length Distribution');

text(3, 20, ['mean = ' num2str(mean_dna,'%.2f') '   ',...
             'sigma = ' num2str(std_dna,'%.2f')], ...
             'HorizontalAlignment', 'left',... 
             'BackgroundColor', [.7 .9 .7],...
             'EdgeColor', 'red');
          
xlabel('DNA Length (um)');
ylabel('Number of Molecules (#)');
f1 = ezfit('gauss');
showfit(f1, 'fitcolor', 'red', 'fitlinewidth', 2);
savefig('dna_length.fig');
print('-dbmp', 'dna_length.bmp'); 
fprintf('\n');
disp('--> DNA length historgam created');
fprintf('\n');

% Read data values in from DNA statistics file.
dna_comp = xlsread('StretchQuant_DNA.csv', 'F:F');

% Mean value and standard deviation of DNA compactness.
mean_comp = mean(dna_comp);
std_comp = std(dna_comp);

% Create histogram of DNA compactness. 
figure(2);
hist(dna_comp, 100);
title('DNA Compactness');

text(3, 20, ['mean = ' num2str(mean_comp,'%.2f') '   ',...
             'sigma = ' num2str(std_comp,'%.2f')], ...
             'HorizontalAlignment', 'left',... 
             'BackgroundColor', [.7 .9 .7],...
             'EdgeColor', 'red');

xlabel('Compacness Value (a.u.)');
ylabel('Number of Occurences (#)');
f2 = ezfit('gauss');
showfit(f2, 'fitcolor', 'red', 'fitlinewidth', 2);
savefig('dna_comp.fig');
print('-dbmp', 'dna_comp.bmp');
fprintf('\n');
disp('--> DNA compactness historgram created');

% Read data values in from DNA statistics file.
corr = xlsread('StretchQuant_DNA.csv', 'U:U');

% Create histogram of dUTP-on-DNA Pearson's correlation. 
figure(3);
hist(corr, 100);
title('dUTP-on-DNA Pearson Correlation');
xlabel('Pearson Correlation Coefficient');
ylabel('Number of Occurences (#)');
xlim([0 1]);
savefig('dna-dutp_corr.fig');
print('-dbmp', 'dna-dutp_corr.bmp');
disp('--> dUTP-on-DNA correlation histogram created');

% Create scatter plot of DNA compactness vs. length. 
figure(4);
scatter(dna_length_m, dna_comp);
title('DNA Compactness vs. Length');
xlabel('DNA Length (um)');
ylabel('Compacness Value (a.u.)');
savefig('dna_comp-length.fig');
print('-dbmp', 'dna-comp-length.bmp'); 
disp('--> DNA compactness vs. length scatter plot created');

% Create 2D scatter plot of dUTP-on-DNA Pearson's correlation coefficient vs. DNA length. 
figure(5);
scatter(dna_length_m, corr);
title('dUTP-on-DNA Pearson Correlation Coefficient vs. DNA Length');
xlabel('DNA Length (um)');
ylabel('Pearson Correlation Coefficient');
savefig('dna_corr-length.fig');
print('-dbmp', 'dna-corr-length.bmp'); 
disp('--> dUTP-on-DNA correlation vs. DNA length 2D scatter plot created');

% Create 2D scatter plot of dUTP-on-DNA Pearson's correlation coefficient vs. DNA compactness. 
figure(6);
scatter(dna_comp, corr);
title('dUTP-on-DNA Pearson Correlation Coefficient vs. DNA compactness');
xlabel('Compacness Value (a.u.)');
ylabel('Pearson Correlation Coefficient');
savefig('dna_corr-comp.fig');
print('-dbmp', 'dna-corr-comp.bmp');
disp('--> dUTP-on-DNA correlation vs. DNA compactness 2D scatter plot created');

% Create 3D scatter plot of DNA compactness vs. DNA length vs. dUTP-on-DNA Pearson's correlation coefficient. 
figure(7);
scatter3(dna_length_m, dna_comp, corr);
title('DNA Compactness vs. DNA Length vs. dUTP-on-DNA Pearson Correlation Coefficient');
xlabel('DNA Length (um)');
ylabel('Compacness Value (a.u.)');
zlabel('Pearson Correlation Coefficient');
savefig('dna_comp-length-corr.fig');
print('-dbmp', 'dna-comp-length-corr.bmp');
disp('--> DNA compactness vs. DNA length vs. dUTP-on-DNA correlation 3D scatter plot created');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        dUTP STATISTICS SECTION                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read data values in from dUTP (on DNA) statistics file.
dutp_length_p = xlsread('StretchQuant_DUTPonDNA.csv', 'K:K');
dutp_length_m = conv * dutp_length_p;

% Create histogram of dUTP length. 
figure(8);
hist(dutp_length_m, 100);
title('dUTP Length Distribution');
xlabel('dUTP length (um)');
ylabel('Number of Molecules (#)');
%f3 = ezfit('gauss');
%showfit(f3, 'fitcolor', 'red', 'fitlinewidth', 2);
savefig('dutp_length.fig');
print('-dbmp', 'dutp_length.bmp');
disp('--> dUTP length historgram created');

% Read data values in from dUTP statistics file.
dutp_comp = xlsread('StretchQuant_DUTPonDNA.csv', 'F:F');

% Create histogram of dUTP compactness. 
figure(9);
hist(dutp_comp, 100);
title('dUTP Compactness');
xlabel('Compacness Value (a.u.)');
ylabel('Number of Occurences (#)');
%f4 = ezfit('gauss');
%showfit(f4, 'fitcolor', 'red', 'fitlinewidth', 2);
savefig('dutp_comp.fig');
print('-dbmp', 'dutp_comp.bmp') 
disp('--> dUTP compactness histogram created');

% Create scatter plot of dUTP compactness vs. length. 
figure(10);
scatter(dutp_length_m, dutp_comp);
title('dUTP Compactness vs. Length');
xlabel('dUTP Length (um)');
ylabel('Compacness Value (a.u.)');
savefig('dutp_comp-length.fig');
print('-dbmp', 'dutp-comp-length.bmp');
disp('--> dUTP compactness vs. length 2D scatter plot created');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       OVERALL STATISTICS SECTION                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total number of single molecule DNA fragments (DFs) in the image.
dna_count = length(dna_length_m);

% Total number of fluorescently-labeled nucleotides (FNs) in the image.
spot = xlsread('StretchQuant_Spots.csv', 'B:B');
spot_count = length(spot);

% Total number of FNs on DFs in the image.
dutp_count = length(dutp_length_m);

% Labeling efficiency, i.e., label-to-baclground ratio.
lb_eff = dutp_count / spot_count * 100;

% Height and width of analysed image in pixels and microns.
image_h_p = xlsread('StretchQuant_Image.csv', 'V2:V2');
image_w_p = xlsread('StretchQuant_Image.csv', 'GL2:GL2');
image_h_m = conv * image_h_p;
image_w_m = conv * image_w_p;

% Determine image area in p^2 and in um^2.
image_area_p = image_h_p * image_w_p;
image_area_m = image_h_m * image_w_m;

% Average DFs/FNs/FN-on-DF per unit area (um^2).
unit_dna = dna_count / image_area_m;
unit_spot = spot_count / image_area_m;
unit_dutp = dutp_count / image_area_m;

% Mean value and standard deviation of dUTP length distribution.
mean_dutp = mean(dutp_length_m);
std_dutp = std(dutp_length_m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          ANALYSIS FINALIZATION                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create text file with statistic summary.
var = {'dna_count'; 'spot_count'; 'dutp_count'; 'labeling_efficiency';
       'mean_dna'; 'std_dna'; 'mean_dutp'; 'std_dutp';
       'unit_dna'; 'unit_spot'; 'unit_dutp';
       'image_heigth_pixel'; 'image_width_pixel'; 'image_area_pixel';
       'image_heigth_micron'; 'image_width_micron'; 'image_area_micron'};
   
val = [dna_count; spot_count; dutp_count; lb_eff;
       mean_dna; std_dna; mean_dutp; std_dutp;
       unit_dna; unit_spot; unit_dutp;
       image_h_p; image_w_p; image_area_p;
       image_h_m; image_w_m; image_area_m];

T = table(var, val);
writetable(T, 'output_data.txt', 'Delimiter', '\t', 'WriteRowNames', false);

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Image analysis end');
fprintf('\n');

fprintf('\n');
