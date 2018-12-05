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

function stats_summary()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           ANALYSIS STARTUP                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
disp('--> Statistical summary start');
fprintf('\n');

% Set default number formatting.
format short;

% Define current working directory.
work_dir = pwd;

% Navigate to raw data directory.
if ~exist('summary', 'dir')
  mkdir('summary');
end

cd 'summary';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        SUMMARY STATISTICS SECTION                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in all summary statistic text files one-by-one.
list = dir('*.txt');

% Create data container for barplot display.
C = zeros(length(list), 2);

for i = 1:length(list)
    
  disp(['--> Processing file: ', list(i).name]);  
  T = tdfread(list(i).name);
  C(i, 1) = T.val(1);
  C(i, 2) = T.val(2);
  C(i, 3) = T.val(3);
  
  labeling_efficiency(i) = T.val(4);
  mean_dna(i) = T.val(5);
  std_dna(i)= T.val(6);
  mean_dutp(i) = T.val(7);
  std_dutp(i) = T.val(8);
  unit_dna(i) = T.val(9);
  unit_spot(i) = T.val(10);
  unit_dutp(i) = T.val(11);
  
end

% Create barplot of total DNA and dUTP-on-DNA count for all experiments.
figure(1);
bar(C, 'grouped');
title('Total DNA and dUTP-on-DNA Object Count');
xlabel('Experiment (#)');
ylabel('Number of Objects');
legend('DNA', 'dUTP spots', 'dUTP-on-DNA', 'Location', 'NorthWest');
xlim([0 length(list)+1]);
savefig('count.fig');
print('-dbmp', 'count.bmp'); 
fprintf('\n');
disp('--> Object count barplot created');

% Create 2D plot of DNA labeling efficiency for all experiments. 
figure(2);
plot(1:i, labeling_efficiency, 'b--o', 'LineWidth', 1,...
                               'MarkerSize', 10,...
                               'MarkerEdgeColor', 'b',...
                               'MarkerFaceColor', 'w');
                           
title('DNA Labeling Effciency');
xlabel('Experiment (#)');
ylabel('Efficiency (%)');
legend('DNA labeling efficiency');
xlim([0 length(list)+1]);
grid;
savefig('labeling_efficiency.fig');
print('-dbmp', 'labeling_efficiency.bmp'); 
disp('--> Labeling efficiency 2D plot created');

% Create 2D plot of mean/std of DNA/dUTP-on-DNA length for all experiments. 
figure(3);
plot(1:i, mean_dna, 'ro', 1:i, mean_dutp, 'bo',...
     1:i, std_dna, 'ks', 1:i, std_dutp, 'gs');
                           
title('Mean and Average of DNA/dUTP-on-DNA Length');
xlabel('Experiment (#)');
ylabel('Length (um)');
legend('Mean DNA', 'Mean dUTP-on-DNA', 'STD DNA', 'STD dUTP-on-DNA', 'Location', 'NorthWest');
xlim([0 length(list)+1]);
grid;
savefig('length.fig');
print('-dbmp', 'length.bmp'); 
disp('--> Mean-STD length 2D plot created');

% Create 2D plot of unit area normalized DNA/dUTP/dUTP-on-DNA count for all experiments. 
figure(4);
plot(1:i, unit_dna, 'b*', 1:i, unit_spot, 'g*', 1:i, unit_dutp, 'r*');
                           
title('Area Normalized DNA/dUTP/dUTP-on-DNA Count');
xlabel('Experiment (#)');
ylabel('Number of Objects per Square Micron');
legend('DNA', 'dUTP', 'dUTP-on-DNA', 'Location','NorthWest');
xlim([0 length(list)+1]);
grid;
savefig('area.fig');
print('-dbmp', 'area.bmp'); 
disp('--> Area normalized object count 2D plot created');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          ANALYSIS FINALIZATION                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Statistical summary end');
fprintf('\n');
