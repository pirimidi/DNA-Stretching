%--------------------------------------------------------------------------
% Author: Mirko Palla, PhD.
% Date: January 9, 2015.
%
% For: Optical mapping via DNA stretching in conjuction with the Synthetic 
% Biology Platform at the Wyss Institite and Church Lab - Genetics
% Department, Harvard Medical School.
%
% Purpose: This program receives a set of EXCEL files containing statistics
% relevant to DNA and some label (Alexa-dUTP, DIG/biotin immunolabeling,
% DNA origami probe, PAINT probe) object identification, then generates
% a text file with label coordinate information on DNA objects.
%
% Requirements: Folder named, 'stats_out', containing 'StretchQuant_DNA.csv' 
% 'StretchQuant_DUTPonDNA.csv' statistics files (output of CellProfiler run). 
%
% Input arguments:      x1: x coordinate of endpoint 1
%                       y1: y coordinate of endpoint 1
%                       x2: x coordinate of endpoint 2
%                       y2: y coordinate of endpoint 2
%
%                       tolerance: % tolerance the dot to be off the line
%                       determined by the DNA object
%
% This software may be used, modified, and distributed freely, but this
% header may not be modified and must appear at the top of this file.
%--------------------------------------------------------------------------

function label_coordinates(x1, y1, x2, y2, tolerance)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           ANALYSIS STARTUP                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
disp('--> Label coordinates start');
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
disp('--> Reading in DNA statistics');

dna_length_p = xlsread('StretchQuant_DNA.csv', 'K:K');
dna_x_centroid = xlsread('StretchQuant_DNA.csv', 'V:V');
dna_y_centroid = xlsread('StretchQuant_DNA.csv', 'W:W');
dna_x_object_id = xlsread('StretchQuant_DNA.csv', 'X:X');

% Iterate through all rows (DNA objects) and create 'DNA_object' data 
% structures with pre-defined fields.
disp('--> Creating DNA object structure');
fprintf('\n');

for i = 1:length(dna_length_p)
      
    DNA_objects(i).major_axis_length = dna_length_p(i);
    DNA_objects(i).x_centroid = dna_x_centroid(i);    
    DNA_objects(i).y_centroid = dna_y_centroid(i); 
    DNA_objects(i).object_number = dna_x_object_id(i); 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       LABEL STATISTICS SECTION                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read data values in from label (on DNA) statistics file.
disp('--> Reading in label statistics');

label_length_p = xlsread('StretchQuant_DUTPonDNA.csv', 'K:K');
label_x_centroid = xlsread('StretchQuant_DUTPonDNA.csv', 'U:U');
label_y_centroid = xlsread('StretchQuant_DUTPonDNA.csv', 'V:V');
label_x_object_id = xlsread('StretchQuant_DUTPonDNA.csv', 'W:W');

% Iterate through all rows (label objects) and create 'label_object' data 
% structures with pre-defined fields.
disp('--> Creating label object structure');

for j = 1:length(label_length_p)
      
    label_objects(j).major_axis_length = label_length_p(j);
    label_objects(j).x_centroid = label_x_centroid(j);    
    label_objects(j).y_centroid = label_y_centroid(j); 
    label_objects(j).object_number = label_x_object_id(j); 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    COORDINATE IDENTIFICATION SECTION                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate alpha - DNA stetching vector angle, relative to horizontal.
alpha = asind((abs(y2 - y1) / ((x1 - x2)^2 - (y2 - y1)^2))^0.5); 

fileID = fopen('coords_pixel.txt','w');
fprintf(fileID, '%6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\r\n', 'al', 'xc', 'yc', 'do', 'ao', 'xl', 'yl', 'lo', 'counter');

% For each 'DNA_object', check the number of dots residing on them.
for k = 1:length(DNA_objects)
    
    % Get current 'DNA_object' major axis length ('al') and centroid
    % coordinates ('xc', 'yc').
    al = DNA_objects(k).major_axis_length;
    xc = DNA_objects(k).x_centroid;
    yc = DNA_objects(k).y_centroid;
      
    % Calculate theoretical y-intercept of line defined by DNA segment.
    b = yc - (tand(alpha) * xc);
      
    % Initialize label counter.
    counter = 0;
    
    % Initialize labels.
    dots = struct;
    
    for h = 1:length(label_objects)
                
        % Get current 'label_object' major axis length ('ao') and centroid
        % coordinates ('xl', 'yl').
        ao = label_objects(h).major_axis_length;
        xl = label_objects(h).x_centroid;
        yl = label_objects(h).y_centroid;
                
        % Evaluate if label is on the (DNA) line; use 5% threshold here (variable).
        b_tilda = yl - (tand(alpha) * xl);          
        
        if abs((b - b_tilda) / b) < (tolerance / 100)  
            
            % Inside the circle determined by the major axis length.
            if (xl - xc)^2 + (yl - yc)^2 < (al / 2)^2

               % Create 'dots', a data structure to hold a set of 'label_objects'. 
               counter = counter + 1;          
                              
               do = DNA_objects(k).object_number;   
               lo = label_objects(h).object_number;
               
               dots(counter).major_axis_length = ao;
               dots(counter).x_centroid = xl;
               dots(counter).y_centroid = yl;
               dots(counter).object_number = lo;

               % PRINT TO FILE.
               fprintf(fileID,'%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\r\n', [al; xc; yc; do; ao; xl; yl; lo; counter]); 

            end          
        end
        
    end
          
    % Place 'lables_object' (list of dots on the same DNA fragment) into
    % current 'DNA_object'. 
    if isfield(dots, 'object_number') == 1        
        
        % Sort dots by descending 'x_centroid' value.
        sorted_dots = nestedSortStruct(dots, 'x_centroid');       
        DNA_objects(k).label_set = sorted_dots;     
        
    end
        
    % Separate DNA objects in the output file by a newline.        
    if counter > 0
        fprintf(fileID,'\r\n');
    end
       
end

fclose(fileID);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          ANALYSIS FINALIZATION                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen('coords_micron.txt','w');

% Convert pixels to microns for final coordinates.
for m = 1:length(DNA_objects)

    % Deal only with labeled DNA objects.
    if isempty(DNA_objects(m).label_set) ~= 1
        if length(DNA_objects(m).label_set) > 1
            
            % Initialize Eucledian distance array.
            euc_dist = [];      

            % Initial label.
            XI = DNA_objects(m).label_set(1).x_centroid;
            YI = DNA_objects(m).label_set(1).y_centroid;

            % Iterate through dots specific for a single DNA fragment.
            for n = 1:length(DNA_objects(m).label_set)-1

                % Current label.
                XC = DNA_objects(m).label_set(n+1).x_centroid;
                YC = DNA_objects(m).label_set(n+1).y_centroid;

                % Calculate Eucledian distance between initial and query point.
                euc_dist(n) = conv * pdist([XI, YI; XC, YC], 'euclidean');        

            end          

            % Print distance data into file (new line).
            fprintf(fileID, '%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\r\n', euc_dist);
            fprintf(fileID,'\r\n');
            
        end        
    end      
end

% Close text file with statistic summary.
fclose(fileID);

% Navigate to working directory.
cd(work_dir);

fprintf('\n');
disp('--> Label coordinates end');
fprintf('\n');
