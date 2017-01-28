% % NuclearMorphology.m
% Code is designed to quantify nuclear and cellular shape factors. It 
% performs nuclear and cellular segmentation from confocal .lsm images, and 
% uses calcien AM segementation as a live/dead mask. 
% 
% Code from Dr. Amy Rowat's Lab
% UCLA Department of Integrative Biology and Physiology
% Los Angeles, CA
% 
% Code originally by Kendra Nyberg (December 2016)
% Code updated by Eliza King Lassman (January 2017)
%   - Added bwselect feature to manually verify segmentation
% 
% Image format: 
% Original filename format: '#.lsm'
% LSM file format from Zen Software
% Images captured on Zeiss 780 confocal
% LSM files have data in two color channels corresponding to:
%   - R channel: Hoechst stain (nuclear)
%   - G channel: Calcein AM (cellular/live stain)
%     LSM files also contain images captured at 4 different locations in the
%     wells.
% 
% Output:
%   - creates output data folder in image parent directory
%   - writes segmented image files
%   - saves .mat file containing output data cell:
% 	- Row 1: cellular data
% 	- Row 2: nuclear data
% 	- Column Indx: corresponds to each image
% 	- Each cell contains a matrix of data:
% 		- Row Indx: each object
% 		- Column 1: Area
% 		- Column 2: Perimeter

%%%% Requires Image Processing Toolbox %%%%

%%
clc; clear all; close all

% Determines the version of MATLAB. This code calls the regionprops
% function to determine the perimeter of objects. In releases 2016+, the
% perimeter calculation can yeild a circularity greater than 1, therefore,
% the 'oldperimeter' calculation is implemented in this code.
v = version('-release');
if str2double(v(3:4)) >= 16
    perim = 'PerimeterOld';
else
    perim = 'Perimeter';
end

% Initializes file locations (could be replaced with a file-selection GUI)
img_count = 28;
folder_name = 'Y:\Eliza\161206 - Osmotic Swelling Blind Analysis\';
% folder_name = '/Volumes/Rowat Lab Data 2/Eliza/161206 - Osmotic Swelling Blind Analysis/';
stack_count = 4; % number of imaged locations per lsm file

% Creates output folder
output_folder = fullfile([folder_name, 'MATLAB Analysis ', datestr(now, 'mm-dd-YY_HH-MM')]);
if ~(exist(output_folder, 'file') == 7)
    mkdir(output_folder);
end

% Creates structure array for image segmentation
for_close = strel('disk', 1); 

% Initializes output data cell
data = cell(2, img_count);

% Determine number of positions on each image
positions = 4;

% Iterate through each lsm image
for img_indx = 1:img_count
    
    % Initializes output data matrices in parent data cell
    data{1, img_indx} = [];
    data{2, img_indx} = [];
    
    % Filename
    filename = [num2str(img_indx) '.lsm'];
    
    % Reads lsm file using the tiffread function (must be in current directory or path)
    img = tiffread([folder_name filename],1);
    resol = img.lsm.VoxelSizeX*10^6; % extracts microns per pixel resolution

    % Iterates through each location captured in the lsm file
    for stack_indx = 1:stack_count
        
        % Reads lsm image
        img = tiffread([folder_name filename], stack_indx*2-1);
        
        % Creates a binary mask for channel #1 (green)
        % Channel #1 corresponds to the Calcien AM tag for live cells
        cell_mask = im2bw(img.green, graythresh(img.green)); % thresholds
        cell_mask = bwareaopen(cell_mask, floor(10000*resol)); % removes small particles
        
        % Performs watershed segmentation
        ws_mask = imextendedmin(-bwdist(~cell_mask),2);
        water_seg = watershed(imimposemin(-bwdist(~cell_mask),ws_mask));
        cell_mask(water_seg == 0) = 0; % watershed segment
        cell_mask = imfill(cell_mask, 'holes'); % fills holes
        cell_mask = imclearborder(cell_mask); %%% Remove cells on border
        cell_mask = bwareaopen(cell_mask, 100); %removes small particles

        % Creates binary image for channel #2 (red)
        % Channel #2 corresponds to the Hoechst tag for cell nuclei
        
        %%% Nuclear Segmentation Method 1: Requires image processing toolbox
        % nucl_mask = edge(img.red, 'canny');
        % nucl_mask = imclose(nucl_mask, for_close);
        % nucl_mask = imfill(nucl_mask, 'holes'); % fills holes
        % nucl_mask = medfilt2(nucl_mask, [1 1]);
        % nucl_mask = bwareaopen(nucl_mask, floor(5000*resol), 4); % removes objects smaller than 5 um(?)
        % nucl_mask = imclearborder(nucl_mask);
        
        %%%% Nuclear Segmentation Method 2: Does not require image processing toolbox, could use
        %%%% improvement!
        nucl_mask = im2bw(img.red, 1.2*graythresh(img.red));
        nucl_mask = imfill(nucl_mask, 'holes'); % fills holes
        nucl_mask = bwareaopen(nucl_mask, floor(2000*resol)); % removes objects smaller than 5 um(?)
        nucl_mask = imclearborder(nucl_mask);
        
        % Checks for dividing cells
        [cell_labels, label_count] = bwlabel(cell_mask);
        cell_circ = regionprops(cell_labels, 'Area', 'Perimeter');
        cell_circ = [cell_circ.Area]*4*pi./[cell_circ.Perimeter].^2;
        
        % Loops through all possible cells
        for cell_indx = 1:label_count
            
            % Locates current cell
            curr_cell = ismember(cell_labels, cell_indx);
            
            % If the current cell has a "fuzzy" perimeter, then remove from the
            % analysis
            if cell_circ(cell_indx) < 0.8
                cell_mask(curr_cell) = 0;
                nucl_mask(curr_cell) = 0;
            else
                
                % Determine the cournt and shape of the current nucleus
                curr_nucl = nucl_mask.*curr_cell;
                curr_nucl_shape = regionprops(curr_nucl, 'Area', 'Perimeter');
                curr_nucl_shape = 4*pi*[curr_nucl_shape.Area]./[curr_nucl_shape.Perimeter].^2; % modifies second column to become circularity
        
                [nucl_labels, nucl_count] = bwlabel(nucl_mask.*curr_cell);
                
                % if the cell does not contain a nucleus, or if it contains
                % more than one nucleus or if the nuclear shape is out of
                % this world (Circularity less than 0.4, which may need to
                % be modified for differentiated cells), then remove it
                % from both the cell and nuclear masks.
                if  nucl_count == 0 || nucl_count > 1 || curr_nucl_shape < 0.4
                    cell_mask(curr_cell) = 0;
                    nucl_mask(curr_cell) = 0;
                end
            end
        end
        
        % clears any outliers (is this nesseary?)
        nucl_mask(~cell_mask) = 0;
        
        % Generates new images
        newimg_green = img.green;
        newimg_green(~cell_mask)=0;
        newimg_red = img.red;
        newimg_red(~cell_mask)=0;
        
        % figure(1); imshow(horzcat(img.green, newimg_green))
        % figure(2); imshow(horzcat(img.red, newimg_red, nucl_mask*255))
        
        uiwait(helpdlg('On binary image, click on nuclei to be unselected. Then press enter.')); % creates help dialogue
        
        fig1 = figure(1);
        set(fig1, 'units', 'inches', 'position', [-2 5 16 3]);
        figure(1), imshow(img.red);
        fig2 = figure(2);
        set(fig2, 'units', 'inches', 'position', [6.5 5 16 3]);
        figure(2), bw_nucl = bwselect(nucl_mask);
        close all
        
        fig3 = figure(3);
        set(fig3, 'units', 'inches', 'position', [-2 5 16 3]);
        figure(3), imshow(img.green);
        fig4 = figure(4);
        set(fig4, 'units', 'inches', 'position', [6.5 5 16 3]);
        figure(4), bw_cell = bwselect(cell_mask);
        close all
        
        % create binary images where the poorly segmented nuclei/cells are
        % deleted
        nucl_mask = nucl_mask & ~bw_nucl;
        cell_mask = cell_mask & ~bw_cell;
   
        
        % Determines cell and nuclear statistics from cell_mask and
        % nucl_mask
        cell_data = regionprops(cell_mask, 'Area', 'Perimeter');
        cell_data = [[cell_data.Area]' 4*pi*[cell_data.Area]'./[cell_data.Perimeter]'.^2]; % modifies second column to become circularity
        nucl_data = regionprops(nucl_mask, 'Area', 'Perimeter');
        nucl_data = [[nucl_data.Area]' 4*pi*[nucl_data.Area]'./[nucl_data.Perimeter]'.^2]; % modifies second column to become circularity
        
        data{1,img_indx} = [data{1,img_indx}; cell_data];
        data{2,img_indx} = [data{2,img_indx}; nucl_data];
        
        output_file = ['processed_' num2str(img_indx) '.tif'];
        output_file2 = ['deleted_' num2str(img_indx) '.tif'];
        
        imwrite(cell_mask, fullfile(output_folder, output_file), 'WriteMode', 'append', 'compression','none');
        imwrite(nucl_mask, fullfile(output_folder, output_file), 'WriteMode', 'append', 'compression','none');
        imwrite(bw_cell, fullfile(output_folder, output_file2), 'WriteMode', 'append', 'compression', 'none');
        imwrite(bw_nucl, fullfile(output_folder, output_file2), 'WriteMode', 'append', 'compression', 'none');
        
     end
   
end


save(fullfile(output_folder, 'analayzed_data.mat'), 'data');
disp('Done.')
