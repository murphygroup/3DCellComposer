%%%%--------ONLY WHEN USED FOR FIRST TIME --------------
% 1. go to folder ../core/mscripts
% 2. type to the command line: makeMex
%%%%----------------------------------------------------

function[] = run_CellX(file_name)

try
% clear; clc; close all;

% paths for core-code
thisFile = mfilename('fullpath');
[folder, name] = fileparts(thisFile);
cd(folder);
addpath('core/mscripts');
addpath('core/mclasses');
addpath('core/mfunctions');
addpath('core/mex');
addpath('core/mex/maxflow');
% path for code of other algos
addpath(genpath('examples_singleImages/extraFunctions'));

    
% set paths
rawImgFolder = [file_name];
resultFolder = rawImgFolder;


fprintf('running file %s \n', file_name)
    
% load files
% inputs: the images
imgFileName = ['membrane.tif'];
imgSegFileNameN = [rawImgFolder filesep imgFileName];
segImage = imgSegFileNameN;
    

% set calibration file
calibrationfilename = append(file_name, '/CellX_config.xml')

config = CellXConfiguration.readXML(calibrationfilename);

    
config.setDebugLevel(1);

config.check();
    
% get file set
frameNumber = 1;    
fileSet = CellXFileSet(frameNumber, segImage);
fileSet.setResultsDirectory(resultFolder);
    
% Run segmentation
seg = CellXSegmenter(config, fileSet);
seg.run();
segmentedCells =seg.getDetectedCells();
    
    
%------SAVE RESULTS----------
% write final images
writeSegmImages(config, fileSet, seg, segmentedCells)

catch
	exit()
end
end
