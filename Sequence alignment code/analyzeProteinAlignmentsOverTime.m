%--------------------------------------------------------------------------
% Analysis of protein alignments at different timesteps.
% First, make protein alignments with the script
% "makeProteinAlignments.mat".
% Then, run this script to find the hamming distance of all (valid)
% sequences, and to find the mismatchBooleans for each protein.

% Sequences with an invalid collection date are excluded from the analysis.

% INPUT FILE:
% -----------

%   proteinAlignment<date>.mat:
%   Contains a data structure with all protein alignments.

% OUTPUT FILES:
% -------------

%   hammingDist<date>OverTime.mat:
%   File with arrays of all collection dates, all Hamming distances, and
%   all countries.

%   mismatchBooleanOverTime<date>.mat:
%   File with the 'mismatchBoolean' datastructure, which contains one
%   boolean array of size (1xLp) per protein (Lp=size of protein).
%   If the boolean array at position x is 1, then the GISAID database 
%   contained at least one sequence with x-1 mutations. 

%   In all filenames, <date> refers to the date on which the sequences
%   were downloaded from the GISAID database.
%--------------------------------------------------------------------------

%clear all
%close all

downloadDate = '20210120';
startT = datetime('01-DEC-2019');
endT = datetime('21-JAN-2021');
dT = 7;   % time between subsequent timepoints in days
timePoints = [startT:dT:endT, endT];

thresholdArray = 3;

%% Load data
tic
file = ['proteinAlignments', downloadDate,'.mat'];
%load(file);
disp(['data loaded after ', num2str(toc),' seconds.'])
%%
nSeqs = length(alignments);
proteinNames = fields(alignments(1).Alignments);
% remove NSP11 (it is already in NSP12):
proteinNames = proteinNames(~ismember(proteinNames, 'NSP11'));

%% Get all dates
dates = NaT(1, length(alignments));
for s = 1:length(alignments)
    dates(s) = alignments(s).LocusModificationDate;
end

%% Loop over timepoints, update mismatches for each time bin, and get the
% mismatchboolean at each timepoint. Do this for different threshold
% levels.

for t = 1:length(thresholdArray)
    threshold = thresholdArray(t);
    disp('\n')
    disp(['STARTING THRESHOLD = ', num2str(threshold)])
    
    % Initialize arrays
    IDs = 1:length(alignments);
    mismatches = initMismatches(proteinNames, alignments, IDs);
    numTimeBins = length(timePoints) - 1;
    mismatchBooleanOverTime = struct();
    
    for i = 1:numTimeBins
        t1 = timePoints(i);
        t2 = timePoints(i+1);
        disp(['Making mismatchboolean between t1 = ', datestr(t1),' and t2 =', datestr(t2)])
        IDs = find(dates >= t1 & dates < t2);

        [mismatches, ~] = findMismatches(proteinNames, alignments, mismatches, IDs);
        mismatchBoolean = fillMismatchBoolean(proteinNames, alignments, mismatches, threshold);
        % Fill in mismatchboolean in mismatchBooleanOverTime
        pNames = fields(mismatchBoolean);
        for p = 1:length(pNames)
            protein = pNames{p};
            mismatchBooleanOverTime(i).(protein) = mismatchBoolean.(protein);
        end
    end

    % Save the result
    save(['../Data/mismatchBooleanOverTime',downloadDate,'_t',num2str(threshold),'.mat'], 'proteinNames',...
          'mismatchBooleanOverTime', 'timePoints', 'dates', 'threshold', '-v7.3')
end

%% Functions 

function mismatches = initMismatches(proteinNames, alignments, IDs)
    %----------------------------------------------------------------------
    % This function finds initializes a mismatch structure with NaNs.
    %----------------------------------------------------------------------
    % INPUTS
    
    %   proteinNames: (1xnP) cell array (nP = number of proteins).
    %   Contains the names of all SARS-CoV-2 proteins.
    
    %   alignments: struct() array of length N (N = number of sequences).
    %   Contains alignments and additional information of all sequences in
    %   the GISIAID database.
    
    %   IDs: int array
    %   Contains the indices of the sequences to be incorporated in the
    %   analysis.
    
    % OUTPUT
    
    %   mismatches: struct() array of size (1xnP)
    %   For each protein, the structure array contains an NxL array, with N
    %   the number of sequences analysed, and L the length of the protein.
    %   The values in the array are NaN.
    %----------------------------------------------------------------------
    
    for p=1:length(proteinNames)
        proteinName = proteinNames{p};
        L = length(alignments(1).Alignments.(proteinName));
        mismatches.(proteinName) = NaN(length(IDs), L);
    end
end

%%
function [mismatches,nSeqPerProtein] = findMismatches(proteinNames, alignments, mismatches, IDs)
    %----------------------------------------------------------------------
    % This function finds the locations of the mismatches per protein in
    % all N alignments of interest.
    % It stores the mismatch locations in the 'mismatches' array. 
    %----------------------------------------------------------------------
    % INPUTS
    
    %   proteinNames: (1xnP) cell array (nP = number of proteins).
    %   Contains the names of all SARS-CoV-2 proteins.
    
    %   alignments: struct() array of length N (N = number of sequences).
    %   Contains alignments and additional information of all sequences in
    %   the GISIAID database.
    
    %   mismatches: struct() array of size (1xnP)
    %   Must be initialised with initMismatches. For each protein, the 
    %   structure array contains an NxL array, with N the number of 
    %   sequences analysed, and L the length of the protein.
    %   The values in the array are NaN.
    
    %   IDs: int array
    %   Contains the indices of the sequences to be incorporated in the
    %   analysis.
    
    % OUTPUT
    
    %   mismatches: struct() array of size (1xnP)
    %   For each protein, the structure array contains an NxL array, with N
    %   the number of sequences analysed, and L the length of the protein.
    %   The values in the array indicate the positions of mismatches in the
    %   nth sequence.
    
    %   nSeqPerProtein: struct() array of length nP.
    %   Contains the number of valid sequences per protein.
    %----------------------------------------------------------------------

    nP = length(proteinNames);
    nSeqPerProtein = struct();
    for p = 1:nP
        proteinName = proteinNames{p};
        nSeqPerProtein.(proteinName) = 0;
    end

    for p=1:length(proteinNames)
        proteinName = proteinNames{p};
        n = 0;
        h = 0;
        for s = IDs  
            
            n = n+1;
            alignmentStructure = alignments(s).Alignments;
            if ~isempty(alignmentStructure)
                % Check if the protein is present
                if ismember(proteinName, fields(alignmentStructure))              
                    alignment = alignmentStructure.(proteinName);
                    nSeqPerProtein.(proteinName) = nSeqPerProtein.(proteinName) + 1;
                    c = 0;
                    for i=1:size(alignment,2)
                        % check if there is an AA difference 
                        if alignment(1,i) ~= alignment(3,i) ... 
                                  && alignment(3,i) ~= 'X'                                     
                              c = c + 1;
                              % add the mismatch location to the array:
                              mismatches.(proteinName)(s,c) = i; 
                        end
                    end
                end
            else
                h = h+1;
            end
        end
    end
end

%%
function mismatchBoolean = fillMismatchBoolean(proteinNames, alignments, mismatches, threshold)
    %----------------------------------------------------------------------
    % This function creates a boolean array for each protein. If the 
    % boolean array at position x is 1, then the GISAID database contained 
    % at least one sequence with x-1 mutations. 
    %------------------------------------------
    % INPUTS
    
    %   proteinNames: cell array 
    %   Contains the names of the SARS-CoV-2 proteins.
    
    %   alignments: struct()
    %   contains alignments and additional information of all sequences in
    %   the GISIAID database.
    
    %   mismatches: struct()
    %   For each protein, the structure array contains an NxL array, with N
    %   the number of sequences analysed, and L the length of the protein.
    %   The values in the array indicate the positions of mismatches in the
    %   nth sequence.
    
    %   threshold: int
    %   Only hamming distances that are observed more than <threshold>
    %   times in the database have infuence on the logistic fits.
    
    % OUTPUTS
    
    %   mismatchBoolean: struct()
    %   A structure array filled with boolean arrays (one for each
    %   protein).
    %----------------------------------------------------------------------
    
    % loop over proteins
    for p=1:length(proteinNames)
        proteinName = proteinNames{p};
        L = length(alignments(1).Alignments.(proteinName));
        % find where the mismatches are
        M = ~isnan(mismatches.(proteinName));
        % sum up the number of mismatches for each sequence
        N = sum(M,2);
        % find the unique number of mismatches
        uniqueNumMismatches = unique(N);
        mismatchBoolean.(proteinName) = zeros(1,L);
        % set the value of the mismatchboolean to 1 if that number of
        % mismatches appears more than <threshold> times in N. Note that
        % zero mismatches means that the first point in the array is still
        % set to 1 (hence the + 1).
        for n=1:length(uniqueNumMismatches)
            if sum(N == uniqueNumMismatches(n)) >= threshold
                l = uniqueNumMismatches(n) + 1;
                mismatchBoolean.(proteinName)(l) = 1;
            end
        end 
    end
end

