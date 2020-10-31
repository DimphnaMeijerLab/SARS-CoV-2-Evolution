function [gRefSeq, pRefSeq, pNames, proteinLocation, genomeLocation] = getRefSeq()

    %----------------------------------------------------------------------
    % Load the coding region of the SARS-CoV-2 reference sequence.
    
    % OUTPUTS
    % -------
    % gRefSeq: char array
    % Full coding SARS-CoV-2 genome sequence.
    
    % pRefSeq: char array
    % Full coding SARS-CoV-2 protein (amino acid) sequence.
    
    % pNames: cell array
    % Names of the proteins.
    
    % proteinLocation: nProtx2 array, where N is the number of proteins.
    % Each nth row represents the starting and ending position of the nth
    % protein in the full protein sequence.
    
    % genomeLocation: Nx2 array, where N is the number of proteins.
    % Each nth row represents the starting and ending position of the nth
    % protein in the full genome sequence.
    %----------------------------------------------------------------------

    % Load data
    pNames = ...
        {'NSP1', 'NSP2', 'NSP3', 'NSP4', 'NSP5', 'NSP6', 'NSP7', 'NSP8', ...
        'NSP9', 'NSP10', 'NSP12', 'NSP13', 'NSP14', 'NSP15', 'NSP16', ...
        'Spike', 'NS3', 'E', 'M', 'NS6', 'NS7a', 'NS7b', 'NS8', 'N'};
    N = length(pNames);                                                     % number of proteins
    gData = load('refSeqGenome.mat');
    gData = gData.refseq;
    
    gRefCell = genome2cell(gData);                                          % cell array with all genome sequences    
    % translate genome into protein sequence:
    pRefCell = nt2aa(gRefCell); 

    % combine the cell arrays into two character arrays
    gRefSeq = horzcat(gRefCell{:});                                       	% char array with all genome sequences 
    pRefSeq = horzcat(pRefCell{:});                                         % char array with all protein sequences
    % Collect starting and ending locations of each protein
    proteinLocation = zeros(N, 2);
    genomeLocation = proteinLocation;
    startP = 1;                                                             % start of protein in protein sequence
    startG = 1;                                                             % start of protein in genome sequence
    for i = 1:N
        L = length(pRefCell{i});                                            % length of protein
        proteinLocation(i, :) = [startP, min(startP + L -1, length(pRefSeq))];  % start and end of protein in amino acid sequence
        genomeLocation(i, :) = [startG, min(startG + 3*L, length(gRefSeq))];    % start and end of protein in nuclotide sequence
        startP = startP + L;
        startG = startG + 3*L;
    end
end