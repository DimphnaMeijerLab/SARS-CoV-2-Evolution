clear all
close all

refSeqID = 'EPI_ISL_402124';
ORFs = {'ORF1ab', 'ORF3a', 'ORF3b', 'ORF6', 'ORF7a', 'ORF7b', ...
        'ORF8', 'ORF9a', 'ORF9b', 'S', 'E', 'M', 'N'};
    
pNames = {'NSP1', 'NSP2', 'NSP3', 'NSP4', 'NSP5', 'NSP6', 'NSP7', ...
          'NSP8', 'NSP9', 'NSP10', 'NSP12', 'NSP13', 'NSP14', 'NSP15', ...
          'NSP16', 'Spike', 'NS3', 'E', 'M', 'NS6', 'NS7a', 'NS7b', ...
          'NS8', 'N', 'NS9b', 'NS9c'};
N = length(pNames); 

ntseq = struct();
for i = 1:length(ORFs)
    orf = ORFs{i};
    filename = [refSeqID, '-', orf, '.fasta'];
    fastafile = fastaread(filename);
    ntseq.(orf) = fastafile.Sequence;
end

gRefCell{1} = ntseq.ORF1ab(1:540);              % NSP1
gRefCell{2} = ntseq.ORF1ab(541:2454);           % NSP2
gRefCell{3} = ntseq.ORF1ab(2455:8289);          % NSP3
gRefCell{4} = ntseq.ORF1ab(8290:9789);          % NSP4
gRefCell{5} = ntseq.ORF1ab(9790:10707);         % NSP5
gRefCell{6} = ntseq.ORF1ab(10708:11577);        % NSP6
gRefCell{7} = ntseq.ORF1ab(11578:11826);        % NSP7
gRefCell{8} = ntseq.ORF1ab(11827:12420);        % NSP8
gRefCell{9} = ntseq.ORF1ab(12421:12759);        % NSP9
gRefCell{10} = ntseq.ORF1ab(12760:13176);       % NSP10
gRefCell{11} = [ntseq.ORF1ab(13177:13203), ntseq.ORF1ab(13203:15971)];  % NSP12
gRefCell{12} = ntseq.ORF1ab(15972:17774);       % NSP13
gRefCell{13} = ntseq.ORF1ab(17775:19355);       % NSP14
gRefCell{14} = ntseq.ORF1ab(19356:20393);       % NSP15
gRefCell{15} = ntseq.ORF1ab(20394:21287);       % NSP16
gRefCell{16} = ntseq.S;                         % Spike
gRefCell{17} = ntseq.ORF3a;                     % NS3
gRefCell{18} = ntseq.E;                         % E
gRefCell{19} = ntseq.M;                         % M
gRefCell{20} = ntseq.ORF6;                      % NS6
gRefCell{21} = ntseq.ORF7a;                     % NS7a
gRefCell{22} = ntseq.ORF7b;                     % NS7b
gRefCell{23} = ntseq.ORF8;                      % NS8
gRefCell{24} = ntseq.N;                         % N
gRefCell{25} = ntseq.ORF9a;                     % NS9b
gRefCell{26} = ntseq.ORF9b;                     % NS9c

% Translate nt sequence
pRefCell = nt2aa(gRefCell); 

% combine the cell arrays into two character arrays
gRefSeq = horzcat(gRefCell{:});
pRefSeq = horzcat(pRefCell{:});

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

save('../Data/refSeq.mat', 'pNames', 'gRefSeq', 'pRefSeq', 'proteinLocation', 'genomeLocation')
