%--------------------------------------------------------------------------
% MAKE PROTEIN ALIGNMENTS from protein sequence data taken from the GISAID
% database. Protein sequences are stored in a .fasta file. This script
% reads the .fasta file, and extracts the following information:
%   1) sequence ID (LocusName);
%   2) specimen collection date (LocusModificationDate);
%   3) specimen collection country (Country);
%   4) amino acid sequence stored per protein (Sequence).

% INPUT FILES
%------------

%   allprot<date>.fasta:
%   fasta file with all protein sequences, directly downloaded from the
%   GISAID database. <date> refers to the download date.

%   refSeqProtein.mat:
%   contains the reference sequence (NCBI NC_045512).

% OUTPUT FILE
%------------

%   proteinAlignment<date>.mat:
%   File contains a datastructure with all protein alignments.
%   <date> refers to the date on which the sequences were downloaded.
%--------------------------------------------------------------------------

clear all
close all
addpath('../Data')

%% Load the GISAID data
tic
file = 'allprot0923.fasta'; % path to FASTA file downloaded from gisaid.org 
data = fastaread(file);  
disp(['Data loaded after ',num2str(toc),' seconds.'])
%% Load the reference protein sequence (NCBI NC_045512)
refData = load('refSeqProtein.mat');
refSeq = refData.proteinref;
proteinNames = fields(refSeq.Translation);
% Extract reference sequence lengths:
L = zeros(length(proteinNames), 1);
for p = 1:length(proteinNames)
    protein = proteinNames{p};
    L(p) = length(refSeq.Translation.(protein).Sequence);
    disp([protein, ' - ', num2str(L(p))])
end

%% Store the sequence information and sequences in 'seqs' datastructure, 
% and store the sequence information in 'alignments' datastructure.

seqs = struct();
alignments = struct();
v = 0;
oldID = 'noIDwillbethesameasthisone';
tic
for i=1:length(data)
    header = split(data(i).Header, '|');
    newID = header{4};
    proteinName = header{1};
    if strcmp(oldID, newID) == false
        v = v+1;
        seqs(v).('LocusName') = header{4};
        alignments(v).('LocusName') = header{4};
        date = header{3};
        if strcmp(date(8:end),'-00')
            seqs(v).('LocusModificationDate') = NaT;
            alignments(v).('LocusModificationDate') = NaT;
        else
            seqs(v).('LocusModificationDate') = datetime(date);
            alignments(v).('LocusModificationDate') = datetime(date);
        end
        info = split(header{2}, '/');
        seqs(v).('Country') = header{end};
        alignments(v).('Country') = header{end};
        seqs(v).Translation.(proteinName).Sequence = data(i).Sequence;
        oldID = newID;
    else
        seqs(v).Translation.(proteinName).Sequence = data(i).Sequence;
    end    
end
disp(['Sequences are stored, this took ',num2str(toc),' seconds.'])


%% Plot coverage
coverage = NaN(length(alignments), length(proteinNames));
seqSize = NaN(length(alignments), length(proteinNames));

for i = 1:length(alignments)

    for p = 1:length(proteinNames)
        protein = proteinNames{p};
        seqStruct = seqs(i).Translation;
        if ismember(protein, fields(seqStruct))
            query = seqStruct.(protein).Sequence;
            Lq = length(query);
         	seqSize(i,p) = Lq;
            coverage(i,p) = length(strfind(query, 'X')) / Lq ;
        end
        
    end
end

sameLength = (seqSize == L');
compLength = (seqSize > 0.9*L' & seqSize < 1.1*L');

figure()
subplot(3,1,1)
histogram(coverage(sameLength))
title('Coverage of sequences with same length as reference')
set(gca, 'YScale', 'log')
xlim([0,1])
drawnow

subplot(3,1,2)
histogram(coverage(compLength))
title('Coverage of sequences with length reference ± 10%')
set(gca, 'YScale', 'log')
xlim([0,1])
drawnow

subplot(3,1,3)
histogram(coverage(~sameLength))
title('Coverage of sequences with different length from reference')
set(gca, 'YScale', 'log')
xlim([0,1])
drawnow

% Plot sizes of sequences
figure()
for p = 1:length(proteinNames)
    protein = proteinNames{p};
    subplot(5,5,p)
    histogram(seqSize(:, p))
    xline(L(p), 'r')
    title(protein)
    set(gca, 'YScale', 'log')
    drawnow
end

% Protein alignments
tic

excludedProteins = zeros(1, length(seqs));

disp('Starting the alignments')
for i = 1:length(seqs)
    proteinList = fields(seqs(i).Translation);
    for p = 1:length(proteinList)
        f = 0;
        protein = proteinList{p};
        
        if ismember(protein, proteinNames)
        
            reference = refSeq.Translation.(protein).Sequence;
            query = seqs(i).Translation.(protein).Sequence;

            Lr = length(reference);
            Lq = length(query);
            
            if length(strfind(query, 'X')) / Lq < 0.05 && Lr == Lq

                [score, alignment] = nwalign(reference,query,'Glocal',true);
                seqs(i).Translation.(protein).Alignment = alignment;
                alignments(i).('Alignments').(protein) = alignment;

            else
                f = f + 1;
            end
            
        else
            disp(['WARNING: Invalid protein in strain ', num2str(i)])
        end
    end
    excludedProteins(i) = f;
    
     % Display progression
    if mod(i, 1e3) == 0
        disp(['Alignment nr ', num2str(i)])
    end
end
disp(['Alignments ready after ',num2str(toc),' seconds.'])
% Save results
tic
save('proteinAlignment0923.mat','alignments', 'excludedProteins','-v7.3')
disp(['Nominal save, elapsed time = ',num2str(toc),' seconds.'])
