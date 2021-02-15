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
%
%% -----------------------------PARAMETERS---------------------------------

clear all
close all

downloadDate = '20210120';
fastaFile = 'allprot0119.fasta'; % path to FASTA file downloaded from gisaid.org
refSeqID = 'EPI_ISL_402124';

%% -----------------------------START CODE---------------------------------

%% Load the GISAID data
addpath('../Data')
tic
data = fastaread(fastaFile);
disp(['Data loaded after ',num2str(toc),' seconds.'])

%% Store the sequence information and sequences in 'seqs' datastructure, 
% and store the sequence information in 'alignments' datastructure.
tic
[seqs, alignments] = store_seqs_in_datastructure(data);
disp(['Sequences are stored, this took ',num2str(toc),' seconds.'])

%% Get reference sequence (should be the 1st entry in seqs, but first check this)
if ~strcmp(seqs(1).LocusName, refSeqID)
    error('Fatal error: The first entry in seqs is not the reference sequence. Please check where the reference EPI_ISL_402124 is stored.')
end
refSeq = seqs(1);
save('refSeqProtein.mat','refSeq','-v7.3');

%% Protein alignments
tic
[excludedProteins, alignments] = perform_nw_alignments(seqs, alignments, refSeq);
disp(['Alignments ready after ',num2str(toc),' seconds.'])

%% Save results
tic
outputFile = ['proteinAlignments', downloadDate, '.mat'];
save(outputFile,'alignments', 'excludedProteins','refSeq','-v7.3')
disp(['Nominal save, elapsed time = ',num2str(toc),' seconds.'])

%% ---------------------------------FUNCTIONS------------------------------

function [refSeq, refProteinNames, refProteinLengths] = load_ref_sequence(refData)

%   -----------------------------------------------------------------------
%   This function extracts the reference sequence from datastructure
%   refData, and returns the sequence, protein names and lengths of the
%   proteins.
%   -----------------------------------------------------------------------

    refSeq = refData.proteinref;
    refProteinNames = fields(refSeq.Translation);
    % Extract reference sequence lengths:
    refProteinLengths = zeros(length(refProteinNames), 1);
    for p = 1:length(refProteinNames)
        protein = refProteinNames{p};
        refProteinLengths(p) = length(refSeq.Translation.(protein).Sequence);
        disp([protein, ' - ', num2str(refProteinLengths(p))])
    end
end

%%

function [seqs, alignments] = store_seqs_in_datastructure(data)

%   -----------------------------------------------------------------------
%   This function takes as input the 'data' variable, which is a structure
%   of all individual (protein) sequences.
%   It converts it into the 'seqs' datastructure array, which has the same
%   length as the number of submission IDs.
%   It stores the protein sequences into
%   seqs(i).Translation.proteinName.Sequence.
%   It also intializes an alignments datastructure array, that will later
%   be used to store the alignments in.
%   -----------------------------------------------------------------------

    seqs = struct();
    alignments = struct();
    v = 0;
    oldID = 'noIDwillbethesameasthisone';
    
    for i=1:length(data)
        header = split(data(i).Header, '|');
        newID = header{4};

        % store all unique protein names
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
            seqs(v).('Country') = header{end};
            alignments(v).('Country') = header{end};
            seqs(v).Translation.(proteinName).Sequence = data(i).Sequence;
            oldID = newID;
            
            % Display progression
            if mod(v, 1e3) == 0
                disp(['Sequence nr ', num2str(v)])
            end
        else
            seqs(v).Translation.(proteinName).Sequence = data(i).Sequence;
        end    
        
    end
end

%%
function [excludedProteins, alignments] = perform_nw_alignments(seqs, alignments, refSeq)

%   -----------------------------------------------------------------------
%   This function performs a Needleman-Wunsch alignment of all protein
%   sequences in seqs with the reference sequence of that protein.
%   Alignments are stored in the alignments datastructure.
%   -----------------------------------------------------------------------

    excludedProteins = zeros(1, length(seqs));
    proteinNames = fields(refSeq.Translation);
    disp(proteinNames)
    
    disp('Starting the alignments')
    for i = 1:length(seqs)
        proteinList = fields(seqs(i).Translation);
        f = 0;
        for p = 1:length(proteinList)
            protein = proteinList{p};
            

            if ismember(protein, proteinNames)

                reference = refSeq.Translation.(protein).Sequence;
                query = seqs(i).Translation.(protein).Sequence;

                Lr = length(reference);
                Lq = length(query);
                %disp(['nr ', num2str(i), ', protein ',protein])

                if length(strfind(query, 'X')) / Lq < 0.05 && Lr == Lq

                    [~, alignment] = nwalign(reference,query,'Glocal',true);
                    alignments(i).('Alignments').(protein) = alignment;
                else
                    f = f + 1;
                end
            
            end
        end
        excludedProteins(i) = f;

         % Display progression
        if mod(i, 1e3) == 0
            disp(['Alignment nr ', num2str(i)])
        end
    end
end