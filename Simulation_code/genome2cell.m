function genomeRefCell = genome2cell(genomeData)
    %----------------------------------------------------------------------
    % Remove the non-coding parts from the reference sequence.
    
    % INPUT
    % -------
    % genomeData: char array
    % Full SARS-CoV-2 genome sequence.
    
    % INPUT
    % -------
    % genomeRefCell: cell array
    % Coding region of the SARS-CoV-2 genome. Each index corresponds to a
    % protein.
    %----------------------------------------------------------------------
    
    genomeRefCell{1} = genomeData(266:805);      % NSP1
    genomeRefCell{2} = genomeData(806:2719);     % NSP2
    genomeRefCell{3} = genomeData(2720:8554);    % NSP3
    genomeRefCell{4} = genomeData(8555:10054);   % NSP4
    genomeRefCell{5} = genomeData(10055:10972);  % NSP5
    genomeRefCell{6} = genomeData(10973:11842);  % NSP6
    genomeRefCell{7} = genomeData(11843:12091);  % NSP7
    genomeRefCell{8} = genomeData(12092:12685);  % NSP8
	genomeRefCell{9} = genomeData(12686:13024);  % NSP9
    genomeRefCell{10} = genomeData(13025:13441); % NSP10
    genomeRefCell{11} = [genomeData(13442:13467), genomeData(13467:16236)]; % NSP12
    genomeRefCell{12} = genomeData(16237:18039); % NSP13
    genomeRefCell{13} = genomeData(18040:19620); % NSP14
    genomeRefCell{14} = genomeData(19621:20658); % NSP15
    genomeRefCell{15} = genomeData(20659:21552); % NSP16
    genomeRefCell{16} = genomeData(21563:25381); % Spike
	genomeRefCell{17} = genomeData(25393:26217); % NS3
    genomeRefCell{18} = genomeData(26245:26469); % E
    genomeRefCell{19} = genomeData(26523:27188); % M
    genomeRefCell{20} = genomeData(27202:27384); % NS6
    genomeRefCell{21} = genomeData(27394:27756); % NS7a
    genomeRefCell{22} = genomeData(27756:27884); % NS7b
    genomeRefCell{23} = genomeData(27894:28256); % NS8
    genomeRefCell{24} = genomeData(28274:29530); % N
end