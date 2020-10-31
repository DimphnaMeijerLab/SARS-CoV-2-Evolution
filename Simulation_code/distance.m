function dist = distance(aseq_loc, proteinLocation)
    %----------------------------------------------------------------------
    % This function calculates the Hamming distance between the reference
    % sequence and another sequence, for all proteins.
    
    % INPUTS
    % -------
    % aseq_loc: array with length nMut (= number of mutations in strain i).
    % Contains the amino acid locations where strain i is mutated w.r.t.
    % the reference.
    
    % proteinLocation: nProtx2 array, where nProt is the number of proteins.
    % Each nth row represents the starting and ending position of the nth
    % protein in the full protein sequence.
    
    % OUTPUT
    % -------
    % dist: nProtx1 array
    % Hamming distance of strain i for all proteins.
    %----------------------------------------------------------------------

    nProt = size(proteinLocation, 1);
    dist = zeros(nProt, 1);
    
    for i = 1:nProt
        startP = proteinLocation(i, 1);
        endP = proteinLocation(i, 2);
                
        % Count number of mutations within protein location boundaries
        dist(i) = sum((aseq_loc >= startP) & (aseq_loc <= endP));
        
    end
end