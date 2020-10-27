function dist = distance(aseq_loc, proteinLocation)
% This function calculates the Hamming distance between the reference
% sequence and another sequence, for all proteins.
% dist = a column vector with distances of all proteins.
    nProt = size(proteinLocation, 1);
    dist = zeros(nProt, 1);
    
    for i = 1:nProt
        startP = proteinLocation(i, 1);
        endP = proteinLocation(i, 2);
                
        % Count number of mutations within protein location boundaries
        dist(i) = sum((aseq_loc >= startP) & (aseq_loc <= endP));
        
    end
end