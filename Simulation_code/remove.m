function [seq_loc, seq_mut, seq_nMut, ...
    aseq_loc, aseq_mut, aseq_nMut, ...
    V, I, r, d, dtot, ntot, nAA, ...
    Tcells, ...
    aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
    aseqUniq_n, aseqUniq_i, aseqUniq_r] = ...
remove(i0, ...
    seq_loc, seq_mut, seq_nMut, ...
    aseq_loc, aseq_mut, aseq_nMut, ...
    V, I, r, d, dtot, ntot, nAA, ...
    Tcells, ...
    aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
    aseqUniq_n, aseqUniq_i, aseqUniq_r)

    %----------------------------------------------------------------------
    % If a viral strain i has gone extinct, this function removes it from
    % all arrays.
    %----------------------------------------------------------------------

    oldAAsequence_loc = aseq_loc{i0};
    oldAAsequence_mut = aseq_mut{i0};
    oldAAsequence_nMut = aseq_nMut(i0);
    seq_loc{i0} = []; seq_mut{i0} = [];
    aseq_loc{i0} = []; aseq_mut{i0} = []; 
    seq_nMut(i0) = 0; aseq_nMut(i0) = 0;
    
    r(i0) = 0;
    d(:,i0) = zeros(size(d,1), 1);
    dtot(i0) = 0;
    ntot = ntot - 1;
    
    % Loop through all unique AAs, if the same unique AA is still present exit,
    % the function and do not decrease AA.
    existing = find(aseqUniq_n ~= 0);
    for iUniqAA = existing 
        
        if oldAAsequence_nMut==aseqUniq_nMut(iUniqAA)
            if isequal(oldAAsequence_loc, aseqUniq_loc{iUniqAA})
                if strcmp(oldAAsequence_mut, aseqUniq_mut{iUniqAA})
                    % The old sequence is identical to this currently
                    % present unique sequence. Reduce the number which
                    % tracks how many strains have this AA by 1.
                    aseqUniq_n(iUniqAA) = aseqUniq_n(iUniqAA) - 1;
                    
                    % If the reduction by 1 causes there to be no more AA
                    % sequences of this type, the unique AA entry can be
                    % completely removed.
                    if aseqUniq_n(iUniqAA)==0
                        nAA = nAA - 1; % AA not the same as any other AA, reduce nAA   
                        Tcells(iUniqAA) = 0; 
                        Tcells_perStrain(aseqUniq_i{iUniqAA}) = 0;
                        aseqUniq_loc{iUniqAA} = [];
                        aseqUniq_mut{iUniqAA} = [];
                        aseqUniq_nMut(iUniqAA) = 0;
                        aseqUniq_i{iUniqAA} = [];
                        aseqUniq_r(iUniqAA) = 0;
                        return
                    end
                    
                    Tcells_perStrain(i0) = 0;
                    aseqUniq_i{iUniqAA}(aseqUniq_i{iUniqAA}==i0) = [];
                    
                    return
                end   
            end
        end
    end
    
end              