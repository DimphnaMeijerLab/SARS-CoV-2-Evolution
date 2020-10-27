function [data, params] = gillespie(nr, U0, mu_array)
        
        tic
        myStream = RandStream('mlfg6331_64', 'Seed', str2num(nr));

        %% Viral evolution with Gillespie algorithm

        % Source of model: 
        % Woo & Reifman (2013). Quantitative Modeling of Virus Evolutionary
        % Dynamics and Adaptation in Serial Passages Using Empirically Inferred
        % Fitness Landscapes. Journal of Virology, V88 N2, p1039 - 1050.

        % Possible events/reactions:
        % Infection: U + Vn -> In     (R1) rate = a, reaction_rate = a*U*Vn 
        % Replication: In -> In + Vm  (R2) rate = r(d), reaction_rate = r*I
        % Death : In -> 0             (R3) rate = b
        % Clearence: Vn -> 0          (R4) rate = b
        
        %% Get protein and genome reference sequences #############################
        [gRefSeq, pRefSeq, pNames, proteinLocation, genomeLocation] = getRefSeq();
        translateCodon = geneticcode();
        [beta, sigma] = logisticRegressionProteins();

        %% Initialize ############################################################# 
        distribution = 'normal';

        if U0 == 1e5
            r0 = 1.5 ;
            a =  4.5e-04 ;
            b =  0.9 ;
            c =  0.9 ;
            V0 = 400;
        elseif U0 == 1e4
            r0 = 1.5 ;
            a =  4.5e-3 ;
            b =  0.9 ;
            c =  0.9 ;
            V0 = 40;
        elseif U0 == 1e6
            r0 = 1.5 ;
            a =  4.5e-5 ;
            b =  0.9 ;
            c =  0.9 ;
            V0 = 4e3;
        end
        alpha = 1;

        params = struct();
        params.('U0') = U0;
        params.('r0') = r0;
        params.('a') = a;
        params.('b') = b;
        params.('c') = c;
        params.('V0') = V0;
        params.('mu') = mu_array;

        T = inf;                    % maximal time (days)
        t_anti = inf;

        antiviral = false;

        % Follows
        N_mu = length(mu_array);                            % number of mutation rates to test
        L = length(gRefSeq);                                % length of genome sequence
        La = length(pRefSeq);  
        d0 = zeros(length(pNames),1);                       % distance of WTseq to fittest strain

        data = struct();
        S = 1e3;                                            % initial size of arrays

        % Start for-loop over mutation rates ######################################
        for k = 1:N_mu
            mu = mu_array(k);                                                     % mutation rate
            % initialize
            disp(['Mutation rate: ', num2str(mu)])
            U = U0;                                                                % number of infected cells.
            seq_loc = cell(1,S); aseq_loc = cell(1,S);                             % number of infected cells.
            seq_mut = cell(1,S); aseq_mut = cell(1,S);
            seq_nMut = zeros(1,S); aseq_nMut = zeros(1,S);                         % cell array to keep track of aa sequences   

            aseqUniq_loc = cell(1,S);
            aseqUniq_mut = cell(1,S);
            aseqUniq_nMut = zeros(1,S);
            aseqUniq_n = zeros(1,S); aseqUniq_n(1) = 1;
            aseqUniq_i = cell(1,S); aseqUniq_i{1} = 1;
            aseqUniq_r = zeros(1,S); aseqUniq_r(1) = r0;

            % counters
            ntot = 1;                                                               % initial number of distinct viral genotypes
            nAA = 1;                                                                % initial number of distinct viral phenotypes
            s = 0;                                                                  % counter for while loop
            m = 1;                                                                  % counter to collect stats
            t = 0.0;                                                                % initial time (days)
            t_collect = 0.0;
            % arrays to keep track of viral strains
            V = zeros(1, S);                                                        % number of free viral particles per strain
            V(1) = V0;
            I = zeros(1, S);                                                        % number of infected cells per strain
            Tcells = zeros(1, S);
            Tcells_perStrain = zeros(1, S);
            r = zeros(1, S);                                                        % fitness per strain
            r(1) = r0;                                                              
            d = zeros(length(pNames), S);                                           % distance to reference per strain per protein
            dtot = zeros(1, S);                                                     % total distance to reference (sum of all proteins) per strain
            % arrays to collect stats
            Y = zeros(S, La + 1);                                                   % matrix for number of viruses with distance d from WT seq.
            Y(1,1) = V(1);      

            data_collect = zeros(S, 6);                                                     % matrix to collect time, # distinct genotypes, # distinct phenotypes, # not-infected cells, # infected cells, # viruses.
            data_collect(1,:) = [t, ntot, nAA, U, sum(I), sum(V)];                          
            meanDistance = zeros(1,S);                                              % mean number of mutations
            meanFitness = zeros(1,S);                                               % mean fitness of population
            meanFitness(1) = r0; 
            maxD = zeros(1,S);
            maxR = zeros(1,S);
            maxR(1) = r0;
            diversity = zeros(1,S);

            % Start while loop ####################################################
            while t < T
                 if t >= t_anti && antiviral == false
                     r0 = alpha*r0;
                     r  = alpha*r;
                     antiviral = true;
                 end    
                 s = s + 1;                                                         
                 % calculate total reaction rate
                 Rtot = (a*U + b)*sum(V) + c*sum(I) + sum(r.*I);                  
                 dt = (1/Rtot) * log(1/rand(myStream));                                       % time step is exponentially distributed and depends on total reaction rate
                 t = t + dt;
                 
                 % choose a random number between 0 and Rtot
                 rnd = Rtot * rand(myStream);
                 z = 0;
                 
                 % Increase the size of arrays if it is necessary
                 existing = find(V+I);
                 if existing(end) + 10 > length(V)
                     V = [V, zeros(1,S)];
                     I = [I, zeros(1,S)];
                     Tcells_perStrain = [Tcells_perStrain, zeros(1, S)];
                     r = [r, zeros(1,S)];
                     d = [d, zeros(length(pNames), S)];
                     dtot = [dtot, zeros(1,S)];
                     seq_loc = [seq_loc, cell(1,S)]; aseq_loc = [aseq_loc, cell(1,S)];
                     seq_mut = [seq_mut, cell(1,S)]; aseq_mut = [aseq_mut, cell(1,S)];
                     seq_nMut = [seq_nMut, zeros(1,S)]; aseq_nMut = [aseq_nMut, zeros(1,S)];
                 end
                 if existing(end) + 10 > length(Tcells)
                     Tcells = [Tcells, zeros(1,S)];
                     aseqUniq_loc = [aseqUniq_loc, cell(1,S)];
                     aseqUniq_mut = [aseqUniq_mut, cell(1,S)];
                     aseqUniq_nMut = [aseqUniq_nMut, zeros(1,S)];
                     aseqUniq_n = [aseqUniq_n, zeros(1,S)];
                     aseqUniq_i = [aseqUniq_i, cell(1,S)];
                     aseqUniq_r = [aseqUniq_r, zeros(1,S)];
                 end
                 % loop over viral strains
                 for i = existing
                     % choose which of the 4 reactions will occur 
                     % and for which strain
                     z = z + a*U*V(i);                                              % R1: infection by strain i
                     if z > rnd                                                         
                         I(i) = I(i) + 1;                                               
                         V(i) = V(i) - 1;
                         U = U - 1;      
                         break
                     end
                     z = z + r(i)*I(i);                                             % R2: replication of strain i
                     if z > rnd
                         i0 = i;
                        [seq_loc, seq_mut, seq_nMut, ...
                            aseq_loc, aseq_mut, aseq_nMut, ...
                            ntot, nAA, V, r, I, d, dtot, ...
                            aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
                            aseqUniq_n, aseqUniq_i, aseqUniq_r] = ...
                        replicate(i0, myStream, ...
                            seq_loc, seq_mut, seq_nMut, ...
                            aseq_loc, aseq_mut, aseq_nMut, ...
                            mu, ntot, nAA, V, r, I, d, dtot, ...
                            sigma, r0, ...
                            gRefSeq, L, pRefSeq, beta, proteinLocation, translateCodon, ...
                            aseqUniq_loc, aseqUniq_mut, aseqUniq_nMut, ...
                            aseqUniq_n, aseqUniq_i, aseqUniq_r, distribution);
                        break
                     end
                     z = z + c*I(i);                                                % R3: clearence of a cell infected by strain i
                     if z > rnd
                         I(i) = I(i) - 1;
                         if V(i) + I(i) == 0    % if the strain has gone extinct, remove it
                             i0 = i;
                            [seq_loc, seq_mut, seq_nMut, ...
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
                                aseqUniq_n, aseqUniq_i, aseqUniq_r);
                         end
                         break % for-loop
                     end
                     
                     z = z + b*V(i);                                                % R4: clearence of a free viral particle of strain i
                     if z > rnd
                         V(i) = V(i) - 1;
                         if V(i) + I(i) == 0    % if the strain has gone extinct, remove it
                             i0 = i;
                            [seq_loc, seq_mut, seq_nMut, ...
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
                                aseqUniq_n, aseqUniq_i, aseqUniq_r);
                         end
                         break % for-loop
                     end
                 end % for-loop 

                % Stop if there are no more viral particles left:
                 if sum(I+V) == 0                                                  	
                    break
                 end

                 % Collect statistics every 0.1 days:
                 if t - t_collect > 0.1
                     t_collect = t;
                     m = m + 1;
                     % Increase size of arrays if necessary
                     if m > length(meanFitness)
                         Y = [Y; zeros(S, La + 1)];
                         data_collect = [data_collect; zeros(S, 6)];
                         meanDistance = [meanDistance, zeros(1,S)];
                         meanFitness = [meanFitness, zeros(1,S)];
                         maxD = [maxD, zeros(1,S)];
                         maxR = [maxR, zeros(1,S)];
                         diversity = [diversity, zeros(1,S)];
                     end
                     % Group strains with the same distance and store them in Y:
                     alive = (V + I ~= 0);
                     shortV = V(alive); % remove zeros of extinct strains
                     shortI = I(alive);
                     for j = 1:max(dtot)+1
                        dist = j - 1;
                        d_logical = (dtot(alive) == dist);
                        sumV = sum(shortV(d_logical));
                        sumI = sum(shortI(d_logical));
                        Y(m,j) = sumV + sumI;
                     end

                     % Collect all data
                     viralTiter = shortV + shortI;
                     meanFitness(m) = sum(r(alive) .* viralTiter) / sum(viralTiter);
                     meanDistance(m) = sum(dtot(alive) .* viralTiter) / sum(viralTiter);
                     data_collect(m,:) = [t, ntot, nAA, U, sum(I), sum(viralTiter)];
                     maxD(m) = max(dtot);
                     maxR(m) = max(r);

                     diversity(m) = 1 - sum((viralTiter/sum(viralTiter)).^2); % Simpson's index as diversity

                     %disp(['t=',num2str(t),', \t ntot=',num2str(ntot), ', \t U=', num2str(U)])
                 end

            end % while-loop

            % Remove zeros at the end of vectors
            maxY = find(sum(Y) > 0, 1, 'last');

            Y = Y(1:m, :);         % discard the last row of zeros   
            data_collect = data_collect(1:m, :);
            meanFitness = meanFitness(1:m);
            meanDistance = meanDistance(1:m);
            maxD = maxD(1:m);
            maxR = maxR(1:m);
            diversity = diversity(1:m);

            % Calculate stationary values:
            time = data_collect(1:m, 1)';
            delta_t = time(2:end) - time(1:end-1);
            sumY = 1 ./ sum(Y,2);
            relativeY = Y .* repmat(sumY, 1, size(Y,2));             
            statY = (delta_t * relativeY(2:end,:)) / time(end);
            statD = sum(meanDistance(2:end) .* delta_t) / time(end);
            statR = sum(meanFitness(2:end) .* delta_t) / time(end); 
            statDiv = sum(diversity(2:end) .* delta_t) / time(end);

            % Collect information for each simulation
            data(k).alpha = alpha;
            data(k).Mu = mu;
            data(k).t = t;
            data(k).ntot = ntot;
            data(k).nAA = nAA;
            data(k).U_sum = sum(U);
            data(k).I_sum = sum(I);
            data(k).V_sum = sum(V);
            data(k).data_collect = data_collect;
            data(k).statY = statY(1:maxY);
            data(k).relativeY = relativeY(:, 1:maxY);
            data(k).statD = statD;
            data(k).statR = statR;
            data(k).maxD = maxD;
            data(k).maxR = maxR;
            data(k).diversity = diversity;
            data(k).statDiv = statDiv;
        end % for-loop

        disp(['Simulation time ', num2str(toc), ' seconds.'])
end