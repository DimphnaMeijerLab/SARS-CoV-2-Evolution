function gillespie_transmission(nr, U0, mu, N_ind, varargin)
        %------------------------------------------------------------------
        % Simulate a viral infection with the Gillespi algorithm.
        %
        % INPUT PARAMETERS
        %-----------------
        %
        %   nr: int
        %   Simulation number, used for random number generator.
        %
        %   U0: int
        %   Initial number of not-infected cells.
        %
        %   mu_array: double or array of doubles.
        %   Array with mutation rates you want to test. If you want to do one
        %   single simulation, mu_array is a single numbers.
        %
        % Optional input parameters:
        %---------------------------
        %
        %   a: double
        %   Infection rate. Default depends on U0.
        %
        %   b: double.
        %   Clearance rate of viral particles. Default is 0.9.
        %
        %   c: double.
        %   Death rate of infected cells. Default is 0.9.
        %
        %   r0: double.
        %   Replication rate of reference sequence. Default is 1.5.
        %
        %   lambda: double
        %   Penalization parameter of ridge regression. Default is 0 (no
        %   rigde regression).
        %
        %   distribution: string; choose from 'normal' or 'gamma' or 'empirical'.
        %   Probability distribution of the fitness of newly emerged strains.
        %   Default is the normal distribution.
        %
        %   sigma: double
        %   Standard deviation of normal / gamma distribution in DFE.
        %   Default is the error of fit of the logistic regression.
        %   
        %   fitnessFunction: function handle
        %   Function to calculate fitness. The first input argument must be b (a
        %   2x1 column vector with b0 and b1, the outputs of the logistc fit).
        %   The second input argument must be the hamming distance, d. 
        %   Default is the multiplicative fitness function.
        %
        %   wholeGenome: logical
        %   Specifies wether the logistic regression is done per protein
        %   (if wholeGenome is false) or on the entire genome. Default is
        %   false.
        %
        % OUTPUT PARAMETERS
        %-----------------
        % data: structure array with the same length as mu_array.
        % Structure to collect all data per mutation rate. The fields are:
        %
        %   data_collect:   T x 6 array containing values that change over 
        %                   time, where T is the number of timepoints in
        %                   the simulation.
        %                   COLUMN 1: time array.
        %                   COLUMN 2: the total number of distinct genome
        %                             sequences per timepoint.
        %                   COLUMN 3: the total number of distinct protein
        %                             sequences per timepoint.
        %                   COLUMN 4: the total number of non-infected
        %                             cells per timepoint (U).
        %                   COLUMN 5: the total number of infected
        %                             cells per timepoint (I).
        %                   COLUMN 6: the total number free viral particles
        %                             per timepoint (V).
        %                            
        %   relativeY:      T x maxY array, where T is the number of
        %                   timepoints and maxY is the maximal hamming
        %                   distance that arose in the simulation + 1. Each
        %                   value (t,i) in the array denotes the relative
        %                   abundance of viruses with hamming distance i-1 
        %                   at timepoint t.
        %
        %   statY:          1 x maxY array with the time-integrated
        %                   relative abundance of viruses with a hamming
        %                   distance i-1.
        %
        %
        %   statR:          double denoting the stationary reproduction
        %                   rate (i.e. reproduction rate integrated over
        %                   time).
        %
        %   maxR:           double, maximal reproduction rate that arose in
        %                   the simulation.
        %
        %   diversity:      1 X T array with de diversity of the viral
        %                   population at all T timepoints.
        %
     	%   statDiv:        double denoting the stationary diversity (i.e.
     	%                   the viral diversity integrated over time).
        %
        %   statD:          double denoting the stationary hamming distance 
        %                   (i.e. the hamming distance integrated over
        %                   time).
        %
        %   maxD:           double, maximal hamming distance that arose in
        %                   the simulation.
        %
        %------------------------------------------------------------------
        
        tic
        myStream = RandStream('mlfg6331_64', 'Seed', nr);
        outputFolder = ['Output/Output_', num2str(nr)];
        mkdir(outputFolder)

        %% Optional parameters ############################################
        
        % specify fitness function
        distribution = 'normal';
        fitnessFunction = @multiplicative_fitness;
        lambda = 0;
        
        %% Get logistic regression results ################################
        
        % specify which mismatchboolean to use to to logistic regression on
        mismatchBoolFileName = 'mismatchBoolean20210120_t3.mat';
        logisticRegressionFunction = @logisticRegressionProteins;
        if any(strcmp(varargin,'wholeGenome'))
            index = find( strcmp(varargin,'wholeGenome') );
            if varargin{index + 1} == true
                mismatchBoolFileName = ['mismatchBoolean_',...
                                        'WholeGenome_',...
                                        '20210120_t3.mat'];
                logisticRegressionFunction = @logisticRegressionWholeGenome;
            end
        end
        % load mismatchboolean and do logistic regression
        mismatchBooleanStructure = load(mismatchBoolFileName);
        mismatchBoolean = mismatchBooleanStructure.mismatchBooleanOverTime;
        [beta, sigma_proteins] = logisticRegressionFunction(mismatchBoolean, lambda);

        % specify default replication dynamics
        if U0 == 1e4
            r0 = 1.5 ;
            a =  4.5e-3 ;
            b =  0.9 ;
            c =  0.9 ;
            V0 = 40;
        elseif U0 == 1e5
            r0 = 1.5 ;
            a =  4.5e-4 ;
            b =  0.9 ;
            c =  0.9 ;
            V0 = 400;
        elseif U0 == 1e6
            r0 = 1.5 ;
            a =  4.5e-5 ;
            b =  0.9 ;
            c =  0.9 ;
            V0 = 4e3;
        end
                
        p = inputParser;
        isFunction = @(f) isa(f,'function_handle');
        addParameter(p, 'a', a, @isnumeric)
        addParameter(p, 'b', b, @isnumeric)
        addParameter(p, 'c', c, @isnumeric)
        addParameter(p, 'r0', r0, @isnumeric)
        addParameter(p, 'lambda', lambda, @isnumeric)
        addParameter(p, 'distribution', distribution, @(s)ischar(s))
        addParameter(p, 'sigma', sigma_proteins, @isnumeric)
        addParameter(p, 'fitnessFunction', fitnessFunction, isFunction)
        addParameter(p, 'wholeGenome', false, @islogical)
        
        parse(p, varargin{:})
        a = p.Results.a;
        b = p.Results.b;
        c = p.Results.c;
        r0 = p.Results.r0;
        lambda = p.Results.lambda;
        distribution = p.Results.distribution;
        sigma = p.Results.sigma;
        fitnessFunction = p.Results.fitnessFunction;
        wholeGenome = p.Results.wholeGenome;
                
        %% Get protein and genome reference sequences #####################
        refSeqStructure = load('refSeq.mat');
        pNames = refSeqStructure.pNames;
        gRefSeq = refSeqStructure.gRefSeq;
        pRefSeq = refSeqStructure.pRefSeq;
        proteinLocation = refSeqStructure.proteinLocation;
        translateCodon = geneticcode();

        %% Initialize #####################################################

        params = struct();
        params.('U0') = U0;
        params.('r0') = r0;
        params.('a') = a;
        params.('b') = b;
        params.('c') = c;
        params.('V0') = V0;
        params.('mu') = mu;
        params.('lambda') = lambda;
        params.('distribution') = distribution;
        params.('sigma') = sigma;
        params.('fitnessFunction') = func2str(fitnessFunction);

        T = inf;
        
        % Follows
        %N_mu = length(mu_array);        % number of mutation rates to test
        L = length(gRefSeq);            % length of genome sequence
        La = length(pRefSeq);           % length of protein sequence

        S = 1e3;                     	% initial size of arrays
        
        t_L = 4.6;                      % latend period (days)
        t_I = 5;                        % infectious period (days)
        
        % initialise transmission arrays
        trans_seq_loc = cell(1,S); trans_aseq_loc = cell(1,S);                            
        trans_seq_mut = cell(1,S); trans_aseq_mut = cell(1,S);
        trans_seq_nMut = zeros(1,S); trans_aseq_nMut = zeros(1,S);                         

        trans_aseqUniq_loc = cell(1,S);
        trans_aseqUniq_mut = cell(1,S);
        trans_aseqUniq_nMut = zeros(1,S);
        trans_aseqUniq_n = zeros(1,S); trans_aseqUniq_n(1) = 1;
        trans_aseqUniq_i = cell(1,S); trans_aseqUniq_i{1} = 1;
        trans_aseqUniq_r = zeros(1,S); trans_aseqUniq_r(1) = r0;

        % arrays to keep track of viral strains
        trans_V = zeros(1, S);                                                        % number of free viral particles per strain
        trans_V(1) = V0;
        trans_r = zeros(1, S);                                                        % fitness per strain
        trans_r(1) = r0;                                                              
        trans_d = zeros(length(pNames), S);                                           % distance to reference per strain per protein
        trans_dtot = zeros(1, S);                                                     % total distance to reference (sum of all proteins) per strain

        % arrays to collect stats
        trans_Y = zeros(S, La + 1);                                                   % matrix for number of viruses with distance d from WT seq.
        trans_Y(1,1) = V0;
        
        trans_Tcells = zeros(1, S);
        trans_Tcells_perStrain = zeros(1, S);

        %% Start for-loop over mutation rates #############################
        for k = 1:N_ind
                        
            fprintf('\n Iteration %d: U0=%d \t V0=%d \t a=%f \t b=%f \t c=%f \t r0=%f \t mu=%e \n', k, U0, V0, a, b, c, r0, mu)
         	
            % initialize
            U = U0;                                                                     % initial number of uninfected cells.
            t = 0.0;                                                                    % initial time (days)
            t_collect = t;
            
            % arrays to keep track of viral strains
            I = zeros(1, S);                                                        % number of infected cells per strain            
            V = trans_V;                                                        % number of free viral particles per strain
            Y = trans_Y;
            r = trans_r;                                                        % fitness per strain
            d = trans_d;                                           % distance to reference per strain per protein
            dtot = trans_dtot; 
            
            Tcells = trans_Tcells;
            Tcells_perStrain = trans_Tcells_perStrain;
                                                   
            seq_loc = trans_seq_loc; 
            aseq_loc = trans_aseq_loc;                            
            seq_mut = trans_seq_mut; 
            aseq_mut = trans_aseq_mut;
            seq_nMut = trans_seq_nMut; 
            aseq_nMut = trans_aseq_nMut;                         

            aseqUniq_loc = trans_aseqUniq_loc;
            aseqUniq_mut = trans_aseqUniq_mut;
            aseqUniq_nMut = trans_aseqUniq_nMut;
            aseqUniq_n = trans_aseqUniq_n;
            aseqUniq_i = trans_aseqUniq_i;
            aseqUniq_r = trans_aseqUniq_r;

            % counters
            ntot = sum(V + I > 0);                                                	% initial number of distinct viral genotypes
            nAA = sum(aseqUniq_n > 0);                                            	% initial number of distinct viral phenotypes
            s = 0;                                                                  % counter for while loop
            m = 1;                                                                  % counter to collect stats
            
            data_collect = zeros(S, 6);                                          	% matrix to collect time, # distinct genotypes, # distinct phenotypes, # not-infected cells, # infected cells, # viruses.
            data_collect(1,:) = [t, ntot, nAA, U, sum(I), sum(V)];                          
            meanDistance = zeros(1,S);                                              % mean number of mutations
            meanDistance(1) = sum(dtot .* V) / sum(V);
            meanFitness = zeros(1,S);                                               % mean fitness of population
            meanFitness(1) = sum(r .* V) / sum(V); 
            maxD = zeros(1,S);
            maxD(1) = max(dtot);
            maxR = zeros(1,S);
            maxR(1) = max(r);
            diversity = 1 - sum((V/sum(V)) .^2);
            
            t_transmission = t_L + t_I + rand;
            transmissionBoolean = false;
            
            transmission = struct;
            titer = struct;

            %% Start while loop ###########################################
            while t < T
  
                 s = s + 1;                                               
                 % Calculate total reaction rate
                 Rtot = (a*U + b)*sum(V) + c*sum(I) + sum(r.*I);  
                 % Time step is exponentially distributed and depends on 
                 % the total reaction rate
                 dt = (1/Rtot) * log(1/rand(myStream));
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
                     
                     % R1: infection by strain i
                     z = z + a*U*V(i);
                     if z > rnd                                                         
                         I(i) = I(i) + 1;                                               
                         V(i) = V(i) - 1;
                         U = U - 1;      
                         break
                     end
                     
                     % R2: replication of strain i
                     z = z + r(i)*I(i);
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
                            aseqUniq_n, aseqUniq_i, aseqUniq_r, distribution, fitnessFunction, wholeGenome);
                        break
                     end
                     
                     % R3: clearance of a cell infected by strain i
                     z = z + c*I(i);
                     if z > rnd
                         I(i) = I(i) - 1;
                         % if the strain has gone extinct, remove it:
                         if V(i) + I(i) == 0
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
                         break
                     end
                     
                     % R4: clearance of a free viral particle of strain i
                     z = z + b*V(i);
                     if z > rnd
                         V(i) = V(i) - 1;
                         % if the strain has gone extinct, remove it:
                         if V(i) + I(i) == 0 
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
                         break
                     end
                 end

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
                     data_collect(m,:) = [t, ntot, nAA, U, sum(I), sum(V)];
                     maxD(m) = max(dtot);
                     maxR(m) = max(r);

                     % Simpson's index as diversity
                     diversity(m) = 1 - sum((viralTiter/sum(viralTiter)).^2); 
                 end

            
            
                % transmission
                 if t >= t_transmission && transmissionBoolean == false
                    transmissionBoolean = true;

                    % get last occupied entry in V and I vectors 

                    tempY = Y(1:m, :);

                    % sample a random pool of 400 viruses.
                    totV = sum(V);
                    pool = randi(totV, [1, V0]);
                    % initiate arrays for transmission
                    transV = zeros(1, S); 
                    transR = zeros(1, S);
                    transD = zeros(length(pNames), S);
                    transDtot = zeros(1, S);

                    transSeq_loc = cell(1,S);
                    transAseq_loc = cell(1,S);
                    transSeq_mut = cell(1,S);
                    transAseq_mut = cell(1,S);
                    transSeq_nMut = zeros(1,S);
                    transAseq_nMut = zeros(1,S);

                    transAseqUniq_loc = cell(1,S);
                    % add 'dummy' 1 values:
                    transAseqUniq_mut = cell(1,S); transAseqUniq_mut(:) = {'1'};
                    transAseqUniq_nMut = zeros(1,S);
                    transAseqUniq_n = zeros(1,S); 
                    transAseqUniq_i = cell(1,S); 
                    transAseqUniq_r = zeros(1,S);

                    g = 0;
                    p = 0;
                    nu = 1;
                    for i = 1:length(V)
                        newg = g + V(i);
                        numTransmitted = sum(pool > g & pool <= newg);
                        if numTransmitted > 0
                            p = p + 1;
                            transV(p) = numTransmitted;
                            transR(p) = r(i);
                            transD(:,p) = d(:,i);
                            transDtot(p) = dtot(i);

                            transSeq_loc{p} = seq_loc{i};
                            transSeq_mut{p} = seq_mut{i};
                            transSeq_nMut(p) = seq_nMut(i);

                            transAseq_nMut(p) = aseq_nMut(i);
                            transAseq_loc{p} = aseq_loc{i};
                            transAseq_mut{p} = aseq_mut{i};

                            % Unique arrays:
                            alreadyPresent = false;
                            for j = 1 : nu
                                if isequal(transAseqUniq_loc{j}, aseq_loc{i})
                                    if isequal(transAseqUniq_mut{j}, aseq_mut{i})
                                        alreadyPresent = true;
                                        transAseqUniq_n(j) = transAseqUniq_n(j) + numTransmitted;
                                        transAseqUniq_i{j} = [transAseqUniq_i{j}, p];
                                        break % for loop over unique AA seqs
                                    end  
                                end
                            end

                            if ~alreadyPresent
                                transAseqUniq_loc{nu} = aseq_loc{i};
                                transAseqUniq_mut{nu} = aseq_mut{i};
                                transAseqUniq_nMut(nu) = aseq_nMut(i);
                                transAseqUniq_n(nu) = numTransmitted;
                                transAseqUniq_i{nu} = p;
                                transAseqUniq_r(nu) = r(i);
                                nu = nu + 1;
                            end                    
                        end
                        g = newg;
                    end
                    % remove the 'dummy' 1 values again
                    transAseqUniq_mut(strcmp(transAseqUniq_mut, '1')) = {[]};

                    % transmission Y array:
                    alive = (transV ~= 0);
                    shortV = transV(alive);
                    transY = zeros(S, La + 1);                                     	% number of free viral particles per strain
                    for j = 1:max(transDtot) + 1
                        dist = j - 1;
                        d_logical = (transDtot(alive) == dist);
                        sumV = sum(shortV(d_logical));
                        transY(1,j) = sumV;
                    end      
                    maxV = find(transV, 1,'last');

                    transmission.('V') = transV(1:maxV);
                    transmission.('r') = transR(1:maxV);
                    transmission.('seq_loc') = transSeq_loc(1:maxV);
                    transmission.('seq_mut') = transSeq_mut(1:maxV);
                    transmission.('tTransmission') = t_transmission;
                    transmission.('transD') = transD(:, 1:maxV);
                    transmission.('transR') = transR(1:maxV);

                    %transmission.('aseq_loc') = transAseq_loc(1:maxV);
                    %transmission.('aseq_mut') = transAseq_mut(1:maxV);
                end % if transmission
                
            end % while loop

            % Remove zeros at the end of vectors
            maxY = find(sum(Y) > 0, 1, 'last');
            Y = Y(1:m, :);          
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
            titer.t = data_collect(:,1);
            titer.ntot = data_collect(:,2);
            titer.nAA = data_collect(:,3);
            titer.U = data_collect(:,4);
            titer.I_sum = data_collect(:,5);
            titer.V_sum = data_collect(:,6);
            titer.Mu = mu;
            titer.t_end = t;
            titer.statY = statY(1:maxY);
            titer.relativeY = relativeY(:, 1:maxY);
            titer.statD = statD;
            titer.statR = statR;
            titer.maxD = maxD;
            titer.maxR = maxR;
            titer.diversity = diversity;
            titer.statDiv = statDiv;
            [titer.V_peak, peakIndex] = max(titer.V_sum);
            titer.V_peakTime = titer.t(peakIndex);
            
            fname = [outputFolder, '/Gillespie', num2str(k), '.mat'];
            save(fname, 'titer', 'params', 'transmission', '-v7.3');
        end % for loop

        disp(['Simulation time ', num2str(toc), ' seconds.'])
end