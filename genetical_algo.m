function [meas_GA,cost_GA,avgCostHist,minCostHist] = genetical_algo(n_measurements,T,pop_size,max_gen,n_part,n_draw)
% SpeedyGA is a vectorized implementation of a Simple Genetic Algorithm in Matlab
% Version 1.3
% Copyright (C) 2007, 2008, 2009  Keki Burjorjee
% Created and tested under Matlab 7 (R14). 

%  Licensed under the Apache License, Version 2.0 (the "License"); you may
%  not use this file except in compliance with the License. You may obtain 
%  a copy of the License at  

%  http://www.apache.org/licenses/LICENSE-2.0 

%  Unless required by applicable law or agreed to in writing, software 
%  distributed under the License is distributed on an "AS IS" BASIS, 
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
%  See the License for the specific language governing permissions and 
%  limitations under the License. 

%  Acknowledgement of the author (Keki Burjorjee) is requested, but not required, 
%  in any publication that presents results obtained by using this script 

%  Without Sigma Scaling, Stochastic Universal Sampling, and the generation of mask 
%  repositories, SpeedyGA faithfully implements the specification of a simple genetic 
%  algorithm given on pages 10,11 of M. Mitchell's book An Introduction to
%  Genetic Algorithms, MIT Press, 1996). Selection is fitness
%  proportionate.

if nargin < 6
    n_part = 250; %number of particles in the particle filter
    n_draw = 100; %number of draws in the MC
end

len=n_measurements;        % The length of the genomes  
popSize=pop_size;          % The size of the population (must be an even number)
maxGens=max_gen;                % The maximum number of generations allowed in a run
probCrossover=1;           % The probability of crossing over. 
probMutation=0.003;        % The mutation probability (per bit)
sigmaScalingFlag=1;        % Sigma Scaling is described on pg 168 of M. Mitchell's
                           % GA book. It often improves GA performance.
sigmaScalingCoeff=1;       % Higher values => less fitness pressure 

SUSFlag=1;                 % 1 => Use Stochastic Universal Sampling (pg 168 of 
                           %      M. Mitchell's GA book)
                           % 0 => Do not use Stochastic Universal Sampling
                           %      Stochastic Universal Sampling almost always
                           %      improves performance

visualizationFlag=1;       % 0 => don't visualize bit frequencies
                           % 1 => visualize bit frequencies

verboseFlag=0;             % 1 => display details of each generation
                           % 0 => run quietly
convergenceFlag=1;         % 1 => plot convergence curve
                           % 0 => does not

useMaskRepositoriesFlag=1; % 1 => draw uniform crossover and mutation masks from 
                           %      a pregenerated repository of randomly generated bits. 
                           %      Significantly improves the speed of the code with
                           %      no apparent changes in the behavior of
                           %      the SGA
                           % 0 => generate uniform crossover and mutation
                           %      masks on the fly. Slower.

% pre-generate two ?repositories? of random binary digits from which the  
% the masks used in mutation and uniform crossover will be picked. 
% maskReposFactor determines the size of these repositories.

maskReposFactor=5;
mutmaskRepos=rand(popSize,(len+1)*maskReposFactor)<probMutation;

% preallocate vectors for recording the average and maximum fitness in each
% generation
avgFitnessHist=zeros(1,maxGens+1);
maxFitnessHist=zeros(1,maxGens+1);

% the population is a popSize by len matrix of randomly generated boolean
% values

pop = zeros(popSize,len);
for i=1:popSize
    pop(i,:) = sort(randperm(T+1,len)-1);
end

% To identify copies in population
pop = sortrows(pop);
gen = 0 ; 
while gen<=maxGens  
    % evaluate the fitness of the population. The vector of fitness values 
    % returned  must be of dimensions 1 x popSize.
    fitnessVals=localFitnessFunction(pop);
     
    [maxFitnessHist(1,gen+1),maxIndex]=max(fitnessVals);
    avgFitnessHist(1,gen+1)=mean(fitnessVals,'omitnan');
     
    % display the generation number, the average Fitness of the population,
    % and the maximum fitness of any individual in the population
    if verboseFlag
        display(['gen=' num2str(gen,'%.3d') '   avgFitness=' ...
            num2str(avgFitnessHist(1,gen+1),'%3.3f') '   maxFitness=' ...
            num2str(maxFitnessHist(1,gen+1),'%3.3f') ]);
    end
    % Conditionally perform bit-frequency visualization
    if visualizationFlag
        figure(1)
        set (gcf, 'color', 'w');
        hold off
        histogram(pop,0:T,'Normalization','countdensity'); hold on;
        plot(pop(maxIndex,:)+0.5,0*pop(maxIndex,:)+popSize,'.','Markersize',25);
        axis([0 T 0 popSize]);
        title(['Generation = ' num2str(gen) ', Average Fitness = ' sprintf('%0.3f', avgFitnessHist(1,gen+1))]);
        ylabel('Frequency of measure in t');
        xlabel('time t');
        drawnow;
    end
    
    % Conditionally perform sigma scaling 
    if sigmaScalingFlag
        sigma=std(fitnessVals,'omitnan');
        if sigma~=0
            fitnessVals=1+(fitnessVals-mean(fitnessVals,'omitnan'))/...
            (sigmaScalingCoeff*sigma);
            fitnessVals(fitnessVals<=0)=0;
        else
            fitnessVals(~isnan(fitnessVals))=1;
        end
    end
    
    % Normalize the fitness values and then create an array with the 
    % cumulative normalized fitness values (the last value in this array
    % will be 1)
    fitnessVals(isnan(fitnessVals))=sum(fitnessVals,'omitnan');
    cumNormFitnessVals=cumsum(fitnessVals/sum(fitnessVals));
    
    % Use fitness proportional selection with Stochastic Universal or Roulette
    % Wheel Sampling to determine the indices of the parents 
    % of all crossover operations
    if SUSFlag
        markers=rand(1,1)+(1:popSize)/popSize;
        markers(markers>1)=markers(markers>1)-1;
    else
        markers=rand(1,popSize);
    end
    
    [~, parentIndices]=histc(markers,[0 cumNormFitnessVals]);
    parentIndices=parentIndices(randperm(popSize)); % shuffle  

    % deterimine the first parents of each mating pair
    firstParents=pop(parentIndices(1:popSize/2),:);
    % determine the second parents of each mating pair
    secondParents=pop(parentIndices(popSize/2+1:end),:);
    
    % CROSSOVER: COUNT PRESERVING CROSSOVER
    % determine which parents will contribute to crossover
    crossoverIndices = rand(popSize/2,1)<probCrossover;
    coupleList = 1:popSize/2;
    coupleList = coupleList(crossoverIndices);
    
    firstKids=firstParents;
    secondKids=secondParents;
    parfor i = coupleList
        [firstKid,secondKid] = slimCrossOver(firstParents(i,:),secondParents(i,:), firstKids(i,:),secondKids(i,:));
        firstKids(i,:)=firstKid;
        secondKids(i,:)=secondKid;
    end
    pop=sort([firstKids; secondKids],2);
    
    % implement mutations
    if useMaskRepositoriesFlag
        temp=floor(rand*len*(maskReposFactor-1));
        masks=mutmaskRepos(:,temp+1:temp+len);
    else
        masks=rand(popSize, len)<probMutation;
    end
    % masks(i,j)==1 iff pop(i,j) has to be mutated (0 elsewhere)
    pop = sort((1-masks).*pop + masks.*(unidrnd(T+1,popSize,len)-1),2);
    
    % Replace duplicates measurements and sort
    parfor i=1:popSize
        pop(i,:) = replace_duplicates(pop(i,:),T+1);
    end
    pop = sort(pop,2);
    pop = sortrows(pop); % to identify copies in population
    
    gen = gen+1; 
end

avgCostHist = -avgFitnessHist;
minCostHist = -maxFitnessHist;

meas_GA = sort(pop(maxIndex,:));
cost_GA = -maxFitnessHist(end);

%% plot and print
if convergenceFlag
    figure
    set(gcf,'Color','w');
    hold off
    plot(0:maxGens,-avgFitnessHist,'k-');
    hold on
    plot(0:maxGens,-maxFitnessHist,'c-');
    title('Minimum and Average Cost');
    xlabel('Generation');
    ylabel('Cost');
end

    function fitness = localFitnessFunction(pop)
       % function to MAXIMIZE
       % is designed to compute only once the fitness in cases of copies of
       % individuals
        [popLocSize,~] = size(pop);
        
        % first time that an individual appears
        indFirstCopy = find(sum( (pop(1:end-1,:)-pop(2:end,:)).^2 ,2)~=0)'+1;
        indFirstCopy = [1 indFirstCopy];
        
        % measurements of these first individuals
        firstMeas = pop(indFirstCopy,:);
        
        firstFitnesses = zeros(1,length(indFirstCopy)-1);
        parfor j = 1:length(indFirstCopy)
            meas = firstMeas(j,:);
            firstFitnesses(j)  = - MC_MSE_estimator(meas,T,n_draw,n_part);
        end
        
        % copy the fitnesses for similar individuals
        fitness = repelem(firstFitnesses,diff([indFirstCopy popLocSize+1]));
    end
end
