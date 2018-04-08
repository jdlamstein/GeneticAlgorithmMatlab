%% Genetic Search Algorithm
% This genetic algorithm finds the global maximum of a function with many
% local extrema. The search uses random and ordered crossover to keep
% searching new possibilities while maintaining previous progress. The boolean,
% moatBool, creates or removes a moat. Depending on the random number generator, the algorithm can find 
% the global max in less than 100 iterations. On average, it takes about 
% 400 iterations to reach the global maximum. 
% The code tracks the center of mass of the chromosomes, which
% indicates the global maximum. When the standard deviation of the 
% center of mass drops below a threshold, the search stops. 

% Initialize Variables

N=100; % The matrix to be searched is NxN
xmax = 4; % arbitrary bounds for the function
ymax = 4;

%Boolean - add or remove a moat
moatBool =1; %Boolean: 1 creates the moat, 0 removes the moat
moatMax = 2; %Bounds for the moat
moatMin = 1;


factor = 20; % number of chromosomes
comLength = 30; % length of center of mass array
stepsize = 30;
cutoff = factor; % number of points sampled from data, 
% does not need to be the same as factor, but it works


% [x y] coordinates searching for the maximum
chromo = [ ceil(factor*rand(factor,1)), ceil(factor*rand(factor,1))]; 
quality = nan(1,factor); % measure the improvement of search
maxQuality = 0; % the highest expectation of quality that we have, this gets updated
com = ones(comLength,2)*(-xmax); % center of mass of the search chromosomes
comInd = zeros(comLength,2); % center of mass of the search chromosomes

k=0; % index
tolerance = 2000; % algorithm will last for 'tolerance' interations
mutateRate = .2; % probability that mutation will occur

h_old=plot(0,0); % initialize point on plot for chromosomes
g_old = plot(0,0); % initialize center of mass marker

hardStep = 30; % baseline stepsize which I found from guess and check
sortWeight = 4.5; % weight found by guess and check
comThresh = 1.25;

x = linspace(-xmax,xmax,N);
y = linspace(-ymax,ymax,N);
f = zeros(N); % function to be maximized.

for i = 1:N
    for j = 1:N
        %Setting function to be maximized. The global maximum is in the center. 
        f(i,j) = ( xmax-abs(x(i)) ) *cos(2*pi*x(i))+( ymax-abs(y(j)) ) *cos(2*pi*y(j));

% Moat to test if the algorithm can leap to the global max
        if moatBool ==1
            if  ( x(i) >-moatMax && x(i) < moatMax && y(j) > -moatMax && y(j) < -moatMin) ||...
                    ( x(i) >-moatMax && x(i) < moatMax && y(j) < moatMax && y(j) > moatMin ) ||...
                     ( x(i) > -moatMax && x(i) < -moatMin && y(j) > -moatMax && y(j) < moatMax) ||...
                     ( x(i) > moatMin && x(i) < moatMax && y(j) > -moatMax && y(j) < moatMin)

                 f(i,j) = 0;
            end
        end

    end
end
% draw function
figure
contour(x,y,f)
hold on;
colorbar

%% Begin Loop

while k< tolerance
%% Analysis of the chromosomes
    
    k = k+1;
    
        % Get quality of chromosomes, which means check the height of 
        % our points searching the function
    for i= 1:factor
         quality(i) = f(chromo(i,1),chromo(i,2)); % Quality is the height
         % our test points (chromosomes)
    end
    
    if min(quality)<0 % make all values of quality positive
        quality = quality - min(quality) +1;
    end
    fitness = quality ./ sum(quality);% fitness normalizes quality

    [sortedFitness, index] = sort(fitness,'descend'); % sort
    
    chromoSorted = chromo(index,:);
    
    if max(quality) > maxQuality
       % If the max we attain, is greater than the prior
       % we raise our standards
       maxQuality = max(quality);
       
       % This chromosome will be the foundation in searches 
       chromoMax = chromoSorted(1,:); 
    end
    
    % Center of mass - We expect the mean of the 
    % center of mass to land on the global maximum
    
    % indices of center of mass in function
    comInd(mod(k,comLength)+1,:) = ceil([mean(chromo(:,1)),mean(chromo(:,2))]); 
    %center of mass
    com(mod(k,comLength)+1,:) = [mean(x (chromo(:,1) ) ), mean( y(chromo(:,2))) ]; 
    % mean of COM
    comIndMean = [mean(comInd(:,1)),mean(comInd(:,2))];
    comMean = [mean(com(:,1)),mean(com(:,2))];
   
    % Measure distance between the test points or chromosomes
    distX = x(chromoSorted(1,1)) - x(chromoSorted(:,1)) ;
    distY = y(chromoSorted(1,2)) - y( chromoSorted(:,2)) ;

    % Rank chromosomes based on their distance from eachother
    % Chromosomes shouldn't be too close or too far and they should be
    % uniformly distributed
    rankDist = (distX.^2 + distY.^2);
    rankNorm = rankDist ./ sum(rankDist); % normalize

    [sortedRank, indexRank] = sort(rankNorm,'descend');
    
    % Weight - I treat the rank as 5x more important than the fitness
    % I prioritize the dispersion of the chromosomes
    % This is especially important in jumping the moat
    RF = (sortedFitness + sortWeight*sortedRank)/2 + 0.1; 
    
    %Sample data with weights RF without replacement. 
    weightedX = datasample(chromoSorted(:,1), cutoff,'Replace', false,'Weights', RF);
    weightedY = datasample(chromoSorted(:,2), cutoff,'Replace', false,'Weights', RF);
    
    % Randomly select a portion of the weighted data sample
    
    clear selectedX % These selections change size, clearing
    clear selectedY % saves us from outdated information
    
    selection = round((rand+1)/2 * cutoff);
    
    selectedX = weightedX(1:selection);
    selectedY = weightedY(1:selection);
    selectedY_2 = selectedY; % initialize gene
    
    % I use two methods of cross over. One switches genes randomly. The
    % other switches with order. The order is weighted, so this may be 
    % seen as the top genes swapping places. 
    
    % CROSS OVER WITH ORDER
    for i = 1:selection
        if mod(i,2) == 0
           holder = selectedY(i-1);
           selectedY(i-1) = selectedY(i);
           selectedY(i) = holder;
        end
    end
    % RANDOM CROSS OVER
    for i = 1:selection
        selectedY_2(i) = weightedY(ceil(rand * selection));       
    end
    
    if length(selectedY) + length(selectedY_2) < factor
        'length problem' % troubleshooting
    end
    
    child1 = horzcat(selectedX,selectedY); % Merge crossover with best order
    child2 = horzcat(selectedX,selectedY_2); % Merge random crossover
    % length of child should match factor
    child = vertcat(child1, child2(1:factor - selection,:)); 
    
    if length(child)~= factor
        'length child does not equal factor'
    end
    
%% Mutations
    
    mutate = child; % rename child to track errors. 
    
    % reproduction and mutation
    for i = 1:factor
        if rand < mutateRate
                hotCoord = ceil(rand*2); % random coordinate to be mutated
                mutate(i,hotCoord) = chromoMax(1,hotCoord) + (rand - 1/2) * stepsize;
        end
    end
    
%% Form new chromosome
    
    chromo = mutate;
    
    chromo = round(chromo); % ensure all chromo values are integers
    % Make sure all the chromosomes fit in the index bounds
    for i = 1:factor
        for j = 1:2
            if chromo(i,j)<1
                chromo(i,j)=hardStep+stepsize*(rand); 
            end
            if chromo(i,j)>N
                chromo(i,j) = N - hardStep+stepsize*(rand-1);
            end
        end
    end

    chromo = round(chromo);
    if (mod(k,50)==0) %show index every 50 iterations
        k 
    end
    %plot
    h = plot(x(chromo(:,1)),y(chromo(:,2)),'k.','MarkerSize',16);
    g = plot(comMean(1),comMean(2),'r.','MarkerSize',20);
    delete(h_old);
    delete(g_old);
    h_old = h;
    g_old = g;
    drawnow;
    
    % Use standard deviation of the indices of the center of mass means
    % for a threshold
    if std(comInd(:,1))<comThresh && std(comInd(:,2))<comThresh
        'Center of mass threshold reached.'
        break;
    end


end
['The search ended after ',num2str(k),' iterations.']
