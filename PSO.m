function solution = PSO(fobj,nPop,nFEs,dim,LB,UB)
%% Parameters
NP = nPop;
MaxIter = round(nFEs/nPop);
c1 = 0.8;
c2 = 0.8;

Dim = dim;

wMax = 1.0;         %Max inirtia weight
wMin = 0.5;         %Min inirtia weight

%% Defined lower bound and upper bound.
LB = repmat(LB,NP,1);
UB = repmat(UB,NP,1);

%% Initialize swarm population randomly
population =  LB+(UB-LB).*rand(NP,Dim);
ibest = Inf;

%% Randomly initialise velocity
vmin = LB.*ones(NP,Dim);
vmax = UB.*ones(NP,Dim);
velocity = zeros(NP,Dim);

%% Evaluate initial population
fvalue = Inf(NP,1);
fbestval = Inf;
nfe = 0;
for i = 1:NP,
    fvalue(i) = fobj(population(i,:));
    nfe = nfe+1;
    % Finding best particle in initial population
    if fvalue(i) <= fbestval
        fbestval = fvalue(i);
        ibest = i;
    end
end

gbest = repmat(population(ibest,:),[NP 1]);
bestPara = population(ibest,:);
pbest = population; % Initializing Best positions matrix
fpbest = fvalue;    % Initializing the corresponding function values

%% Main loop
iteration = 0;
while nfe < nFEs
    
    R1 = rand(NP,Dim);
    R2 = rand(NP,Dim);
    
    iteration = iteration+1;
    w = (MaxIter-iteration+1)/MaxIter * (wMax - wMin) + wMin;
    
    velocity = w*(velocity + c1*R1.*(pbest-population) + c2*R2.*(gbest-population));
    velocity = ( (velocity <= vmin).*vmin ) + ( (velocity > vmin).*velocity );
    velocity = ( (velocity >= vmax).*vmax ) + ( (velocity < vmax).*velocity );
    
    %% Update the swarm particle
    population = population + velocity;
    
    population(population>UB)=UB(population>UB);
    population(population<LB)=LB(population<LB);
    
    %% Evaluate the new swarm
    for i = 1:NP,
        fvalue(i) = fobj(population(i,:));
        nfe = nfe+1;
    end
    
    %% Updating the pbest for each particle
    changeRows = fvalue < fpbest;
    pbest(changeRows,:) = population(changeRows,:);
    fpbest = fpbest.*~changeRows + fvalue.*changeRows;
    
    %% Updating best particle gbest
    [fbestval_, ibest] = min(fpbest);
    %gbest = repmat(population(ibest,:),[NP 1]);
    %bestPara = population(ibest,:);
    if fbestval_ < fbestval
        fbestval = fbestval_;
        gbest = repmat(population(ibest,:),[NP 1]);
        bestPara = population(ibest,:);
    end
    
    if mod(iteration,10)==0
        fprintf('iteration: %d, fBest: %f \n',iteration,fbestval);
    end
end
solution = bestPara;
end
