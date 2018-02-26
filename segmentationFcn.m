function [ output ] = segmentationFcn( img, optimization, level )
MAX = 256;

NP = 100;   % the population size
nFEs = 50000; % the maximum number of function evaluations


prob = probCalculator(img,MAX);

fobj = @(solution) fitnessFnc( solution, img, prob, level );

dim = (level-1)*size(img,3);
LB = ones(1,dim);
UB = MAX*ones(1,dim);

solution = optimization(fobj,NP,nFEs,dim,LB,UB);
solution = fix(solution);

if size(img,3) == 1 % grayscale image
    solution = sort(solution);
    output = imageGRAY(img,solution);
else % RGB image
    nPar = level-1;
    
    xR = solution(1:nPar);
    xG = solution(nPar+1:nPar*2);
    xB = solution(nPar*2+1:end);
    
    xR = sort(xR);
    xG = sort(xG);
    xB = sort(xB);
    output = imageRGB(img,xR,xG,xB);
end

end

function fitness = fitnessFnc( solution, img, prob, level )
MAX = 256;

solution = fix(solution);

if size(img,3) == 1 % grayscale image
    fitness = 0;
    for i = 1:level
        if i == 1
            I = 1 : solution(1);
        elseif i == level
            I = solution(level-1)+1 : MAX;
        else
            I = solution(i-1)+1 : solution(i);
        end
        tmp = sum(prob(I)) * (sum(I.*prob(I)/sum(prob(I))) - sum((1:MAX).*prob(1:MAX)))^2;
        fitness = fitness + tmp;
    end
else % RBG image
    nPar = level-1;
    
    xR = solution(1:nPar);
    xG = solution(nPar+1:nPar*2);
    xB = solution(nPar*2+1:end);
    
    I = cell(1,3);
    fitness = zeros(1,3);
    for i = 1:level
        if i == 1
            I{1} = 1 : xR(1);
            I{2} = 1 : xG(1);
            I{3} = 1 : xB(1);
        elseif i == level
            I{1} = xR(level-1)+1 : MAX;
            I{2} = xG(level-1)+1 : MAX;
            I{3} = xB(level-1)+1 : MAX;
        else
            I{1} = xR(i-1)+1 : xR(i);
            I{2} = xG(i-1)+1 : xG(i);
            I{3} = xB(i-1)+1 : xB(i);
        end
        for j = 1:3
            tmp = sum(prob(j,I{j})) * (sum(I{j}.*prob(j,I{j})/sum(prob(j,I{j}))) - sum((1:MAX).*prob(j,1:MAX)))^2;
            fitness(j) = fitness(j) + tmp;
        end
    end
    fitness = sum(fitness);
end
fitness = -fitness;
end

function prob = probCalculator(img,MAX)

if size(img,3) == 1 % grayscale image
    [N,~] = imhist(img(:,:,1));
elseif size(img,3) == 3 % RGB image
    [nR,~] = imhist(img(:,:,1));
    [nG,~] = imhist(img(:,:,2));
    [nB,~] = imhist(img(:,:,3));
end
nPixel = size(img,1)*size(img,2);

prob = zeros(size(img,3),MAX);
if size(img,3) == 1 % grayscale image
    prob = N'./nPixel;
else % RGB image
    prob(1,:) = nR'./nPixel;
    prob(2,:) = nG'./nPixel;
    prob(3,:) = nB'./nPixel;
end
end

function imgOut=imageRGB(img,Rvec,Gvec,Bvec)
imgOutR=img(:,:,1);
imgOutG=img(:,:,2);
imgOutB=img(:,:,3);

Rvec=[0 Rvec 256];
for iii=1:size(Rvec,2)-1
    at=find(imgOutR(:,:)>=Rvec(iii) & imgOutR(:,:)<Rvec(iii+1));
    imgOutR(at)=Rvec(iii);
end

Gvec=[0 Gvec 256];
for iii=1:size(Gvec,2)-1
    at=find(imgOutG(:,:)>=Gvec(iii) & imgOutG(:,:)<Gvec(iii+1));
    imgOutG(at)=Gvec(iii);
end

Bvec=[0 Bvec 256];
for iii=1:size(Bvec,2)-1
    at=find(imgOutB(:,:)>=Bvec(iii) & imgOutB(:,:)<Bvec(iii+1));
    imgOutB(at)=Bvec(iii);
end

imgOut=img;

imgOut(:,:,1)=imgOutR;
imgOut(:,:,2)=imgOutG;
imgOut(:,:,3)=imgOutB;
end

function imgOut=imageGRAY(img,Rvec)
limites=[0 Rvec 255];
tamanho=size(img);
imgOut(:,:)=img*0;
k=1;
for i= 1:tamanho(1,1)
    for j=1:tamanho(1,2)
        while(k<size(limites,2))
            if(img(i,j)>=limites(1,k) && img(i,j)<=limites(1,k+1))
                imgOut(i,j,1)=limites(1,k);
            end
            k=k+1;
        end
        k=1;
    end
end
end