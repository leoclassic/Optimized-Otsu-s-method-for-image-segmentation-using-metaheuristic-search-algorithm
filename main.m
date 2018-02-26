close all

level = 3;

IMG = imread('peppers.png');

%RGB umage
IM = IMG;
figure,imshow(IM);
I = segmentationFcn( IM, @PSO, level );
figure,imshow(I)

%Grey image
IM = rgb2gray(IMG);
figure,imshow(IM);
I = segmentationFcn( IM, @PSO, level );
figure,imshow(I)










