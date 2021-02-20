% blankimage3 = imread('I2.png');
% I = rgb2gray(blankimage3);
% BW1 = edge(I,'sobel');
% imwrite(BW1,'imagem11.png');

I = imread('SobelScan.png');
J = imcomplement(I);
imwrite(J,'imagem11.png');