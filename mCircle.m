function [Roundness, Area, Centroid, Solidity, ERoundness] = mCircle(Filename)
%MCIRCLE Measure circularity of an object in an image
%   Read in the image
%   Filter it down
%   Apply the measurements
%   Return the results

I=imread(Filename);
%Convert to Binary
%BB=im2bw(I);
%blue channel
BB = I(:, :, 3);
%Invert the image - white circles on black
B=imcomplement(BB);                                                 
%Fill the holes
C=imfill(B,'holes');

%threshold
BW = im2bw(C, 0.9);

 %Label the image
[Label,Total]=bwlabel(BW,8);
%Object Number
num=1;

Sdata=regionprops(Label,'all');

Area = Sdata(num).Area;

Centroid = Sdata(num).Centroid;

Roundness=(4*Sdata(num).Area*pi)/Sdata(num).Perimeter.^2;

%estimate perimeter using an ellipse
thisBB=Sdata(num).BoundingBox;
a = thisBB(3)/2;
b = thisBB(4)/2;
p = pi *(3*(a+b)- sqrt((3*a+b)*(a+3*b)));
ERoundness = (4*Sdata(num).Area*pi)/p.^2;

Solidity = Sdata(num).Solidity;

%display(Roundness);