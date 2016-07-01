%clear all; clc;

function out = Center_largestblob(im)

[center_r,center_c] = size(im);    %size of the image
Imgr = center_r/2;   %center row
Imgc = center_c/2;   %center column
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% First largest object %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = bwlabel(im, 8);     %labeling of the image with L being the matrix 
stats = regionprops(L, 'Area');     %find area of each labelled object in L
allArea = [stats.Area];     %single vectore with all object area
[maxno largestBlobNo] = max(allArea);    %maxno = max area, largestBlobNo = object number with large area
outim = zeros(size(im),'uint8');    %creat image with all zeros
outim(find(L==largestBlobNo)) = 1;  %find largest object and place it in outim with all 1 values
im(find(L==largestBlobNo)) = 0;     %remove largest object from the original image
stats = regionprops(outim, 'Centroid');     %find centroid of the largest object
centroid = [stats.Centroid];     %single vectore with row and column value of the largest object
c = centroid(:,1);       %center of the row for largest object
r = centroid(:,2);       %center of the column for largest object
ED(1) = sqrt((Imgr - r).^2 + (Imgc - c).^2);     %distance between the center of image and the center of the largest object and stored in a vector
 

% figure; 
% % subplot(1,2,1); 
% imshow(im);     %original image after removing the largest object
% % subplot(1,2,2); 
% imshow(255*outim);      %largest object
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Second largest object %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %use original image after removing the largest object
L = bwlabel(im, 8);    %label rest of the object
stats = regionprops(L, 'Area');     %find area
allArea = [stats.Area];     
[maxno1 largestBlobNo1] = max(allArea);      %find second largest area in the image
outim1 = zeros(size(im),'uint8');
outim1(find(L==largestBlobNo1)) = 1;
im(find(L==largestBlobNo1)) = 0;
stats = regionprops(outim1, 'Centroid');
centroid1 = [stats.Centroid];
c1 = centroid1(:,1);
r1 = centroid1(:,2);
ED(2) = sqrt((Imgr - r1).^2 + (Imgc - c1).^2);
 
% figure; 
% % subplot(1,2,1); 
% imshow(im);
% % subplot(1,2,2); 
% imshow(255*outim1);
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Third largest object %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = bwlabel(im, 8);
stats = regionprops(L, 'Area');
allArea = [stats.Area];
[maxno2 largestBlobNo2] = max(allArea);      %find third largest area in the image
outim2 = zeros(size(im),'uint8');
outim2(find(L==largestBlobNo2)) = 1;
im(find(L==largestBlobNo2)) = 0;
stats = regionprops(outim2, 'Centroid');
centroid2 = [stats.Centroid];
c2 = centroid2(:,1);
r2 = centroid2(:,2);
ED(3) = sqrt((Imgr - r2).^2 + (Imgc - c2).^2);
%  
% figure; 
% subplot(1,2,1); imshow(im);
% subplot(1,2,2); imshow(255*outim2);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Forth largest object %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L = bwlabel(im, 8);
% stats = regionprops(L, 'Area');
% allArea = [stats.Area];
% [maxno3 largestBlobNo3] = max(allArea)      %find forth largest area in the image
% outim3 = zeros(size(im),'uint8');
% outim3(find(L==largestBlobNo3)) = 1;
% im(find(L==largestBlobNo3)) = 0;
% stats = regionprops(outim3, 'Centroid');
% centroid3 = [stats.Centroid]
% c3 = centroid3(:,1)
% r3 = centroid3(:,2)
% ED(4) = sqrt((Imgr - r3).^2 + (Imgc - c3).^2)
%  
% figure; 
% subplot(1,2,1); imshow(im);
% subplot(1,2,2); imshow(255*outim3);
 
Dist = min(ED);%% min distance from center of the image to the center of the object
 
% As the shortest distance is between the center of the first object to the
% center of the image. Then first object image i.e. outim is the correct
% output and show that image as the final result. 
Image(1)= struct('ED', ED(1), 'Image', outim);
Image(2)= struct('ED', ED(2), 'Image', outim1);
Image(3)= struct('ED', ED(3), 'Image', outim2);
% Image(4)= struct('ED', ED(4), 'Image', outim3);
 
location=find(Dist==(ED));
 
%largest Blob
BiggestBlob= Image(location).Image;

out = 255*BiggestBlob;

 
% figure; 
% imshow(255*BiggestBlob);
%imshow(255*outim);


