%% Remarks
% In this example case, there are two  challenges: 1. extreme scaling (<0.2); 2. high outlier ratio (>90%). 
% Written by:   Liang Shen,   12,31,2020
% Citation: to be added

im1 = imread('00000005.jpg');
im2 = imread('00000006.jpg');
if size(im1,3)==3;  im1 = rgb2gray(im1);   end;
if size(im2,3)==3;  im2 = rgb2gray(im2);   end;
I1gray = single(im1);
I2gray = single(im2);

%% detection, extraction
keypoints1 = vl_sift(I1gray);   N1 = length(keypoints1);
keypoints2 = vl_sift(I2gray);   N2 = length(keypoints2);
[keypoints1, descriptors1] = vl_covdet(I1gray, 'Frames', keypoints1,...
                                      'descriptor', 'SIFT' ,...
                                      'estimateAffineShape', true, ...
                                      'estimateOrientation', true) ;
[keypoints2, descriptors2] = vl_covdet(I2gray, 'Frames', keypoints2,...
                                      'descriptor', 'SIFT' ,...
                                      'estimateAffineShape', true, ...
                                      'estimateOrientation', true) ; 
keypoints1 = keypoints1(:, 1:N1);   descriptors1 = descriptors1(:, 1:N1);
keypoints2 = keypoints2(:, 1:N2);   descriptors2 = descriptors2(:, 1:N2);

loc1=keypoints1([2,1,6,5,4,3],:)'; 
loc2=keypoints2([2,1,6,5,4,3],:)'; 

des1 = descriptors1';
des2 = descriptors2';

[loc1, des1] = delSamePoint(loc1, des1);                [loc2, des2] = delSamePoint(loc2, des2);

%% Pre-match
K_nei=3; Thr_sr=0.92;
[Simi_Nei_Ind_ori, Df_Table] = knnsearch(des2, des1, 'K', K_nei, 'NSMethod', 'exhaustive');
Corresp_original=Simi_Nei_Ind_ori(:,1);  

dist_ratios = Df_Table(:,1)./Df_Table(:,2);
[dist_in_sort,~]=sort(dist_ratios);   midInd=min(round(length(dist_in_sort)/2),10);
threshold=max(dist_in_sort(midInd),Thr_sr);
matchTable = find(dist_ratios < threshold);

XX1 = loc1(matchTable,:);   N=length(XX1);
XX2 = loc2(Corresp_original(matchTable),:);
Simi_Nei_Ind = Simi_Nei_Ind_ori(matchTable,:);
Sift_Ratio=dist_ratios(matchTable); %         2002Df_Table=Df_Table(matchTable,:);
Dist_desc=Df_Table(matchTable,1);
X_ind = matchTable;   
Y2_ind= Corresp_original(matchTable);   
% Visualization of the pre-matched results
Matching_Plot_No_global(im1,im2,loc1,loc2,[X_ind, Y2_ind]); 

%% GAFM Match
[C, sigma, iter, Ti] = GAFM(XX1, XX2, 150, 1e-4, 1, 0.9, Sift_Ratio, 0.75);
% Visualization of the GAFM-matched results
Matching_Plot_No_global(im1,im2,loc1,loc2,[X_ind(C(1,:),:), Y2_ind(C(2,:),:)]);
 