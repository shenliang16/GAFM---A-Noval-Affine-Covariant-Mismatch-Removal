%% Remarks
% In this example case, there are two  challenges: 1. extreme scaling (<0.2); 2. high outlier ratio (>90%). 
% Written by:   Liang Shen,   12,31,2020
% Citation: to be added


%% Note: we can't upload the compiled vlfeat since it is larger than 25M, therefore do the following two steps before running this demo 
% Step 1: please download at https://www.vlfeat.org/download/vlfeat-0.9.21-bin.tar.gz
% Step 2: Compile vlfeat 

%% add path of VL_feat   
addpath(genpath('.\vlfeat-0.9.21'));   

%% load image
im1 = imread('00000005.jpg');   im1_color = im1;
im2 = imread('00000006.jpg');   im2_color = im1;
if size(im1,3)==3;  im1 = rgb2gray(im1);   end;
if size(im2,3)==3;  im2 = rgb2gray(im2);   end;
I1gray = single(im1);
I2gray = single(im2);

keypoints1 = vl_sift(I1gray);   N1 = length(keypoints1);
keypoints2 = vl_sift(I2gray);   N2 = length(keypoints2);

%% Feature Frame Detection (local affine frame);         and Descriptor Extraction
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

% Check the correctness of the putative correspondence using the optimal transformation (given by the WxBS dataset);
threshold = 5;
B_opt = reshape([0.154203165467626	0.000188674532374101	652.316546762590	-0.00470356834532374	0.153508489208633	420.932949640288	-1.55845755395683e-06	-4.54489208633094e-06	1], 3,3)';
[Crr_GT,d1,d2] = VrfCrrsp(B_opt,  loc1(X_ind,:), loc2(Y2_ind,:), threshold);
% Visualization of the pre-matched results
Matching_Plot_No_global(im1,im2,loc1,loc2,[X_ind, Y2_ind],Crr_GT); 
title('Putative Matches')


%% GAFM Mismatch Removal
[C, sigma, iter, Ti] = GAFM(XX1, XX2, 150, 1e-4, 1, 0.9, Sift_Ratio, 0.75);
C_crrect = intersect(C, find(Crr_GT));          Num_c = length(C_crrect);
C_false = setdiff(C, find(Crr_GT));             Num_f = length(C_false);
C_sort = [C_crrect; C_false];
GAFM_flag = [ones(Num_c, 1); zeros(Num_f, 1)];
Matching_Plot_No_global(im1,im2,loc1,loc2,[X_ind(C_sort,:), Y2_ind(C_sort,:)], GAFM_flag);
title('Mismatch Removed')
Precision = Num_c / (Num_c + Num_f)
Recall = Num_c / sum(Crr_GT)