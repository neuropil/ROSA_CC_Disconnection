function [perDC , perABL] = UngETal_Analysis_06152021(caseSel , dataRoot)
%
% Code and data will live on SanDisk - NeuroWork
% All that is needed is the Directory Root
% dataRoot = 'G';
dataPATH = [dataRoot , '\' , 'ROSA_CC_Visualase\ROSA_Laser_CaseReport\AnalysisJAT\'];
path1 = [dataRoot , '\ROSA_CC_Visualase\Code_repository\'];
path2 = [dataRoot , '\ROSA_CC_Visualase\Code_repository\NIfTI_tools'];
% OUTPUTS - 



cd(dataPATH);
addpath(path1)
addpath(path2)
%% Corpus callosum ROI
if caseSel == 1
    preOP = 'CC2.nii.gz'; % Case 1
else
    preOP = 'BCC_Preop_SegBilateral_06112021.nii.gz'; % Case 2
end

ccnii = niftiread(preOP);
ccinfo = niftiinfo(preOP);
ccimg = ccnii;

ccimg1 = ccimg;
ccimg1(ccimg == 4) = 1;

ccimgbin = Process_MRI_maskBinary_CC(ccimg1);
cc_volume = sum(sum(sum(ccimgbin)));

ccpxINFO = ccinfo.PixelDimensions;

cc_NormVol = round(cc_volume * (ccpxINFO(1) * ccpxINFO(2) * ccpxINFO(3)));

%% Ablation ROI

if caseSel == 1
    postOP = 'ablation2.nii.gz'; % Case 1
else
    postOP = 'Postop_seg.nii.gz'; % Case 2
end

ablii = niftiread(postOP);
ablinfo = niftiinfo(postOP);
ablimg = ablii;

ablimg1 = ablimg;
ablimg1(ablimg == 1) = 1;

ablimgbin = Process_MRI_maskBinary_CC(ablimg1);
% abl_volume = sum(sum(sum(ablimgbin)));
abpxINFO = ablinfo.PixelDimensions;
% abNorm = round(abl_volume * (abpxINFO(1) * abpxINFO(2) * abpxINFO(3)));

%% Overlap ROI
% ablF = flipud(fliplr(ablimgbin));
% ccF = flipud(fliplr(ccimgbin));
if caseSel == 1
    ccF = rot90(ccimgbin,2); % 2 for case 1
    ablF = rot90(ablimgbin,2);
else
    ccF = rot90(ccimgbin,3); % 2 for case 1
    ablF = rot90(ablimgbin,3);
end
% overlapROI = ccimgbin & ablimgbin;
overlapROI = ccF & ablF;

%% Compute disconnection
ABcoordLocs = array2table(nan(size(overlapROI,1),6),'VariableNames',...
    {'MaxLx','MaxLy','MaxLz','MinLx','MinLy','MinLz'});

CCcoordLocs = array2table(nan(size(overlapROI,1),6),'VariableNames',...
    {'MaxLx','MaxLy','MaxLz','MinLx','MinLy','MinLz'});

if caseSel == 1
    sliceCOUNT = size(overlapROI,1);
else
    sliceCOUNT = size(overlapROI,3);
end

for oi = 1:sliceCOUNT
    
    if caseSel == 1
        tmpSlice = squeeze(overlapROI(oi,:,:)); % Both combined
    else
        tmpSlice = squeeze(overlapROI(:,:,oi));
    end
    
    close all
    if ~any(tmpSlice)
        continue
    else
        % Extract Ablation Segmentation
        if caseSel == 1
            tmpSr = rot90(tmpSlice,1);
        else
            tmpSr = tmpSlice;
        end
        
        % Extract CC Segmentation
        if caseSel == 1
        ccimgbinT = squeeze(ccF(oi,:,:)); % Corpus Callosum
        
        % Rotate CC Segmentation such that it aligns with Ablation
        ccimgR = rot90(ccimgbinT,1);
        else
            ccimgbinT = squeeze(ccF(:,:,oi)); 
            ccimgR = ccimgbinT;
        end
        
        % Show overlay
        imshowpair(tmpSr,ccimgR,'falsecolor')
        hold on
        text(10, 10, 'Corpus Callosum', 'Color', 'Magenta')
        text(10, 25, 'Lesion', 'Color', 'White')
        
        % Extract Blobs
        bOver = bwboundaries(tmpSr); % Ablation
        bCC = bwboundaries(ccimgR); % CC 
        
        % Combine blobs for Ablation
        [overall] = combineCELLS(bOver);
        % Combine blobs for CC
        [ccAll] = combineCELLS(bCC);

        % Maximum anterior point of Ablation - Find most anterior point
        % Use as index to extract the superior location
        [ABcoordLocs.MaxLx(oi),maxInd] = max(overall(:,2));
        % Superior index for maximum anterior point of Ablation
        ABcoordLocs.MaxLz(oi) = overall(maxInd,1);
        % Minimum posterior point of Ablation - Find most posterior point
        [ABcoordLocs.MinLx(oi),minInd] = min(overall(:,2));
        ABcoordLocs.MinLz(oi) = overall(minInd,1);
        
        % Plot points
        hold on
        % Anterior Point
        plot(ABcoordLocs.MaxLx(oi),ABcoordLocs.MaxLz(oi),...
            '.', 'MarkerSize',30,'MarkerEdgeColor','blue',...
            'MarkerFaceColor',[0.6 0.6 1]);
        % Posterior Point
        plot(ABcoordLocs.MinLx(oi),ABcoordLocs.MinLz(oi),...
            '.', 'MarkerSize',30,'MarkerEdgeColor','green',...
            'MarkerFaceColor',[0.6 1 0.6]);
        % Slice Number
        ABcoordLocs.MaxLy(oi) = oi;
        ABcoordLocs.MinLy(oi) = oi;
           
        % CC Anterior and Posterior extents
        [CCcoordLocs.MaxLx(oi),~] = max(ccAll(:,2));
        [CCcoordLocs.MinLx(oi),~] = min(ccAll(:,2));

    end

end

% MOST ANTERIOR POINT
[overlapMAXx, ~] = max(ABcoordLocs.MaxLx);
% overlapMAXy = ABcoordLocs.MaxLy(MaxxInd);
% overlapMAXz = ABcoordLocs.MaxLz(MaxxInd);

% MOST POSTERIOR POINT
[overlapMINx, ~] = min(ABcoordLocs.MinLx);
% overlapMINy = ABcoordLocs.MinLy(MinnInd);
% overlapMINz = ABcoordLocs.MinLz(MinnInd);

ccMAX = max(CCcoordLocs.MaxLx);
ccMIN = min(CCcoordLocs.MinLx);
%% for C1 84 to 145
%% for C2 
% Used in CC Disconnection Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overlapDist = overlapMAXx - overlapMINx;
ccDist = ccMAX - ccMIN;
perDC.perDiscect = overlapDist/ccDist;
perDC.Ablmm = overlapDist;
perDC.CCmm = ccDist;

%% 
coverageBin = ccimgbin & ablimgbin;
% Used in CC Ablation Volume Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cc_volume;
coverageVol = sum(sum(sum(coverageBin)));
covNorm = round(coverageVol * (abpxINFO(1) * abpxINFO(2) * abpxINFO(3)));
perABL.percentCCabl = covNorm/cc_NormVol;
perABL.AblVol = covNorm;
perABL.CCVol = cc_NormVol;


%%
newOver = zeros(size(overlapROI),'logical');
newCC = zeros(size(overlapROI),'logical');
for oi = 1:sliceCOUNT
    
    if caseSel == 1
        tmpSlice = squeeze(overlapROI(oi,:,:));
    else
        tmpSlice = squeeze(overlapROI(:,:,oi));
    end
    
    
    if ~any(tmpSlice)
        continue
    else
        
        if caseSel == 1
            ccimgbinT = squeeze(ccF(oi,:,:));
            newOver(oi,:,:) = tmpSlice;
            newCC(oi,:,:) = ccimgbinT;
        else
            ccimgbinT = squeeze(ccF(:,:,oi));
            newOver(:,:,oi) = tmpSlice;
            newCC(:,:,oi) = ccimgbinT;
        end
        
    end
end


[ CC ] = extract3dobject_CC(newCC,1,0);
[ blobPointsCC , blobBoundsCC , ~ ] = FreeSurf_Extract_CC( CC , 0.5 );

[ LAZ ] = extract3dobject_CC(newOver,1,0);
[ blobPointsL , blobBoundsL , ~ ] = FreeSurf_Extract_CC( LAZ , 0.5 );

trisurf(blobBoundsCC,blobPointsCC(:,2),blobPointsCC(:,1),blobPointsCC(:,3),...
    'Facecolor',[0 0 0],...
    'Edgecolor','none',...
    'FaceAlpha',0.10);

hold on

trisurf(blobBoundsL,blobPointsL(:,2),blobPointsL(:,1),blobPointsL(:,3),...
    'Facecolor',[1 0 0],...
    'Edgecolor','none',...
    'FaceAlpha',1);

% Column 1 of blob is the medial/lateral -> Need mid point
% Column 2 of blob is anterior/posterior -> Greater is anterior
% Column 3 of blob is the supeior/inferior -> Greater is superior
antPostPoints = blobPointsL(:,2);
[anterNum, ~] = max(antPostPoints);
anteriorINDs = blobPointsL(antPostPoints == anterNum,:);
AntzPoint = round(mean(anteriorINDs(:,3)));
AntmidPoint = median(anteriorINDs(:,1));

[posterNum, ~] = min(antPostPoints);
posteriorINDs = blobPointsL(antPostPoints == posterNum,:);
PostzPoint = round(mean(posteriorINDs(:,3)));
PostmidPoint = median(posteriorINDs(:,1));

scatter3(posterNum, PostmidPoint, PostzPoint , 200,'b','filled')
scatter3(anterNum, AntmidPoint, AntzPoint , 200,'b','filled')
% Plot Disconnection line
% plot3([anterNum,posterNum],[AntmidPoint,PostmidPoint],[AntzPoint,PostzPoint],'k')

zticklabels([])
xticklabels([])
yticklabels([])

end





function [combinedCell] = combineCELLS(inputCEll)

numCells = length(inputCEll);
sizeCells = cellfun(@(x) size(x,1), inputCEll, 'UniformOutput', true);
totalVec = sum(sizeCells);

startI = 1;
stopI = sizeCells(1);
combinedCell = zeros(totalVec, 2);
for ci = 1:numCells
    
    tmpCell = inputCEll{ci};
    combinedCell(startI:stopI,:) = tmpCell;
    startI = stopI + 1;
    if ci == numCells
        continue
    else
        stopI = startI + (sizeCells(ci + 1)) - 1;
    end

end


end