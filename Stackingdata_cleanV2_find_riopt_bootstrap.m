% totallength variable controls bootstrap sample number; 
% Line:"mod sample length" controls sample length (number of simulation)
% mod Lines:"Mod ri_opt given varying N_cluster range" (there are six) to change where in max N cluster for which N_cluster \in [1,N] can Ri_opt occur
clear
close all
clc
% importing
f = 1;

farray_all = {'201','208','214','215','216','220','221','222','223','224','225', '226','227','228',...
    '229','230','231','232','233','234','235', '236','237','238','239','240','241','242', '243', '244',...
    '245','246','247','248','249','250','251','254','255'};


ijk = 1;
totallength = 15; 
% msSSE_output = cell(length(farray_all) - ijk + 1, 1);
% msSSE_changepercentage_output = cell(length(farray_all) - ijk + 1, 1);
% numC_opt = zeros(length(farray_all) - ijk + 1, 1);
% numC_opt_fs = zeros(length(farray_all) - ijk + 1, 1);
% ri_opt = zeros(length(farray_all) - ijk + 1, 1);
% ri_opt_fs = zeros(length(farray_all) - ijk + 1, 1);
% ri_opt_fsms = zeros(length(farray_all) - ijk + 1, 1);

numC_opt = zeros(totallength, 1);
numC_opt_fs = zeros(totallength, 1);
ri_opt = zeros(totallength, 1);
ri_opt_fs = zeros(totallength, 1);
ri_opt_fsms = zeros(totallength, 1);

for iii = 1:totallength%ijk:length(farray_all)
    farray = farray_all(randsample(length(farray_all), 35)); %mod sample length
    for i = 1:length(farray)
        load(strcat('S',farray{i},'jump_cat'),'jump_cat');
        jump_Cat(i) = jump_cat;
        %jump over jump_coef(slope of lin fit to)
        load(strcat('S',farray{i},'ratio_jump'),'ratio_jump');
        fnR = fieldnames(ratio_jump);
        jump_Rat(i) = ratio_jump;

        load(strcat('S',farray{i},'ladder_jump'),'ladder_jump');
        fn = fieldnames(ladder_jump);
        for k=1:numel(fn)
            if( isnumeric(ladder_jump.(fn{k})) )
                ladder_jump.(fn{k}) = ladder_jump.(fn{k})';
            end
        end
        jump_Ladder(i) = ladder_jump;


        load(strcat('S',farray{i},'sbs'),'sbs');
        Sbs{i} = sbs;

        load(strcat('S',farray{i},'bcoord'),'bcoord');
        Bcoord{i} = bcoord;

        load(strcat('S',farray{i},'seg_sjcrk'),'seg_sjcrk');
        seg_Sjcrk{i} = seg_sjcrk;

        load(strcat('S',farray{i},'sj_ss_ini'),'sj_stickslip_ini');
        ss_Ini{i} = sj_stickslip_ini;

        load(strcat('S',farray{i},'sj_kernel2000'),'ySix');
        sj_Kernel2000{i} = ySix;

        load(strcat('S',farray{i},'kernelpks'),'kernelpeaks');
        kn_Pks{i} = kernelpeaks;

        load(strcat('S',farray{i},'kernelpks_leverarm'),'kernelpks_leverarm');
        kn_leverarm{i} = kernelpks_leverarm;

        load(strcat('S',farray{i},'sj_spatdist'),'sj_spatdist');
        sj_spacdist{i} = sj_spatdist;

        load(strcat('S',farray{i},'pf_ss_sj'),'pf_ss_sj');
        sj_ss_polyfit{i} = pf_ss_sj;

        load(strcat('S',farray{i},'ss_stressdrop'),'ss_stressdrop');
        ss_sd{i} = ss_stressdrop; %stickslip stressdrop

        load(strcat('S',farray{i},'seg_ms_sj'),'seg_ms_sj');
        seg_AE{i} = seg_ms_sj; %stickslip stressdrop
    end

    % cross covariance matrix
    %remove the NaN rows
    shock = cell(length(farray),1);
    FS = cell(length(farray),1);
    for i = 1:length(farray)
        shock{i} = [kn_leverarm{i},kn_Pks{i}',sj_spacdist{i}];
        shock{i}(isnan(shock{i}(:,3)), : ) = [];
        FS{i} = shock{i}((1:find(shock{i}(:,1)==0)-1 ),: );
    end

    for i = 1:length(farray)
        FSbar(i,:) = mean(FS{i},1);
        FSvar(i,:) = var(FS{i},1);
        FScov{i} = ((FS{i} - FSbar(i,:))' * (FS{i} - FSbar(i,:)))/ (length(FS{i}(:,1))-1); %same as cov(FS{i})    
        [FScorr{i},s(i,:)] = corrcov(FScov{i}); %Use matlab built in function to calculate correlation matrix
    end
    %
    writematrix(FScorr{end},'FScorr.csv');

    % calculate average and variance of FScorr with respect to number of
    % simulations (to see if it converges)
    for ii = 1:length(farray)
        for i = 1:size(FS{1},2) %1:5, 5 dimensions
            for j = 1:size(FS{1},2) %1:5, 5 dimensions
                for k = 1:ii%length(farray)
                    temp(k) = FScorr{k}(i,j);
                end
                FScorravg(i,j, ii) = mean(temp);
                FScorrvar(i,j, ii) = var(temp);
            end
        end
    end
    %{
    figure1 = figure(f); f=f+1;
    for i = 1:size(FS{1},2)
        for j = 1:size(FS{1},2)
            temp = zeros(length(farray),1);
            subplot(size(FS{1},2), size(FS{1},2), (size(FS{1},2)*i-(size(FS{1},2)-j)) );
            temp(:,1) = FScorravg(i, j, :);
            plot(temp);
            clear temp
        end
    end
    figurename = 'FSmean.png'; saveas(figure1,figurename);
    figure1 = figure(f); f=f+1;
    for i = 1:size(FS{1},2)
        for j = 1:size(FS{1},2)
            temp = zeros(length(farray),1);
            subplot(size(FS{1},2), size(FS{1},2), (size(FS{1},2)*i-(size(FS{1},2)-j)) );
            temp(:,1) = FScorrvar(i, j, :);
            plot(temp);
            clear temp
        end
    end
    figurename = 'FScovariance.png'; saveas(figure1,figurename);
    %}
    %{
    % PCA based on each foreshock
    FSnorm = cell(length(farray),1);
    for i = 1:length(farray)
        FSnorm{i} = normalize(FS{i}, 'range');
    end
    fsscore = cell(length(farray),1);
    fstsquared = cell(length(farray),1);
    for i = 1:length(farray)
        [fscoeff(:,:,i),fsscore{i}(:,:),fslatent(:,i),fstsquared{i}(:,:),fsexplained(:,i)] = pca(FSnorm{i});
    end
    %}
    %{
    % Calculate the total mean and variance of all foreshocks
    Ltotal = size(FS{1},1);
    for i = 2:length(farray)
        Ltotal = Ltotal + size(FS{i},1);
    end
    FStotal = zeros( Ltotal, 6);
    pltl = size(FStotal,1);
    j = 1;
    for i = 1:length(farray)
        FStotal(j:j+length(FS{i}(:,1))-1, 1:5) = FS{i};
        FStotal(j:j+length(FS{i}(:,1))-1, 6) = i;
        FStotal(j:j+length(FS{i}(:,1))-1, 7) = 1:size(FS{i},1);
        j = j+length(FS{i}(:,1));
    end
    FStotalmean = mean(FStotal,1);
    FStotalvar = var(FStotal,1);

    FStotalcov = cov(normalize(FStotal(:,1:5), 'range')');
    FStotalcor = corrcov(FStotalcov);
    %}
    % Find the prominant peaks in foreshock kernel density function. First
    % reorder the peaks, and then pick the top 5
    FSreorder = cell(length(farray),1);
    for i = 1:length(farray)
        FSreorder{i} = sortrows(FS{i}, 2, 'descend');
    end
    % Compile the prominant peaks in FS kernel density function
    tempn = 5;
    FShigh = zeros(tempn, 5, length(farray));
    for i = 1:length(farray)
        FShigh(:,:,i) = FSreorder{i}(1:tempn,:);
    end
    FShightotal = zeros(tempn*length(farray),5);
    j = 1;
    for i = 1:length(farray)
        FShightotal(j:(j+tempn-1), 1:5) = FShigh(:,:,i);
        FShightotal(j:(j+tempn-1), 6) = i;
        j = j+tempn;
    end
    %
    %{
    figure1 = figure(f);f=f+1;
    for i = 1:4
        subplot(2,2,i)
        histogram(FShightotal(:,i));
    end
    figurename = 'FShistogram.png'; saveas(figure1,figurename);
    %}
    %
    % Boxcox transformation on prominant foreshock peaks to guarantee normality
    FShightotalboxcox = zeros(tempn*length(farray), 5);
    FShightotalboxcox(:,1) = (FShightotal(:,1).^(0.5) - 1)/(0.5);
    FShightotalboxcox(:,2) = log(FShightotal(:,2));
    FShightotalboxcox(:,3:4) = FShightotal(:,3:4);
    FShightotalboxcox(:,5) = FShightotal(:,6);
    %{
    figure1 = figure(f);f=f+1; 
    for i = 1:4
        subplot(2,2,i)
        histogram(FShightotalboxcox(:,i));
    end
    figurename = 'FShighhistogram.png'; saveas(figure1,figurename);
    figure1 = figure(f);f=f+1;
    for i = 1:4
        for j = 1:4
            subplot(4, 4, 4*i-(4-j));
            scatter(FShightotalboxcox(:,i),FShightotalboxcox(:,j));
        end
    end
    figurename = 'FSmultivariatenormality.png'; saveas(figure1,figurename);
    %}
    % cross correlation based on FShightotalbc
    FShtbccov = zeros(4,4);
    FShtbccor = zeros(4,4);
    FShtbccov(:,:) = cov( FShightotalboxcox(:, 1:4) );
    FShtbccor(:,:) = corrcov( FShtbccov(:,:) );
    %
    FShbccov = zeros(4,4, length(farray));
    FShbccor = zeros(4,4, length(farray));
    for i = 1:length(farray)
        FShbccov(:,:,i) = cov( FShightotalboxcox((tempn*i-(tempn-1)):(tempn*i), 1:4) );
        FShbccor(:,:,i) = corrcov( FShbccov(:,:,i) );
    end

    FShightotalbc = zeros(tempn,4,length(farray)); %after boxcox transformation, split them back into various simulations
    for i = 1:length(farray)
        FShightotalbc(:,:,i) = FShightotalboxcox((tempn*i-(tempn-1)):(tempn*i), 1:4) ;
    end
    FShighmahal = zeros(length(farray),length(farray),tempn);
    FShighmahaltr = zeros(length(farray),length(farray));
    iFShtbccov = inv(FShtbccov);
    for j = 1:length(farray)
        for i = 1:length(farray)
            temp1 = FShightotalbc(:,:,j);
            temp2 = FShightotalbc(:,:,i); 
            temp3 = [temp1 - temp2];
            FShighmahal(j,i,:) = diag( temp3 *iFShtbccov*temp3' ); %The 3rd dimension is from 1:tempn
            FShighmahaltr(j,i) = sum(FShighmahal(j,i,:));
        end
    end
    %
    fsZ = linkage(FShighmahaltr); %FShighmahal is calculated after boxcox transformation
    %figure1 = figure(f);f=f+1; dendrogram(fsZ);figurename = 'FSlinkage.png'; saveas(figure1,figurename);
    %
    I = inconsistent(fsZ);
    %
    fsbcC = cluster(fsZ, 'cutoff', 0.8);
    %
    %{
    % Foreshocks closest to mainshocks (not necessarily the ones with highest magnitudes, just those closest)
    tempn = 5;
    FSclose = zeros(tempn, 5, length(farray));
    for i = 1:length(farray)
        FSclose(:,:,i) = FS{i}(end-tempn+1:end,:);
    end
    FSclosetotal = zeros(tempn*length(farray),5);
    j = 1;
    for i = 1:length(farray)
        FSclosetotal(j:(j+tempn-1), 1:5) = FSclose(:,:,i);
        FSclosetotal(j:(j+tempn-1), 6) = i;
        j = j+tempn;
    end
    %}
    
    %{
    figure(f);f=f+1;
    for i = 1:5
        subplot(2,3,i)
        histogram(FSclosetotal(:,i));
    end
    %}
    %{
    % Boxcox transformation on foreshocks closest to mainshock
    FSctbc = zeros(tempn*length(farray), 5);
    FSctbc(:,1) = (FSclosetotal(:,1).^(0.5) - 1)/(0.5);
    FSctbc(:,2) = log(FSclosetotal(:,2));
    FSctbc(:,3:4) = FSclosetotal(:,3:4);
    FSctbc(:,5) = FSclosetotal(:,6);
    %}
    %{
    figure1 = figure(f);f=f+1; 
    for i = 1:4
        subplot(2,2,i)
        histogram(FSctbc(:,i));
    end
    figurename = 'FSclosehistogram.png'; saveas(figure1,figurename);
    %}
%{
    % cross correlation based on FSctbc
    FSctbccov = zeros(4,4);
    FSctbccor = zeros(4,4);
    FSctbccov(:,:) = cov( FSctbc(:, 1:4) );
    FSctbccor(:,:) = corrcov( FSctbccov(:,:) );
    % mahalanobis distances and clustering based on FSclosetotalboxcox
    FScbc = zeros(tempn,4,length(farray));
    for i = 1:length(farray)
        FScbc(:,:,i) = FSctbc((tempn*i-(tempn-1)):(tempn*i), 1:4) ;
    end
    FSclosemahal = zeros(length(farray),length(farray),tempn);
    FSclosemahaltr = zeros(length(farray),length(farray));
    for j = 1:length(farray)
        for i = 1:length(farray)
            temp1 = FScbc(:,:,j);
            temp2 = FScbc(:,:,i);
            temp3 = [temp1 - temp2];
            FSclosemahal(j,i,:) = diag( temp3 *inv(FSctbccov)*temp3' ); %The 3rd dimension is from 1:tempn
            FSclosemahaltr(j,i) = sum(FSclosemahal(j,i,:));
        end
    end
    %
    fscZ = linkage(FSclosemahaltr);
    %figure1 = figure(f);f=f+1; dendrogram(fscZ);%figurename = 'FScloselinkage.png'; saveas(figure1,figurename);
%}
    % (keyword M) Mainshock characteristics
    cc = cell(length(farray),1);
    xx = cell(length(farray),1);
    yy = cell(length(farray),1);
    zz = cell(length(farray),1);
    for i = 1:length(farray)
        temp = zeros(size(seg_Sjcrk{i},1), 1);
        for j = 1:size(seg_Sjcrk{i},1)
            temp(j) = size(seg_Sjcrk{i}{j},1);
        end
        window_stickslip(i) = find(temp == max(temp));
        clear temp
    end

    % AE analysis
    AE_ss_each = zeros(length(farray),4);
    AE_mag = seg_AE{1}{window_stickslip(1)}(:,8);
    for i = 1:length(farray)
        AE = seg_AE{i}{window_stickslip(i)};
        AE_ss_each(i,1) = size(AE,1);
        AE_ss_each(i,2) = mean(AE(:,8));
        AE_ss_each(i,3) = var(AE(:,8));
        if i>1
            AE_mag_new = seg_AE{i}{window_stickslip(i)}(:,8); 
            AE_mag = cat(1, AE_mag, AE_mag_new);
        end
    end
    %{
    figure(f);f=f+1;histogram(AE_mag);
    figure(f);f=f+1;qqplot(AE_mag);
    %}
    % AE analysis contd, identify how far each AE is away from the distribution
    AE_mahal_fromeach_toall = cell(length(farray),1);
    for i = 1:length(farray)
        AE_mahal_fromeach_toall{i} = mahal(seg_AE{i}{window_stickslip(i)}(:,8) , rmoutliers(AE_mag));
        AE_ss_each(i, 4) = mean(AE_mahal_fromeach_toall{i}(:,:));
        if AE_ss_each(i, 2) < mean(AE_mag)
            AE_ss_each(i,4) = -abs(AE_ss_each(i,4));
        end
    end

    % The algorithm above manages to describe each MS AE as a 1x4 vector, which
    % allows it to compute with all kinds of linear algebra algorithms (no long varying in dimension)
    % AS_ss_each(:,4) is how much each MS AE is away from the distribution
    % (based on population). It is all positive for now, but adding a sign
    % based on the average of each MS AE mag to the average of the population
    % AE mag will allow us to describe all AE characterstics with one single
    % number, with the exception of number of AE, but that should already be
    % accounted for in the MS initiation sites (at least close in number)

    %
    for i = 1:length(farray)
        cc{i} = seg_Sjcrk{i}{window_stickslip(i)}(:,1);
        xx{i} = seg_Sjcrk{i}{window_stickslip(i)}(:,2);
        yy{i} = seg_Sjcrk{i}{window_stickslip(i)}(:,3);
        zz{i} = seg_Sjcrk{i}{window_stickslip(i)}(:,4);
    end

    porder = 4;
    % Decompose in Chebyshev polynomial space
    porder = 4;
    shapefcn = polyBasis('chebyshev',porder);
    ccn = cell(length(farray),1);
    xxn = cell(length(farray),1);
    zzn = cell(length(farray),1);
    yyn = cell(length(farray),1);
    xxncheby = cell(length(farray),1);
    zzncheby = cell(length(farray),1);
    yyncheby = cell(length(farray),1);
    for i = 1:length(farray)
        xxn{i} = normalize(xx{i}, 'range')*2 - 1;
        yyn{i} = normalize(yy{i}, 'range')*2 - 1;
        zzn{i} = normalize(zz{i}, 'range')*2 - 1;

        temp = linspace(-1 , 1, length(xx{i}))'; 
        ccn{i} = temp;
        G = zeros(length(xx{i}), porder);
        for j = 1:length(xx{i})
            G(j,:) = shapefcn(temp(j));
        end

        xxncoef(:,i) = inv(G'*G)*G'*xxn{i};
        yyncoef(:,i) = inv(G'*G)*G'*yyn{i};
        zzncoef(:,i) = inv(G'*G)*G'*zzn{i};

        xxncheby{i} = G*xxncoef(:,i);
        yyncheby{i} = G*yyncoef(:,i);
        zzncheby{i} = G*zzncoef(:,i);
    end
    %{
    figure1 = figure(f);f=f+1;ah1 = ceil(sqrt(length(farray)));
    set(gcf,'Position',get(0,'Screensize'));
    for i = 1:length(farray)
        subplot(ah1, ah1, i);
        plot( ccn{i} , xxn{i} ); hold on; 
        plot( ccn{i} , xxncheby{i} ); hold off
        %plot(1:length(cc{i}), zzp{i});
    end
    han=axes(figure1,'visible','off');han.Title.Visible='off';han.XLabel.Visible='on';han.YLabel.Visible='on';
    ylabel(han,'x coordinates SJ cracks, normalized to within [-1,1]');
    xlabel(han,'PFC cycles, normalized to within [-1,1]');
    figurename = 'MSchebyshevx.png'; saveas(figure1,figurename);

    figure1 = figure(f);f=f+1;ah1 = ceil(sqrt(length(farray)));
    set(gcf,'Position',get(0,'Screensize'));
    for i = 1:length(farray)
        subplot(ah1, ah1, i);
        plot( ccn{i} , yyn{i} ); hold on; 
        plot( ccn{i} , yyncheby{i} ); hold off
    end
    han=axes(figure1,'visible','off');han.Title.Visible='off';han.XLabel.Visible='on';han.YLabel.Visible='on';
    ylabel(han,'y coordinates SJ cracks, normalized to within [-1,1]');
    xlabel(han,'PFC cycles, normalized to within [-1,1]');
    figurename = 'MSchebyshevy.png'; saveas(figure1,figurename);
    %}
    % Use chebyshev polynomial coefficients to construct mainshock characteristic vector
    msL = porder * 2 + 1 + 1 + 1 + 1; %4 coefficients for x and y propagation, 1 for strain, 1 for stress drop, 1 for number of initiation sites, 1 for AE
    Mshock = zeros(length(farray),msL);
    Mshock(:, 1:porder*2) = [xxncoef' , yyncoef'];
    for i = 1:length(farray)
        Mshock(i,porder*2+1) = ss_sd{i}(end,1) - ss_sd{i}(1,1);
        Mshock(i,porder*2+2) = ss_sd{i}(1,2) - ss_sd{i}(end,2);
        Mshock(i,porder*2+3) = size(ss_Ini{i}, 1);
        Mshock(i,porder*2+4) = AE_ss_each(i,4);
    end
    %{
    figure1 = figure(f);f=f+1;
    ah1 = ceil(sqrt(msL));
    for i = 1:msL
        subplot(ah1,ah1,i); histogram(Mshock(:,i));
    end
    han=axes(figure1,'visible','off');han.Title.Visible='on';han.XLabel.Visible='off';han.YLabel.Visible='off';
    title('Histogram of each mainshock variable');
    figurename = 'MShistogram.png'; saveas(figure1,figurename);
    %
    figure1 = figure(f);f=f+1;
    for i = 1:msL
        for j = 1:msL
            subplot(msL,msL, msL*i-(msL-j));
            scatter(Mshock(:,i), Mshock(:,j),5,'filled');
        end
    end
    han=axes(figure1,'visible','off');han.Title.Visible='on';han.XLabel.Visible='off';han.YLabel.Visible='off';
    title('Scatter plot of each pair of mainshock variables');
    figurename = 'MSmultivariatenormal.png'; saveas(figure1,figurename);
    %}
    MScov = zeros(msL, msL, length(farray)-1);
    MScor = zeros(msL, msL, length(farray)-1);

    for i = 1:length(farray)-1
        MScov(:,:,i) = cov(Mshock(1:i+1,:));
        MScor(:,:,i) = corrcov(MScov(:,:,i)); %Use matlab built in function to calculate correlation matrix
    end
    MScoravg = zeros(msL,msL,length(farray)-1);
    MScorvar = zeros(msL,msL,length(farray)-1);
    for i = 1:length(farray)-1
        MScoravg(:,:,i) = mean(MScor(:,:,1:i),3);
        MScorvar(:,:,i) = var(MScor(:,:,1:i),0, 3);
    end
    %{
    figure1 = figure(f);f=f+1;
    set(gcf,'Position',get(0,'Screensize'));
    for i = 1:msL
        for j = 1:msL
            temp = zeros(length(farray)-1,1);
            subplot(msL,msL, msL*i-(msL-j));
            temp(:,1) = MScoravg(i,j,:);
            plot(temp);
        end
    end
    han=axes(figure1,'visible','off');han.Title.Visible='off';han.XLabel.Visible='on';han.YLabel.Visible='on';
    ylabel(han,'mainshock parameters: averages of correlation');xlabel(han,'increasing number of simulations');
    figurename = 'MScoravg.png'; saveas(figure1,figurename);

    figure1 = figure(f);f=f+1;
    set(gcf,'Position',get(0,'Screensize'));
    for i = 1:msL
        for j = 1:msL
            temp = zeros(length(farray)-1,1);
            subplot(msL,msL, msL*i-(msL-j));
            temp(:,1) = MScorvar(i,j,:);
            plot(temp);
        end
    end
    han=axes(figure1,'visible','off');han.Title.Visible='off';han.XLabel.Visible='on';han.YLabel.Visible='on';
    ylabel(han,'mainshock parameters: variances of correlation');xlabel(han,'increasing number of simulations');
    figurename = 'MScorvar.png'; saveas(figure1,figurename);
    %}
    MScorstd = MScorvar(:,:,end).^0.5;
    %
    iMScov = inv(MScov(:,:,end));
    for j = 1:length(farray)
        for i = 1:length(farray)
            temp = Mshock(j,:)-Mshock(i,:);
            %msdist(j,i) = temp*inv(MScov(:,:,end))*temp';
            msdist(j,i) = temp*iMScov*temp';
        end
    end
    msZ = linkage( msdist); 
    %figure1=figure(f);f=f+1;dendrogram(msZ);
    %figurename = 'MSlinkage.png'; saveas(figure1,figurename);
    % Cluster based on msZ
    msC = cluster(msZ, 'Maxclust', 4); %Maxclust, maximum number of clusters

    % Mainshock correlation between KDF peak and stress drop. very strong 0.92
    MSkdfpk = zeros(length(farray),1);
    for i = 1:length(farray)
        MSkdfpk(i) = max(shock{i}(:,2));
    end
    MSkdfsdcov = corrcov(cov([MSkdfpk, Mshock(:,10)]));
    %MSkdfsdcov = cov([MSkdfpk, Mshock(:,10)]);
    %MSkdfsdcor = corrcov(MSkdfsdcov);

    % Mainshock initiation sites vs stickslip strain change. moderate 0.4
    %figure(f);f=f+1;
    MSinilength = zeros(length(farray),1);
    for i = 1:length(farray)
        MSinilength(i) = ss_Ini{i}(end,1)-ss_Ini{i}(1,1);
    %     plot(ss_Ini{i}(:,1)-ss_Ini{i}(1,1),1:1:length(ss_Ini{i}(:,1)));
    %     hold on
    end
    MSsimul = corrcov(cov([MSinilength,Mshock(:,11)]));
    % MSsimulcov = cov([MSinilength,Mshock(:,11)]);
    % MSsimulcor = corrcov(MSsimulcov);
    %hold off

    % Box cox transform on Mshock
    MSbc = Mshock;
    MSbc(:,1) = ((Mshock(:,1)-min(Mshock(:,1))).^(0.625) - 1)/(0.625);  
    MSbc(:,2) = ((Mshock(:,2)-min(Mshock(:,2))).^(0.875) - 1)/(0.875); 
    MSbc(:,4) = ((Mshock(:,4)-min(Mshock(:,4))).^(0.5) - 1)/(0.5);
    MSbc(:,5) = ((Mshock(:,5)-min(Mshock(:,5))).^(1.5) - 1)/(1.5);
    MSbc(:,6) = ((Mshock(:,6)-min(Mshock(:,6))).^(0.4) - 1)/(0.4);
    MSbc(:,7) = ((Mshock(:,7)-min(Mshock(:,7))*1.5).^(-0.02) - 1)/(-0.02);%log((Mshock(:,7)-min(Mshock(:,7))));
    MSbc(:,12) = ((Mshock(:,12)-min(Mshock(:,12))).^(2.75) - 1)/(2.75);
    %{
    figure1=figure(f);f=f+1;ah1 = ceil(sqrt(msL));
    for i = 1:msL
        subplot(ah1,ah1,i); histogram(MSbc(:,i));
    end
    title('Histogram of each mainshock variable (after Box-Cox transformation)');
    figurename = 'MSbchistogram.png'; saveas(figure1,figurename);
    %}
    % %% Code for figure out the lambda in box cox transformation
    % temp = linspace(-3, 3, 49); %a range of lamdba for testing
    % aa1 = zeros(length(Mshock(:,12)), 49);
    % tempn = 12; %1-12, each corresponds to a column in MShock
    % for i = 1:49
    %     aa1(:,i) = ((Mshock(:,tempn)-min(Mshock(:,tempn))).^(temp(i)) - 1)/(temp(i));
    % end
    % figure(f);f=f+1; 
    % for i = 1:49
    %     subplot(7,7,i);histogram(aa1(:,i)); %Find the lambda (temp(i)) that
            % % results in most normal destribution and plug in back in to the
            % %code above
    % end
    %{
    figure1 = figure(f);f=f+1;
    for i = 1:msL
        for j = 1:msL
            subplot(msL,msL, msL*i-(msL-j));
            scatter(MSbc(:,i), MSbc(:,j),5,'filled');
        end
    end
    han=axes(figure1,'visible','off');han.Title.Visible='on';han.XLabel.Visible='off';han.YLabel.Visible='off';
    title('Scatter plot of each pair of mainshock variables (after Box-Cox transformation)');
    figurename = 'MSbcmultivariatenormal.png'; saveas(figure1,figurename);
    %}
    MSbccov = zeros(msL, msL, length(farray)-1);
    MSbccor = zeros(msL, msL, length(farray)-1);
    for i = 1:length(farray)-1
        MSbccov(:,:,i) = cov(MSbc(1:i+1,:));
        MSbccor(:,:,i) = corrcov(MSbccov(:,:,i)); %Use matlab built in function to calculate correlation matrix
    end

    MSbccoravg = zeros(msL,msL,length(farray)-1);
    MSbccorvar = zeros(msL,msL,length(farray)-1);
    for i = 1:length(farray)-1
        MSbccoravg(:,:,i) = mean(MSbccor(:,:,1:i),3);
        MSbccorvar(:,:,i) = var(MSbccor(:,:,1:i),0, 3);
    end
    %{
    figure1 = figure(f);f=f+1;
    set(gcf,'Position',get(0,'Screensize'));
    for i = 1:msL
        for j = 1:msL
            temp = zeros(length(farray)-1,1);
            subplot(msL,msL, msL*i-(msL-j));
            temp(:,1) = MSbccoravg(i,j,:);
            plot(temp);
        end
    end
    han=axes(figure1,'visible','off');han.XLabel.Visible='on';han.YLabel.Visible='on';
    ylabel(han,'mainshock (after Box-Cox transformation) parameters: averages of correlation');xlabel(han,'increasing number of simulations');
    figurename = 'MSbccoravg.png'; saveas(figure1,figurename);

    figure1 = figure(f);f=f+1;
    set(gcf,'Position',get(0,'Screensize'));
    for i = 1:msL
        for j = 1:msL
            temp = zeros(length(farray)-1,1);
            subplot(msL,msL, msL*i-(msL-j));
            temp(:,1) = MSbccorvar(i,j,:);
            plot(temp);
        end
    end
    han=axes(figure1,'visible','off');han.XLabel.Visible='on';han.YLabel.Visible='on';
    ylabel(han,'mainshock (after Box-Cox transformation) parameters: variances of correlation');xlabel(han,'increasing number of simulations');
    figurename = 'MSbccorvar.png'; saveas(figure1,figurename);
    
    MSbccorstd = MSbccorvar(:,:,end).^0.5;
    %}
    %
    iMSbccov = inv(MSbccov(:,:,end));
    for j = 1:length(farray)
        for i = 1:length(farray)
            temp = MSbc(j,:)-MSbc(i,:);
            %msbcdist(j,i) = temp*inv(MSbccov(:,:,end))*temp';
            msbcdist(j,i) = temp*iMSbccov*temp';
        end
    end
    msbcZ = linkage( msbcdist); 
    %figure1=figure(f);f=f+1;dendrogram(msbcZ);
    %figurename = 'MSbclinkage.png'; saveas(figure1,figurename);
    % Cluster based on msZ
    I = inconsistent(msbcZ);
    %
    msbcC = cluster(msbcZ, 'cutoff', 0.8); %Maxclust, maximum number of clusters
    
    %{
    % Check mainshock clustering based solely on strain drop
    temp = pdist(Mshock(:,10));
    msstrainZ = linkage( temp ); 
    %figure(f);f=f+1;dendrogram(msstrainZ);
    msstrainC = cluster(msstrainZ, 'Maxclust', 4);
    % Check mainshock clustering based solely on strain stress (disregard propagation)
    MSss = Mshock(:,9:10);
    MSsscov = cov(MSss);
    for j = 1:length(farray)
        for i = 1:length(farray)
            temp = MSss(j,:)-MSss(i,:);
            MSssmahal(j,i) = temp*inv(MSsscov)*temp';
        end
    end
    msssZ = linkage( MSssmahal ); 
    %figure(f);f=f+1;dendrogram(msssZ); 
    msssC = cluster(msssZ, 'Maxclust', 4);
    % Check mainshock clustering based solely on AE
    MSae = Mshock(:,11:12);
    MSaecov = cov(MSae);
    for j = 1:length(farray)
        for i = 1:length(farray)
            temp = MSae(j,:)-MSae(i,:);
            MSaemahal(j,i) = temp*inv(MSaecov)*temp';
        end
    end
    msaeZ = linkage( MSaemahal ); 
    %figure(f);f=f+1;dendrogram(msaeZ);
    msaeC = cluster(msaeZ, 'Maxclust', 4);
    %}
    
    MSssae = MSbc(:,9:12);
    MSssaecov = cov(MSssae);
    iMSssaecov = inv(MSssaecov); 
    for j = 1:length(farray)
        for i = 1:length(farray)
            temp = MSssae(j,:)-MSssae(i,:);
            MSssaemahal(j,i) = temp*iMSssaecov*temp';
        end
    end
    msssaeZ = linkage( MSssaemahal ); 
%     figure(f);f=f+1;dendrogram(msssaeZ); 
%     msssaeC = cluster(msssZ, 'Maxclust', 4);
    
    % (keyword MF) Foreshock and mainshock relationship - part I. Check if the two partitions based on FS and MS match each other
    % cluster based on maxclust
    % Find percentage change (decrease) in SSE with increasing number of
    % clusters. see https://hlab.stanford.edu/brian/number_of_clusters_.html It
    % is found at 8 clusters and 11 clusters the most drastic change of dSSE,
    % so either 8 clusters or 11 clusters for 30 sets of simulations.
    clear numC
    n = floor(length(farray)/3)-3 ;%length(farray); % Mod ri_opt given varying N_cluster range
    msSSE = zeros(n, 1);
    for i = 1:n
        numC = i+3;%i % Mod ri_opt given varying N_cluster range
        msbcC_elbow = cluster(msbcZ, 'Maxclust', numC);
        tempmahaltotal = zeros(numC, 1);
        for j = 1:numC %given a cluster number, foreach cluster. Say 5 clusters, for j = first, second,...fifth clusters
            temp = find(msbcC_elbow == j);
            tempmean = mean( MSbc(temp,:), 1);
            templength = length(temp);
            tempmahal = zeros(templength, 1);
            for jj = 1:templength
                %find distance between each Mshock(temp) and tempmean
                tempsub = MSbc(temp(jj),:) - tempmean;
                tempmahal(jj) = pdist2(MSbc(temp(jj),:),tempmean,'euclidean');
                %tempmahal(jj) = tempsub * iMSbccov * tempsub';
            end
            tempmahaltotal(j) = sum(tempmahal); %distances from each MS to the centroid of each cluster, given number of clusters. For instance, 5 clusters, 5 distances
        end
        msSSE(i) = sum(tempmahaltotal);
    end
    temp = flip(msSSE);
    %msSSE_output{iii-ijk+1} = msSSE;
    
    % msSSE_changepercentage = ( (temp(2:end)-temp(1:end-1))./temp(1:end-1) ) * 100;
    % figure(f);f=f+1; scatter(n-1:-1:1, msSSE_changepercentage);
    msSSE_changepercentage = abs( (msSSE(2:end) - msSSE(1:end-1))./msSSE(1:end-1) )*100;
    %msSSE_changepercentage_output{iii-ijk+1} = msSSE_changepercentage;
    %figure(f);f=f+1; scatter(1:n-1, msSSE_changepercentage);
    ddmsSSE = msSSE_changepercentage(2:end) - msSSE_changepercentage(1:end-1);
    numC_opt(iii) = find(ddmsSSE == min(ddmsSSE)) + 1+3;%numC_opt(iii-ijk+1) = find(ddmsSSE == min(ddmsSSE)) + 1; % Mod ri_opt given varying N_cluster range

%     % Foreshock and mainshock relationship - part I.1
%     % With increasing amount of simulations, check if rand index between two
%     % clusters converges at various number of clusters
%     numC = 1:n;%16;
%     ms_sil = zeros(length(numC),1);
%     ms_gap = zeros(length(numC),1);
%     %figure(f);f=f+1;
%     ah1 = ceil(sqrt(length(numC)));
%     for j = 1:length(numC) %with increasing number of clusters
%         ri = zeros(length(farray)-numC(j) , 1);
%         for i = numC(j):length(farray) %Given a number of clusters, keep increasing number of simulations and calculate rand index with respect of number of simulations. Check to see if ri converges
%             FShighmahaltr_stack = FShighmahaltr(1:i,1:i);
%             fsZ_stack = linkage(FShighmahaltr_stack);
%             fsbcC_stack = cluster(fsZ_stack, 'Maxclust', numC(j));
% 
%             msbcdist_stack = msbcdist(1:i,1:i);
%             msbcZ_stack = linkage(msbcdist_stack);
%             msbcC_stack = cluster(msbcZ_stack, 'Maxclust', numC(j));
% 
%             ri(i-(numC(j)-1)) = rand_index(fsbcC_stack, msbcC_stack);
%         end
%         if numC(j) == numC_opt(iii)%numC_opt(iii-ijk+1)
%             ri_opt(iii) = ri(end);%ri_opt(iii-ijk+1) = ri(end);
%         end
%     %     subplot(ah1, ah1, j); %subfigure number indicates various number of clusters. In each figure, horizontal axis indicates increasing number of simulations used to calculate ri 
%     %     plot(ri);
%     end
    
    % (keyword MF) Foreshock and mainshock relationship - part I.2. Num cluster based on FS. Check if the two partitions based on FS and MS match each other
    clear numC
    n = floor(length(farray)/3)-3 ;%16;%length(farray)-1; % Mod ri_opt given varying N_cluster range
    fsSSE = zeros(n, 1);
    for i = 1:n
        numC = i+3; %i; % Mod ri_opt given varying N_cluster range
        fsbcC_elbow = cluster(fsZ, 'Maxclust', numC);
        tempmahaltotal = zeros(numC, 1);
        for j = 1:numC %given a cluster number, foreach cluster. Say 5 clusters, for j = first, second,...fifth clusters
            temp = find(fsbcC_elbow == j);
            tempmean = mean(FShightotalbc(:,:,temp),3);
            templength = length(temp);
            tempmahal = zeros(templength, 1);
            for jj = 1:templength
                tempsub = FShightotalbc(:,:,temp(jj)) - tempmean;
                tempmahal(jj) = trace(tempsub * iFShtbccov  * tempsub'); 
            end
            tempmahaltotal(j) = sum(tempmahal); %distances from each MS to the centroid of each cluster, given number of clusters. For instance, 5 clusters, 5 distances
        end
        fsSSE(i) = sum(tempmahaltotal);
    end
    temp = flip(fsSSE);
    fsSSE_changepercentage = abs( (fsSSE(2:end) - fsSSE(1:end-1))./fsSSE(1:end-1) )*100;
    ddfsSSE = fsSSE_changepercentage(2:end) - fsSSE_changepercentage(1:end-1);
    numC_opt_fs(iii) = find(ddfsSSE == min(ddfsSSE)) + 1+3;%numC_opt_fs(iii-ijk+1) = find(ddfsSSE == min(ddfsSSE)) + 1; % Mod ri_opt given varying N_cluster range
    %figure(f);f=f+1; scatter(1:n-1, fsSSE_changepercentage);
    %
    numC = 1:n;%16;
    %figure(f);f=f+1;
    %ah1 = ceil(sqrt(length(numC)));
    for j = 1:length(numC) %with increasing number of clusters
        ri = zeros(length(farray)-numC(j) , 1);
        for i = numC(j):length(farray) %Given a number of clusters, keep increasing number of simulations and calculate rand index with respect of number of simulations. Check to see if ri converges
            FShighmahaltr_stack = FShighmahaltr(1:i,1:i);
            fsZ_stack = linkage(FShighmahaltr_stack);
            fsbcC_stack = cluster(fsZ_stack, 'Maxclust', numC(j));

            msbcdist_stack = msbcdist(1:i,1:i);
            msbcZ_stack = linkage(msbcdist_stack);
            msbcC_stack = cluster(msbcZ_stack, 'Maxclust', numC(j));

            ri(i-(numC(j)-1)) = rand_index(fsbcC_stack, msbcC_stack);
        end
        
        if numC(j) == numC_opt(iii)%numC_opt(iii-ijk+1)
            ri_opt(iii) = ri(end);%ri_opt(iii-ijk+1) = ri(end);
        end
        
        if numC(j) == numC_opt_fs(iii)%numC_opt_fs(iii-ijk+1)
            ri_opt_fs(iii) = ri(end);%ri_opt_fs(iii-ijk+1) = ri_fs(end);
        end
        %subplot(ah1, ah1, j); %subfigure number indicates various number of clusters. In each figure, horizontal axis indicates increasing number of simulations used to calculate ri 
        %plot(ri_fs);
    end
    
    fsbcC_opt = cluster(fsZ, 'Maxclust', numC_opt_fs(iii));%numC_opt_fs(iii-ijk+1));
    msbcC_opt = cluster(msbcZ, 'Maxclust', numC_opt(iii));%numC_opt(iii-ijk+1));
    ri_opt_fsms(iii) = rand_index(fsbcC_opt, msbcC_opt);%ri_opt_fsms(iii-ijk+1) = rand_index(fsbcC_opt, msbcC_opt);
    
end


%%
figure(f);f=f+1;ah1 = ceil(sqrt(length(msSSE_changepercentage_output)));
for i = 1:length(msSSE_changepercentage_output)
    subplot(ah1, ah1, i); plot(msSSE_changepercentage_output{i});
end
%%
figure(f);f=f+1;ah1 = ceil(sqrt(length(msSSE_output)));
for i = 1:length(msSSE_changepercentage_output)
    subplot(ah1, ah1, i); plot(msSSE_output{i});
end


