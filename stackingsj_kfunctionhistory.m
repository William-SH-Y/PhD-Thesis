before_after = 1;
f = 1;
%coordinate transformation
fault_angle = pi/6;
coord_trans = [sin(fault_angle), 0, -cos(fault_angle); ...
    0, 1, 0; ...
     cos(fault_angle), 0, sin(fault_angle)];

%ellipse correction factor
theta = linspace(0, 2*pi, 200);
phi = atan(0.5*tan(theta));
dtheta = theta(2) - theta(1);
dS = sqrt( ( 4*(sin(theta(2:end))).^2 ) + ( (cos(theta(2:end))).^2 ) ) * dtheta;
Sprob = dS/sum(dS);
% figure(f);f=f+1;plot(theta(2:end),Sprob);
Sval = sqrt( 4*(cos(theta)).^2 + ((sin(theta)).^2)    );
% figure(f);f=f+1;plot(theta(2:end),Sval(2:end));
EI = sum(Sprob.*Sval(2:end));
%
mv_W = 40e-3;
mv_H = 107e-3;
ns = 14;
box(1,:) = [-mv_W/2, mv_W/2, -mv_W/2, mv_W/2];
box(2,:) = [0, mv_W/2, 0, mv_W/2];
box(3,:) = [-mv_W/2, 0, 0, mv_W/2];
box(4,:) = [-mv_W/2, 0, -mv_W/2, 0];
box(5,:) = [0, mv_W/2, -mv_W/2, 0];
box(6,:) = [0, mv_W/2, -mv_W/4, mv_W/4];
box(7,:) = [-mv_W/4, mv_W/4, 0, mv_W/2];
box(8,:) = [-mv_W/2, 0, -mv_W/4, mv_W/4];
box(9,:) = [-mv_W/4, mv_W/4, -mv_W/2, 0];
box(10,:) = [-mv_W/4, mv_W/4, -mv_W/4, mv_W/4];
box(11,:) = [0, mv_W/4, 0, mv_W/4];
box(12,:) = [-mv_W/4,0, 0, mv_W/4];
box(13,:) = [-mv_W/4, 0 , -mv_W/4, 0];
box(14,:) = [0, mv_W/4, -mv_W/4, 0];
ori = [ mean([box(:,1),box(:,2)],2), mean([box(:,3), box(:,4)],2) ];
boxrad = (box(:,2) - box(:,1))/2;

% xK = 0:1e-3:40e-3;
% dxK = xK(2) - xK(1);
% L = size(xK, 1);
% Kref = (pi*xK.^2)';
xKL = 40;
xK = zeros(xKL, 14); dxK = zeros(xKL, 1);
for i = 1:ns
    xK(:, i) = linspace(0, boxrad(i)*2, xKL)';
    dxK(i) = xK(2) - xK(1);
end
L = size(xK, 1);
Kref = (pi*xK.^2);
%
writewhen = zeros(26,1);
writewhen(1) = 7e5;
for i = 2:26
    writewhen(i) = writewhen(i-1) + 2e4;
end

mw = 0;
%
load('sd_happen', 'sd_happen'); %when stress drop happen at and after stickslip. Must run Stackingdata_cleanV1, section keyword SJ, first!

farray = {'201','208','214','215','216','220','221','222','223','224','225', '226','227','228'...
    '229','230','231','232','233','234','235', '236','237','238','239','240','241'};

%
n = length(farray);%1;
sd_surround = cell(n , 1);
for i = 1:n %for each simulation
    temparray = zeros(length(sd_happen{i}(:)) , 2);
    for j = 1:length(sd_happen{i})
        temp = sd_happen{i}(j);
        temp1 = find(writewhen < temp, 1, 'last' ); %biggest writewhen less than the j-th stress drop cycle, or writewhen right before the j-th stress drop cycle
        temp2 = find(writewhen >= temp, 1 ); %writewhen right before the j-th stress drop cycle
        temparray = [temp1 temp2];
        sd_surround{i}(j,:) = temparray;
    end
end
%
n = length(farray);
num_sd = 1;
for i = 1:n
    temp = size(sd_happen{i}, 1);
    num_sd = max(num_sd, temp); %number of stress drop.
end
sj_info = cell(n, num_sd);
sj_loc = cell(n, num_sd);
sj_loc_sqz = cell(n, num_sd);
sj_loc_trans = cell(n, num_sd);
for k = 1:n
    %sj_info_surround = cell(size(sd_surround{i}, 1) , 1);
    for j = 1:size(sd_surround{k}, 1)
        jj = sd_surround{k}(j,before_after);
        filename = fopen(strcat('k function test/',farray{k},'/sj',num2str(jj),'.txt'));
        C = textscan(filename, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f');
        sj_temp = cell2mat(C);
        sj_info{k,j} = sj_temp;
        %sj_info_surround{j} = sj_temp;
        %sj_loc{k,i}(:,1:3) = sj_info{k}(sj_info{k}(:,14)==3,1:3,i);
    end
    
    
      
%{  
    filename = fopen(strcat('k function test/',farray{k},'/sj1.txt'));
    C = textscan(filename, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    sj_temptemp = cell2mat(C); %figure out the size for each sj#.txt file
    tempsize = size(sj_temptemp);%allocate memory
    sj_temp = zeros(tempsize(1), tempsize(2), nwrite);
    sj_temp(:,:,1) = sj_temptemp;
    for i = 2:nwrite %after the sj1.txt is read, and its file size is determined, we can import the rest with a loop
        % DOES NOT WORK! sj2.txt may have different length compared to
        % sj1.txt because some contacts might have been deleted already.
        filename = fopen(strcat('k function test/',farray{k},'/sj',num2str(i),'.txt'));
        C = textscan(filename, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f');
        sj_temp(:,:,i) = cell2mat(C);
    end
    sj_info{k} = sj_temp;
%}    
    for i = 1:size(sd_surround{k}, 1)
        sj_loc{k,i}(:,1:3) = sj_info{k,i}(sj_info{k,i}(:,14)==3,1:3);
        sj_loc_trans{k,i} = (coord_trans * sj_loc{k,i}(:,1:3)')'; %new coordinates after transformation
        sj_loc_sqz{k,i} = sj_loc_trans{k,i}(:,1:2);
        sj_loc_sqz{k,i}(:,1) = sj_loc_sqz{k,i}(:,1)/2;
    end
    

%     figure1 = figure(f); f=f+1;
%     scatter3(sj_loc_trans{k}(:,1), sj_loc_trans{k}(:,2), sj_loc_trans{k}(:,3));%%, 50, sj_kn{5});
%     view([90, 90]); xlabel('x(m)'); ylabel('y(m)');
    %figurename = strcat(farray{1}, 'SJ.png'); saveas(figure1,figurename);
    %close all
    

    
end
%
Dop_all = cell(n, num_sd, ns);
for i = 1:n
    for k = 1:size(sd_surround{i}, 1)
    for j = 1:ns
        Dop_all{i,k,j} = pdist2([sj_loc_sqz{i,k}(:,1), sj_loc_sqz{i,k}(:,2)], ori(j,:)); %all the points in simulation i that are within the radius of each jth center
    end
    end
end

%
dataXY = cell(n, num_sd, ns);
dataXY_stretch = cell(n, num_sd, ns);
K = zeros(size(xK, 1), n, num_sd, ns);
LK = zeros(size(xK, 1), n, num_sd, ns);
Kxcorr = zeros(2*size(xK, 1)-1, n, num_sd, ns);
Kxcorrmax = zeros(n, num_sd, ns);
for j = 1:ns %j controls which patch of the fault
    for i = 1:n %i controls which simulation
        for k = 1:size(sd_surround{i}, 1)
        dataXY{i,k,j} = sj_loc_sqz{i,k}( Dop_all{i,k,j} < boxrad(j) , 1:2 );
        dataXY_stretch{i,k,j} = sj_loc_trans{i,k}( Dop_all{i,k,j} < boxrad(j) , 1:2 );
        sj_kn = sj_info{i,k}(:,5).*sj_info{i,k}(:,6); %markweight
        K(:,i,k,j) = kfunction3(dataXY{i,k,j}, dataXY_stretch{i,k,j}, xK(:,j), box(j,:), 1, ...
            sj_kn(Dop_all{i,k,j} < boxrad(j)) ...
            ,mw, EI); %mw: mark weight, 0 for no weight, 1 for weighted
        LK(:,i,k,j) = sqrt(K(:,i,k,j)/pi);
        Kxcorr(:,i,k,j) = xcorr(K(:,i,k,j), Kref(:,j));
        Kxcorrmax(i,k,j) = max(Kxcorr(:,i,k,j));
        end
    end
end
%
load('sizedataXY', 'sizedataXY');
sizePatch = zeros(n, num_sd, ns);
densityPatch = zeros(n, num_sd, ns);
for j = 1:ns
    for k = 1:num_sd
        temp = dataXY(:,k,j);
        for i = 1:n
            sizePatch(i, k, j) = size(cell2mat(temp(i)),1);
        end
        densityPatch(:,k,j) = sizePatch(:, k, j)./sizedataXY(:,j);
    end
end
%%
save('densityPatch','densityPatch');
% figure(f); f=f+1;
% scatter3(sj_loc(:,1), sj_loc(:,2), sj_loc(:,3));%%, 50, sj_kn{5});
% %     colormap(jet);
% %     colorbar;
% view([90, 30]);
% xlabel('x(mm)');
% ylabel('y(mm)');

%% (keyword: plt1) Plot K for the first three simulations, with subplot showing fault heterogeneity prior to each stress drop (including mainshock and aftershocks)
% sampled over the entire fault
a = 1;
for i = 1:3%n %just test the first three simulations
    figure(f);f=f+1; ah1=ceil(sqrt(size(sd_surround{i}, 1)));
    for j = 1:size(sd_surround{i}, 1) %before each stress drop
        subplot(ah1, ah1, j);
        plot(xK(1:L, a), K(1:L,i,j,a)-Kref(:,a),xK(1:L,a), Kref(:,a)-Kref(:,a), 'k');
    end
end
% Plot LK
for i = 1:3%n %just test the first three simulations
    figure(f);f=f+1; ah1=ceil(sqrt(size(sd_surround{i}, 1)));
    for j = 1:size(sd_surround{i}, 1) %before each stress drop
        subplot(ah1, ah1, j);
        plot(xK(1:L, a), LK(1:L,i,j,a), xK(1:L, a),xK(1:L, a),'k');
    end
end
%% (keyword: plt5) Plot K for the first three simulations, with subplot showing fault heterogeneity prior to each stress drop (including mainshock and aftershocks)
% sampled over box(5,:); 
a = 5;
for i = 1:3%n %just test the first three
    figure(f);f=f+1; ah1=ceil(sqrt(size(sd_surround{i}, 1)));
    for j = 1:size(sd_surround{i}, 1)
        subplot(ah1, ah1, j);
        plot(xK(1:L, a), K(1:L,i,j,a)-Kref(:,a),xK(1:L, a), Kref(:,a)-Kref(:,a), 'k');
    end
end
%% (keyword: plt10) Plot K for the first three simulations, with subplot showing fault heterogeneity prior to each stress drop (including mainshock and aftershocks)
% sampled over box(5,:); 
a = 10;
for i = 1:3%n %just test the first three
    figure(f);f=f+1; ah1=ceil(sqrt(size(sd_surround{i}, 1)));
    for j = 1:size(sd_surround{i}, 1)
        subplot(ah1, ah1, j);
        plot(xK(1:L, a), K(1:L,i,j,a)-Kref(:,a),xK(1:L, a), Kref(:,a)-Kref(:,a), 'k');
    end
end
%% (keyword: int) find zeros for K sampled over entire fault, for various files. Integrate the area of each section of K-Kref above or below 0
Kresidual = zeros(size(K,1),size(K,2),size(K,3),size(K,4));
for i = 1:ns
    Kresidual(:,:,:,i) = K(:,:,:,i) - Kref(:,i);
end
a = 1;
Kzeros = cell(n, num_sd, ns);
Kresidualsum = cell(n,num_sd, ns); %integrate the each area above and below 0 for each K-Kref. 
xKzeros = cell(n,num_sd, ns);
for a = 1:ns
    for i = 1:n
        for j = 1:size(sd_surround{i}, 1)
            Kzeros{i,j,a} = find((Kresidual(1:end-1,i,j,a).*Kresidual(2:end, i,j, a)) < 0);
            temp = [1;Kzeros{i,j,a};L]; %boundaries where K-Kref goes between 0
            temp2 = zeros(length(Kzeros{i,j,a})+1, 1);
            for k = 2:length(Kzeros{i,j,a})+2 %within each pair of boundaries
                temp2(k-1) = sum(Kresidual(temp(k-1):temp(k),i,j,a)) * dxK(a); %integrate the each zone above and below 0
            end
            Kresidualsum{i,j,a} = temp2;
            xKzeros{i,j,a} = xK(temp);
        end  
    end
end
% (keyword: r squared)
%K dimension order: number of simulation, stress drop number (1 for stick slip, larger than one for aftershocks),sample patch number (box number)
LKrsq = zeros(n, num_sd, ns);
for k = 1:ns
    for j = 1:num_sd
        for i = 1:n
            temp = corrcoef(LK(:,i,j,k),xK(:,k));
            LKrsq(i,j,k) = temp(1,2);
        end
    end
end
save('LKrsq','LKrsq');
%%
% figure(f); f=f+1;
% for i = 1:ns
%     subplot(5, 3, i);
%     plot(xK, K(:,1,i)-Kref, 'b', xK, K(:,2,i)-Kref, 'r', xK, K(:,3,i)-Kref, 'g',xK,  K(:,4,i)-Kref, 'k');
%     xlabel('r(m)');
%     ylabel('Ripley K function L estimator');
%     title('Ripley K function L estimator for each simulation; each subfigure represents one fault patch');
% %     hold on
% %     %plot(xK, zeros(size(xK, 1),1));
% %     plot(xK, pi*xK.^2);
% %     hold off
% end
%%
Kaverage = zeros(L, n);
for i = 1:n %first dimension is the K function values; 2nd dimension, i, controls which simulation; 3 dimension controls fault patch.
    Kaverage(:,i) = mean( K(:,i,2:10), 3 );
end
%% Plot K for the first file, sampled over entire fault, before stickslip
figure1 = figure(f);f=f+1;plot(xK(1:L), K(1:L,1,1,1)-Kref(:,1),xK(1:L), Kref(:,1)-Kref(:,1), 'k');
xlabel('xK(m)'); ylabel('K');%figurename = strcat(farray{1}, 'K.png'); saveas(figure1,figurename);

%% Plot K for the various file, sampled over entire fault
figure1 = figure(f);f=f+1;
ah1 = ceil(sqrt(n));
pltl = size(xK, 1);
for i = 1:n
    subplot(ah1, ah1, i); 
    plot(xK(1:pltl), K(1:pltl,i,1)-Kref(1:pltl),xK(1:pltl),  Kref(1:pltl)-Kref(1:pltl), 'k');
    %plot(xK(1:pltl), K(1:pltl,i,1)-Kref(1:pltl),xK(1:pltl), zeros(1,pltl), 'k');
    xlabel('r(m)'); ylabel('K');
end
figurename = strcat('K.png'); saveas(figure1,figurename);
%% Plot K for the various file, average of samples over various fault patches
figure(f);f=f+1;
ah1 = ceil(sqrt(n));
for i = 1:n
    subplot(ah1, ah1, i); plot(xK, Kaverage(:,i)-Kref, xK, Kref-Kref);
    xlabel('r(m)'); ylabel('Ripley K function');
end
%
% figure(f); f=f+1;
% subplot(1, 2, 1);
% plot(xK, K(:,1,1)-Kref,xK, Kref-Kref, 'k');%, xK, Kref);
% %plot(xK, K(:,1,1), xK, Kref);
% xlabel('r(m)');
% ylabel('Ripley K function');
% title('Ripley K function average over entire fault');
% subplot(1, 2, 2);
% plot(xK, Kaverage(:,1)-Kref, xK, Kref-Kref);%, xK, Kref);
% xlabel('r(m)');
% ylabel('Ripley K function');
% title('Ripley K function averaged over various sample fault patch');


%
% Kaverage = zeros(L, n);
% for i = 1:n %first dimension is the K function values; 2nd dimension, i, controls which simulation; 3 dimension controls fault patch.
%     Kaverage(:,i) = mean( K(:,i,2:10), 3 );
% end
% figure(f); f=f+1;
% 
% subplot(1, 2, 1);
% %plot(xK, K(:,1,1)-Kref, 'b',xK, K(:,2,1)-Kref, 'r',xK, K(:,3,1)-Kref, 'g',xK, K(:,4,1)-Kref, 'k');%, xK, Kref);
% plot(xK, K(:,1,1), 'b',xK, K(:,6,1), 'r',xK, K(:,7,1), 'k', xK, Kref);
% xlabel('r(m)');
% ylabel('Ripley K function');
% title('Ripley K function average over entire fault');
% % subplot(1, 2, 2);
% % plot(xK, Kaverage(:,1)-Kref, 'b',xK, Kaverage(:,2)-Kref, 'r',xK, Kaverage(:,3)-Kref, 'g',xK, Kaverage(:,4)-Kref, 'k');%, xK, Kref);
% % xlabel('r(m)');
% % ylabel('Ripley K function');
% % title('Ripley K function averaged over various sample fault patch');








