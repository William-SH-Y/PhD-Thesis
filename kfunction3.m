%function [K, DIST, W, Dop, rbox, ori] = kfunction(dataXY,xK,box,method)
function K = kfunction3(dataXY,dataXY_stretch,xK,box,method,mark,markweighted, EI) %mark weights based on April 16, 2020 derived results
    % KFUNCTION calculates Ripleys K function 
    % K = kfunction(dataXY,xK,box,method) - returns vector K containing value
    % of Ripley's K-function of dataXY in the distances in xK.
    % dataXY - N-by-2 vector where N is number of datapoints. Each row
    % corresponds to x and y coordinates of each datapoint
    % xK - corresponds to the distances where K function should be computed.
    % K is the same size as xK...
    % box - rectangular boudnary of the data: box = [xlim1, xlim2, ylim1,
    % ylim2]
    % method - switch between edge correction. If method=0, no edge correction
    % is applied. If method=1, datapoint is used for estimation of K(h) only if
    % it is at least h units away from the box (aka rbox>xK(k))


    if nargin<4 method=1; end
    [N,k] = size(dataXY);
    if k~=2 error('dataXY must have two columns'); end

%     rbox = min([    dataXY(:,1)'-box(1);
%                     box(2)-dataXY(:,1)';
%                     dataXY(:,2)'-box(3); 
%                     box(4)-dataXY(:,2)']);
    % rbox is the nearest distance of each datapoint to the box. It's
    % calculating the distance from the point to each side of the box, and
    % find the minimum one in each four, so you end up with 4 values.
    % In a circular box, it will be ( radius minus (distance between the origin
    % and the point)).    
    radius = abs((box(2) - box(1)) / 2);
    ori = [ mean([box(1), box(2)]), mean([box(3), box(4)]) ];
    Dop = sqrt((dataXY(:,1)-ori(1)).^2 + (dataXY(:,2)-ori(2)).^2); %distance from the point to the origin of the study area
    rbox = radius - Dop;
    
    mu = mean(mark);
    mi = meshgrid(mark, mark);
    mj = mi';

    DIST = squareform(pdist(dataXY,'euclidean')); % calculated based on squeezed-in-x-direction coordinates
    DIST_stretch = squareform(pdist(dataXY_stretch,'euclidean'));
    DIST_E = DIST(:,:)*EI;

    if method==1 % edge correction...
    K = zeros(length(xK),1);
    Nk = length(K);
    %wb = waitbar(0,'Computing Ripley''s K-function...');
   % for k=1:Nk
        %waitbar(k/Nk,wb);    
%         I = find(rbox>xK(k)); %If xK is larger than rbox, the K function is not calculated
%         if ~isempty(I)
%             K(k) = sum(sum(DIST(2:end,I)<xK(k)))/length(I);
%         end
        %Instead of just ignoring all the K for rbox>xK(k), I should apply
        %a weight function instead
        
        for k=1:length(K)            
            W = zeros(N, N); %weight factor
            for i = 1:N
                for j = 1:N
                    if i == j
                        W(i, j) = 0;
                    else
                        if DIST(i,j) <= rbox(i)
                            W(i, j) = 1;
                        elseif DIST(i, j) > rbox(i)
%                             if DIST(i,j) <= Dop(i)                              
                                %W(i, j) = 1/( acos( (Dop(i)^2+DIST(i,j)^2-radius^2) / (2*Dop(i)*DIST(i,j)) ) / pi );
                                beta = acos((Dop(i)^2 + DIST(i,j)^2 - radius^2)/(2*Dop(i)*DIST(i,j)));
                                %alpha = pi - beta;
                                W(i,j) = pi/beta;
%                             elseif DIST(i,j) > Dop(i)
%                                 beta = acos((radius^2+Dop(i)^2-DIST(i,j)^2)/(2*radius*Dop(i)));
%                                 alpha = pi - beta;
%                                 W(i,j) = pi/alpha;
%                             end
                        end
                    end
                end
            end
            
            temp = (DIST_E(:,:)<xK(k));
            if markweighted == 0 %Mark-weighted K function
                K(k) = sum(sum(...
                    temp.*W ...
                    ))/N; %DIST(2:end,:)<xK(k) is essentially the I(x) indicator function. (Ripley’s K function, Philip M. Dixon, Volume 3, pp 1796–1803)
            elseif markweighted == 1
                K(k) = sum(sum(...
                temp.*W.*mi.*mj ...
                ))/N/(mu^2); %mi mj are weights for i and j; mu is the average
            end
        end
    %end
    %close (wb);

    elseif method==0 % no correction
        K = zeros(length(xK),1);
        for k=1:length(K)
            %K(k) = sum(sum(DIST(2:end,:)<xK(k)))/N; %DIST(2:end,:)<xK(k) is essentially the I(x) indicator function. (Ripley’s K function, Philip M. Dixon, Volume 3, pp 1796–1803)
            if markweighted == 0
                K(k) = ( sum(sum(DIST(:,:)<xK(k))) - N )/N;
            elseif markweighted == 1
                K(k) = ( sum(sum(DIST(:,:)<xK(k)).*mi.*mj) - N*mu^2 )/N/(mu^2);
            end
        end
    end

    % fprintf ('\b'); fprintf ('\b'); fprintf ('\b'); fprintf ('\b');
    % fprintf ('100%%\n');

    %lambda = N/((box(2)-box(1))*(box(4)-box(3))); %lambda is an estimation
    %of number over area
    if markweighted == 0
        lambda = N/(pi*radius^2);
    elseif markweighted == 1
        lambda = N^2 /(pi*radius^2);
    end
    lambda = N/(2*pi*radius^2);
        K = K/lambda;

end