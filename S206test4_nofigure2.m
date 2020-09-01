for input = 1:27  
    run(input);
    close all
    clear
    clc
end

function [varargout] = run(input)
farray = {'201','208','214','215','216','220','221','222','223','224','225', '226','227','228',...
    '229','230','231','232','233','234','235', '236','237','238','239','240','241','242', '243', '244',...
    '245','246','247','248','249','250','251','254','255'};
sim_workingonthis = farray{input};
n_simulation = length(farray);
% sim_workingonthis = '231';
n_workingonthis = find(strcmp(farray,sim_workingonthis));
fnumber = farray{n_workingonthis};

%fnumber = '201'; %Non dynamic weakening. fa = 30
%fnumber = '208';
%fnumber = '214';
%fnumber = '215';
%fnumber = '216';
%fnumber = '220';
%fnumber = '221';
%fnumber = '222';
%fnumber = '223';
%fnumber = '224';
%fnumber = '225';
%fnumber = '226';
%fnumber = '227';
%fnumber = '228';
%fnumber = '229';
%fnumber = '230';
%fnumber = '231';
%fnumber = '232';
%fnumber = '233';
%fnumber = '234';
%fnumber = '235';
%fnumber = '236';
%fnumber = '237';
%fnumber = '238';
%fnumber = '239';
%fnumber = '240';
%fnumber = '241';

sh = csvread(strcat('S',fnumber,'-sh.csv'),1);
shh = sh(1:(end-1), 2);

ss = csvread(strcat('S',fnumber,'-ss.csv'), 1);
ss = ss(1:(end-1), :);

%num_sensor = 30;
if strcmp(fnumber, '224')==1 || strcmp(fnumber, '225')==1 || strcmp(fnumber,'226')==1 || strcmp(fnumber,'227')==1 ...
        || strcmp(fnumber,'228')==1 || strcmp(fnumber,'229')==1 || strcmp(fnumber,'230')==1 || strcmp(fnumber,'231') == 1 ...
        || strcmp(fnumber,'232') == 1 || strcmp(fnumber,'233') == 1 || strcmp(fnumber,'234') == 1 || strcmp(fnumber,'235') == 1 ...
        || strcmp(fnumber,'236') == 1 || strcmp(fnumber,'237') == 1 || strcmp(fnumber,'238') == 1 || strcmp(fnumber,'239') == 1 ...
        || strcmp(fnumber,'240') == 1 || strcmp(fnumber,'241') == 1 || strcmp(fnumber,'242') == 1|| strcmp(fnumber,'243') == 1 ...
        || strcmp(fnumber,'244') == 1 || strcmp(fnumber,'245') == 1 || strcmp(fnumber,'246') == 1 || strcmp(fnumber,'247') == 1 ...
        || strcmp(fnumber,'248') == 1 || strcmp(fnumber,'249') == 1 || strcmp(fnumber,'250') == 1 || strcmp(fnumber,'251') == 1 ...
        || strcmp(fnumber,'254') == 1 || strcmp(fnumber,'255') == 1
    num_sensor = 80;
else
    num_sensor = 30;
end
num_rec = 6e4;
num_record = length(shh); %technically it is the total nonzero elements of bv, but here for simplicity sake (since bv is not read), we temporarily set it to total number of recordings per trial

if strcmp(fnumber, '201') == 1
    ki = 5:9;
elseif strcmp(fnumber, '202') == 1
    ki = 2:2;
elseif strcmp(fnumber, '208') == 1
    ki = [1 3:7];
elseif strcmp(fnumber, '214') == 1 || strcmp(fnumber, '216')==1
    ki = 1:6;
elseif strcmp(fnumber, '215') == 1 || strcmp(fnumber,'220')==1 || strcmp(fnumber,'221')==1 || strcmp(fnumber,'222')==1
    ki = 1:5;
elseif strcmp(fnumber, '223') ==1
    ki = 1:4;
else
    ki = 1:1;
end
num_trials = length(ki);

cyc_offset_array = zeros(n_simulation, 1);
cyc_offset_array(1) = 1.31e4;
cyc_offset_array(2) = 18200;
cyc_offset_array(3) = 29800;
cyc_offset_array(4) = 21880;
cyc_offset_array(5) = 19520;
cyc_offset_array(6) = 18040;
cyc_offset_array(7) = 13540;
cyc_offset_array(8) = 13640;
cyc_offset_array(9) = 11740;
cyc_offset_array(10) = 11920;
cyc_offset_array(11) = 12600;
cyc_offset_array(12) = 13520;
cyc_offset_array(13) = 12740;
cyc_offset_array(14) = 10900;
cyc_offset_array(15) = 11880;
cyc_offset_array(16) = 10180;
cyc_offset_array(17) = 12120;
cyc_offset_array(18) = 17980;
cyc_offset_array(19) = 10400;
cyc_offset_array(20) = 11920;
cyc_offset_array(21) = 10940;
cyc_offset_array(22) = 12780;
cyc_offset_array(23) = 10780;
cyc_offset_array(24) = 15540;
cyc_offset_array(25) = 10540;
cyc_offset_array(26) = 15260;
cyc_offset_array(27) = 10420;
cyc_offset_array(28) = 17560;
cyc_offset_array(29) = 11180;
cyc_offset_array(30) = 14020;
cyc_offset_array(31) = 11000;
cyc_offset_array(32) = 14080;
cyc_offset_array(33) = 10560;
cyc_offset_array(34) = 11200;
cyc_offset_array(35) = 12400;
cyc_offset_array(36) = 11360;
cyc_offset_array(37) = 23000; %251
cyc_offset_array(38) = 12340; %254
cyc_offset_array(39) = 12560; %255
cyc_offset = cyc_offset_array(n_workingonthis);
save('cyc_offset_array','cyc_offset_array');

%{
if strcmp(fnumber, '206') == 1
    cyc_offset = 1.31e4;
elseif strcmp(fnumber, '210') == 1
    cyc_offset = 1.31e4;
elseif strcmp(fnumber, '201') == 1
    cyc_offset = 1.31e4;
elseif strcmp(fnumber, '207') == 1
    cyc_offset = 81760;
elseif strcmp(fnumber, '205') == 1
    cyc_offset = 81760;
elseif strcmp(fnumber, '202') == 1
    cyc_offset = 9520;
elseif strcmp(fnumber, '204') == 1
    cyc_offset = 18200;
elseif strcmp(fnumber, '208') == 1
    cyc_offset = 18200;
elseif strcmp(fnumber, '214') == 1
    cyc_offset = 29800;
elseif strcmp(fnumber, '215') == 1
    cyc_offset = 21880; 
elseif strcmp(fnumber, '216') == 1
    cyc_offset = 19520; 
elseif strcmp(fnumber, '220') == 1
    cyc_offset = 18040;
elseif strcmp(fnumber, '221') == 1
    cyc_offset = 13540; 
elseif strcmp(fnumber, '222') == 1
    cyc_offset = 13640;
elseif strcmp(fnumber, '223') == 1
    cyc_offset = 11740;
elseif strcmp(fnumber, '224') == 1
    cyc_offset = 11920;
elseif strcmp(fnumber, '225') == 1
    cyc_offset = 12600;
elseif strcmp(fnumber, '226') == 1
    cyc_offset = 13520;
elseif strcmp(fnumber, '227') == 1
    cyc_offset = 12740;
elseif strcmp(fnumber, '228') == 1
    cyc_offset = 10900;
elseif strcmp(fnumber, '229') == 1
    cyc_offset = 11880;
elseif strcmp(fnumber, '230') == 1
    cyc_offset = 10180;
elseif strcmp(fnumber, '231') == 1
    cyc_offset = 12120;
elseif strcmp(fnumber, '232') == 1
    cyc_offset = 17980;
elseif strcmp(fnumber, '233') == 1
    cyc_offset = 10400;
elseif strcmp(fnumber, '234') == 1
    cyc_offset = 11920;
elseif strcmp(fnumber, '235') == 1
    cyc_offset = 10940;
elseif strcmp(fnumber, '236') == 1 
    cyc_offset = 12780;
elseif strcmp(fnumber, '237') == 1
    cyc_offset = 10780;
elseif strcmp(fnumber, '238') == 1
    cyc_offset = 15540;
elseif strcmp(fnumber, '239') == 1
    cyc_offset = 10540;
elseif strcmp(fnumber, '240') == 1
    cyc_offset = 15260;
elseif strcmp(fnumber, '241') == 1
    cyc_offset = 10420;
end
%}
dt = 1.0e-7;
sample_rate = 20;
cyc = (0:1:num_record-1)* sample_rate + cyc_offset;
t = cyc * dt;

cn = zeros(num_record, num_sensor*num_trials);
cs = zeros(num_record, num_sensor*num_trials);
bdz = zeros(num_record, num_sensor*num_trials);
bvz = zeros(num_record, num_sensor*num_trials);
bdy = zeros(num_record, num_sensor*num_trials);
bvy = zeros(num_record, num_sensor*num_trials);
bdx = zeros(num_record, num_sensor*num_trials);
bvx = zeros(num_record, num_sensor*num_trials);
bdzoe = zeros(num_record, num_sensor*num_trials);
bvzoe = zeros(num_record, num_sensor*num_trials);
bdyoe = zeros(num_record, num_sensor*num_trials);
bvyoe = zeros(num_record, num_sensor*num_trials);
bdxoe = zeros(num_record, num_sensor*num_trials);
bvxoe = zeros(num_record, num_sensor*num_trials);

popcnind = 1;

for k = ki 
    fname = strcat('c_model',fnumber,'-',int2str(k),'.txt');
    filename = fopen(fname);
    c_model = textscan(filename, '%s');
    
    is_sj = zeros(num_sensor, 1);
    for i = 1:num_sensor
        if strcmp(c_model{1}(i), 'smoothjoint') == 1
            is_sj(i) = 1;
        end
    end

    isnt_sj_ind = find(is_sj == 0);
    fname = strcat('cn',fnumber,'-',int2str(k),'.txt');
    cn_temp = openfile_cncs(fname, num_sensor,num_rec, num_record);
    fname = strcat('cs',fnumber,'-',int2str(k),'.txt');
    cs_temp = openfile_cncs(fname, num_sensor,num_rec, num_record);
    fname = strcat('bvz',fnumber,'-',int2str(k),'.txt');
    [num_record, bvz_temp] = openfile(fname, num_sensor,num_rec);
    fname = strcat('bdz',fnumber,'-',int2str(k),'.txt');
    [num_record, bdz_temp] = openfile(fname, num_sensor,num_rec);
    fname = strcat('bvy',fnumber,'-',int2str(k),'.txt');
    [num_record, bvy_temp] = openfile(fname, num_sensor,num_rec);
    fname = strcat('bdy',fnumber,'-',int2str(k),'.txt');
    [num_record, bdy_temp] = openfile(fname, num_sensor,num_rec);
    fname = strcat('bvx',fnumber,'-',int2str(k),'.txt');
    [num_record, bvx_temp] = openfile(fname, num_sensor,num_rec);
    fname = strcat('bdx',fnumber,'-',int2str(k),'.txt');
    [num_record, bdx_temp] = openfile(fname, num_sensor,num_rec);
    
    fname = strcat('bvz_oe',fnumber,'-',int2str(k),'.txt');
    [num_record, bvzoe_temp] = openfile(fname, num_sensor,num_rec);
    fname = strcat('bdz_oe',fnumber,'-',int2str(k),'.txt');
    [num_record, bdzoe_temp] = openfile(fname, num_sensor,num_rec);
    fname = strcat('bvy_oe',fnumber,'-',int2str(k),'.txt');
    [num_record, bvyoe_temp] = openfile(fname, num_sensor,num_rec);
    fname = strcat('bdy_oe',fnumber,'-',int2str(k),'.txt');
    [num_record, bdyoe_temp] = openfile(fname, num_sensor,num_rec);
    fname = strcat('bvx_oe',fnumber,'-',int2str(k),'.txt');
    [num_record, bvxoe_temp] = openfile(fname, num_sensor,num_rec);
    fname = strcat('bdx_oe',fnumber,'-',int2str(k),'.txt');
    [num_record, bdxoe_temp] = openfile(fname, num_sensor,num_rec);
%     fname = strcat('bdz_oe',fnumber,'-',int2str(k),'.txt');
%     [num_record, bdz_oe_temp] = openfile(fname, num_sensor,num_rec);
    
    
    num_ball_sets = 2;
    fname = strcat('b_coord',fnumber,'-',int2str(k),'.txt');
    bcoord_temp = openfile_bcoord(fname, num_sensor, num_ball_sets);
    
    num_usefulsensor = num_sensor - length(isnt_sj_ind);
    if length(cn_temp(:, 1)) > num_usefulsensor
        bcoord_temp(isnt_sj_ind, :) = [];
        cn_temp(:, isnt_sj_ind) = [];
        cs_temp(:, isnt_sj_ind) = [];
        bvz_temp(:, isnt_sj_ind) = [];
        bdz_temp(:, isnt_sj_ind) = [];
        bvy_temp(:, isnt_sj_ind) = [];
        bdy_temp(:, isnt_sj_ind) = [];
        bvx_temp(:, isnt_sj_ind) = [];
        bdx_temp(:, isnt_sj_ind) = [];
        
        bvzoe_temp(:, isnt_sj_ind) = [];
        bdzoe_temp(:, isnt_sj_ind) = [];
        bvyoe_temp(:, isnt_sj_ind) = [];
        bdyoe_temp(:, isnt_sj_ind) = [];
        bvxoe_temp(:, isnt_sj_ind) = [];
        bdxoe_temp(:, isnt_sj_ind) = [];
%         bdz_oe_temp(:, isnt_sj_ind) = [];
    else
        ;
    end
    

    
    if strcmp(sim_workingonthis,'255')==1 
        cn(:, ( popcnind: popcnind+num_usefulsensor-1) ) = cn_temp;
        cs(:, ( popcnind: popcnind+num_usefulsensor-1) ) = cs_temp;
        bvz(1:end-1, ( popcnind: popcnind+num_usefulsensor-1) ) = bvz_temp;
        bdz(1:end-1, ( popcnind: popcnind+num_usefulsensor-1) ) = bdz_temp;
        bvy(1:end-1, ( popcnind: popcnind+num_usefulsensor-1) ) = bvy_temp;
        bdy(1:end-1, ( popcnind: popcnind+num_usefulsensor-1) ) = bdy_temp;
        bvx(1:end-1, ( popcnind: popcnind+num_usefulsensor-1) ) = bvx_temp;
        bdx(1:end-1, ( popcnind: popcnind+num_usefulsensor-1) ) = bdx_temp;
        bvz(end, ( popcnind: popcnind+num_usefulsensor-1) ) = bvz_temp(end,:);
        bdz(end, ( popcnind: popcnind+num_usefulsensor-1) ) = bdz_temp(end,:);
        bvy(end, ( popcnind: popcnind+num_usefulsensor-1) ) = bvy_temp(end,:);
        bdy(end, ( popcnind: popcnind+num_usefulsensor-1) ) = bdy_temp(end,:);
        bvx(end, ( popcnind: popcnind+num_usefulsensor-1) ) = bvx_temp(end,:);
        bdx(end, ( popcnind: popcnind+num_usefulsensor-1) ) = bdx_temp(end,:);
    else
    cn(:, ( popcnind: popcnind+num_usefulsensor-1) ) = cn_temp;
    cs(:, ( popcnind: popcnind+num_usefulsensor-1) ) = cs_temp;
    bvz(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvz_temp;
    bdz(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdz_temp;
    bvy(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvy_temp;
    bdy(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdy_temp;
    bvx(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvx_temp;
    bdx(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdx_temp;
    end
    bvzoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvzoe_temp;
    bdzoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdzoe_temp;
    bvyoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvyoe_temp;
    bdyoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdyoe_temp;
    bvxoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvxoe_temp;
    bdxoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdxoe_temp;
    
    bcoord(( popcnind: popcnind+num_usefulsensor-1), :) = bcoord_temp;
    popcnind = popcnind+num_usefulsensor;
    %cn(:, (num_sensor*(k-2)+1):(num_sensor*(k-2)+1)+num_usefulsensor-1 )= cn_temp;
    %cs(:, (num_sensor*(k-2)+1):(num_sensor*(k-2)+1)+num_usefulsensor-1 )= cs_temp;
    
    
end

num_usefulsensor = popcnind - 1;

[cn, ia, ic] = unique(cn', 'rows', 'stable');
[cs, ia2, ic2] = unique(cs', 'rows', 'stable');
[bdz, ia3, ic3] = unique(bdz', 'rows', 'stable');
[bvz, ia4, ic4] = unique(bvz', 'rows', 'stable');
[bdy, ia5, ic5] = unique(bdy', 'rows', 'stable');
[bvy, ia6, ic6] = unique(bvy', 'rows', 'stable');
[bdx, ia7, ic7] = unique(bdx', 'rows', 'stable');
[bvx, ia8, ic8] = unique(bvx', 'rows', 'stable');

%bdzoe = bdzoe(:, )

% [bdzoe, ia9, ic9] = unique(bdzoe', 'rows', 'stable');
% [bvzoe, ia10, ic10] = unique(bvzoe', 'rows', 'stable');
% [bdyoe, ia11, ic11] = unique(bdyoe', 'rows', 'stable');
% [bvyoe, ia12, ic12] = unique(bvyoe', 'rows', 'stable');
% [bdxoe, ia13, ic13] = unique(bdxoe', 'rows', 'stable');
% [bvxoe, ia14, ic14] = unique(bvxoe', 'rows', 'stable');

cn = cn';
cs = cs';
bdz = bdz';
bvz = bvz';
bdy = bdy';
bvy = bvy';
bdx = bdx';
bvx = bvx';
num_usefulsensor = size(cn, 2) - 1;

if size(bdyoe, 2) > length(ia3)
    bdzoe = bdzoe(:, ia3);
    bvzoe = bvzoe(:, ia4);
    bdyoe = bdyoe(:, ia5);
    bvyoe = bvyoe(:, ia6);
    bdxoe = bdxoe(:, ia7);
    bvxoe = bvxoe(:, ia8);
else
    ;
end

cdz = bdz - bdzoe;
cdy = bdy - bdyoe;
cdx = bdx - bdxoe;


filename = fopen(strcat('sj',fnumber,'.txt'));
C = textscan(filename, '%f %f %f %f %f %f %f %f %f %f %f %f %f');
sj_info = cell2mat(C);
for j = 1:3
    sj_loc(:, j) = sj_info(:, j);
    sj_bloc1(:, j) = sj_info(:, j+7);
    sj_bloc2(:, j) = sj_info(:, j+10);
end
sj_radius = sj_info(:, 4);
sj_area = sj_info(:, 5);
sj_kn = sj_info(:, 6).*sj_area;
sj_ks = sj_info(:, 7).*sj_area;

%Use particle coordinates to figure out contacts coordinates and
%orientations
temp = ia(1:(end-1));
bcoord = bcoord(temp, :);
ccoord = [ (bcoord(:, 1)+bcoord(:, 4))/2 , (bcoord(:, 2)+bcoord(:, 5))/2 , (bcoord(:, 3)+bcoord(:, 6))/2]; 
[ccoord] = around_decimals(ccoord, 4);
cvector = [bcoord(:, 1) - bcoord(:, 4) , ...
    bcoord(:, 2) - bcoord(:, 5) , ...
    bcoord(:, 3) - bcoord(:, 6)];
%contacts orientation 
%The fault has an orientation of 30 degrees. To figure out how much
%contact coordinates are perpendicular to the fault, project the contact
%vector to the new rotated coordinate
fault_angle = pi/6;
coord_trans = [sin(fault_angle), 0, -cos(fault_angle); ...
    0, 1, 0; ...
    cos(fault_angle), 0, sin(fault_angle)];
cvector_project = (coord_trans * cvector')';
%hist(cvector_project(:, 1));
%rearrange the contacts based on how much it is perpendicular to z' axis (aka smallest abs(x'))
cvector_project_rearrange = sortrows([abs(cvector_project), (1:1:num_usefulsensor)']);
ro = cvector_project_rearrange(:, 4); %cvector_project_recorder, give it an easy name to call in the future.

% After rearrangement of the contacts based on orientation, try plotting cn based on new order
f = 1;
%{

figure(f); f=f+1;
scatter(1:1:num_usefulsensor, cvector_project_rearrange(:, 1));
ylabel('contact vector projection onto the fault plane');

sp = ceil(sqrt(num_usefulsensor));
figure1 = figure(f); f=f+1;
for i = 1:num_usefulsensor
    subplot(sp, sp, i);
    plot(cyc, cn(:, ro(i)));
%      ylabel('CN (N)');
%     xlabel('time');
end
han=axes(figure1,'visible','off');han.XLabel.Visible='on';han.YLabel.Visible='on';
ylabel(han,'CN (N)');xlabel(han,'PFC cycles');
figurename = strcat(fnumber, 'CN.png'); saveas(figure1,figurename);

figure1 = figure(f); f=f+1;
for i = 1:num_usefulsensor
    subplot(sp, sp, i);
    plot(cyc, cs(:, ro(i)));
%      ylabel('CS (N)');
%     xlabel('time');
end
han=axes(figure1,'visible','off');han.XLabel.Visible='on';han.YLabel.Visible='on';
ylabel(han,'CS (N)');xlabel(han,'PFC cycles');
figurename = strcat(fnumber, 'CS.png'); saveas(figure1,figurename);

figure1 = figure(f); f=f+1;
for i = 1:num_usefulsensor
    subplot(sp, sp, i);
    plot(cyc, bdz(:, ro(i)));
%      ylabel('bdz (m)');
%     xlabel('time');
end
han=axes(figure1,'visible','off');han.XLabel.Visible='on';han.YLabel.Visible='on';
ylabel(han,'bdz (m)');xlabel(han,'PFC cycles');
figurename = strcat(fnumber, 'bdz.png'); saveas(figure1,figurename);

figure1 = figure(f); f=f+1;
for i = 1:num_usefulsensor
    subplot(sp, sp, i);
    plot(cyc, bvz(:, ro(i)));
%     ylabel('bvz (m/s)');
%     xlabel('time');
end
han=axes(figure1,'visible','off');han.XLabel.Visible='on';han.YLabel.Visible='on';
ylabel(han,'bvz (m/s)');xlabel(han,'PFC cycles');
figurename = strcat(fnumber, 'bvz.png'); saveas(figure1,figurename);

figure1 = figure(f); f=f+1;
for i = 1:num_usefulsensor
    subplot(sp, sp, i);
    plot(cyc, bdy(:, ro(i)));
%      ylabel('bdy (m)');
%     xlabel('time');
end
han=axes(figure1,'visible','off');han.XLabel.Visible='on';han.YLabel.Visible='on';
ylabel(han,'bdy (m)');xlabel(han,'PFC cycles');
figurename = strcat(fnumber, 'bdy.png'); saveas(figure1,figurename);

figure1 = figure(f); f=f+1;
for i = 1:num_usefulsensor
    subplot(sp, sp, i);
    plot(cyc, bvy(:, ro(i)));
%     ylabel('bvy (m/s)');
%     xlabel('time');
end
han=axes(figure1,'visible','off');han.XLabel.Visible='on';han.YLabel.Visible='on';
ylabel(han,'bvy (m/s)');xlabel(han,'PFC cycles');
figurename = strcat(fnumber, 'bvy.png'); saveas(figure1,figurename);

figure1 = figure(f); f=f+1;
for i = 1:num_usefulsensor
    subplot(sp, sp, i);
    plot(cyc, bdx(:, ro(i)));
%      ylabel('bdx (m)');
%     xlabel('time');
end
han=axes(figure1,'visible','off');han.XLabel.Visible='on';han.YLabel.Visible='on';
ylabel(han,'bdx (m)');xlabel(han,'PFC cycles');
figurename = strcat(fnumber, 'bdx.png'); saveas(figure1,figurename);

figure1 = figure(f); f=f+1;
for i = 1:num_usefulsensor
    subplot(sp, sp, i);
    plot(cyc, bvx(:, ro(i)));
%     ylabel('bvx (m/s)');
%     xlabel('time');
end
han=axes(figure1,'visible','off');han.XLabel.Visible='on';han.YLabel.Visible='on';
ylabel(han,'bvx (m/s)');xlabel(han,'PFC cycles');
figurename = strcat(fnumber, 'bvx.png'); saveas(figure1,figurename);
%}
%
% Crack, AE, MS
if strcmp(fnumber, '201') == 1 
    [sjcrk, crk_cyc, is_sj_crk, sj_strength, is_shear] = openfile_crk2(strcat('crk',fnumber,'-5.txt'));
elseif strcmp(fnumber, '202') == 1 
    [sjcrk, crk_cyc, is_sj_crk, sj_strength, is_shear] = openfile_crk2(strcat('crk',fnumber,'-2.txt'));
else
    [sjcrk, crk_cyc, is_sj_crk, sj_strength, is_shear] = openfile_crk2(strcat('crk',fnumber,'-1.txt'));
end
temp = 1:length(sjcrk(:,1))-1;
d_sjcrk = diff(sjcrk(:,1));
%f = 1;
%{
figure(f);f=f+1;
plot(temp, d_sjcrk);
%}
%
%{
if strcmp(fnumber, '206')
    d_sjcrk_threshold = 6000;
elseif strcmp(fnumber, '210')
    d_sjcrk_threshold = 1.6e4;
elseif strcmp(fnumber, '201')
    d_sjcrk_threshold = 7900;
elseif strcmp(fnumber, '202')
    d_sjcrk_threshold = 7900;
elseif strcmp(fnumber, '207')
    d_sjcrk_threshold = 15000;
elseif strcmp(fnumber, '205')
    d_sjcrk_threshold = 8400;
elseif strcmp(fnumber, '204')
    d_sjcrk_threshold = 6000;
elseif strcmp(fnumber, '208')
    d_sjcrk_threshold = 6000;
elseif strcmp(fnumber, '214')
    d_sjcrk_threshold = 7500;
elseif strcmp(fnumber, '215')
    d_sjcrk_threshold = 800; %5000;
elseif strcmp(fnumber, '216')
    d_sjcrk_threshold = 1400; %5000;
elseif strcmp(fnumber, '220')
    d_sjcrk_threshold = 3000;
elseif strcmp(fnumber, '221')
    d_sjcrk_threshold = 2500;
elseif strcmp(fnumber, '222')
    d_sjcrk_threshold = 2000;
elseif strcmp(fnumber, '223')
    d_sjcrk_threshold = 350; %5000;
elseif strcmp(fnumber, '224')
    d_sjcrk_threshold = 1000;%6500; 
elseif strcmp(fnumber, '225')
    d_sjcrk_threshold = 1400; 
elseif strcmp(fnumber, '226')
    d_sjcrk_threshold = 1120; 
elseif strcmp(fnumber, '227')
    d_sjcrk_threshold = 3500; 
elseif strcmp(fnumber, '228')
    d_sjcrk_threshold = 3000; 
end

[pks, locs] = findpeaks(d_sjcrk, 'MinPeakHeight', d_sjcrk_threshold);
%}
%{
figure(f);f=f+1;
plot(temp, d_sjcrk);
hold on
scatter(temp(locs), d_sjcrk(locs));
hold off
%}
%
% AE
% if strcmp(fnumber, '201')
%     filename = fopen(strcat('ae',fnumber,'-5.txt'));
% elseif strcmp(fnumber, '202')
%     filename = fopen(strcat('ae',fnumber,'-2.txt'));
% else
%     filename = fopen(strcat('ae',fnumber,'-1.txt'));
% end
% %id   cyc0    x      y      z     M11  M12  M22  M31  M32  M33  cluster_id
% C = textscan(filename, '%f %f %f %f %f %f %f %f %f %f %f %f');
% ae = cell2mat(C);
% ae = sortrows(ae, 12);
% 
% clear temp
% for i = 1:length(sjcrk(:,1))
%     temp = find(ae(:,2) == sjcrk(i,1));
%     ae(temp, 14) = 1;
% end
% 
% [ae_unique, ia_ae, ic_ae] = unique(ae(:,12), 'rows', 'stable');
% clear ae_unique
% ae_unique = ae(ia_ae,:);
% clear temp
% for i = 1:size(ae_unique, 1)
%     %find all the ae whose 12th column is one of the ae_unique(:,12)
%     %find if within this group, if there is at least one 1 in the 18th
%     %column
%     %if there is at least one 1 in 18th column, calculate the average
%     %coordinate, calculate total moment-tensor, write that into a new array
%     temp = find(ae(:, 12) == ae_unique(i, 12));
%     ae(temp, 13) = i;
%     if nnz(ae(temp, 14)) > 0
%         ae(temp, 15) = mean(ae(temp,2), 1); %cycle
%         ae(temp, 16) = mean(ae(temp,3), 1); %average x location
%         ae(temp, 17) = mean(ae(temp,4), 1); %average x location
%         ae(temp, 18) = mean(ae(temp,5), 1); %average x location
%         ae(temp, 19) = sum(ae(temp,6));
%         ae(temp, 20) = sum(ae(temp,7));
%         ae(temp, 21) = sum(ae(temp,8));
%         ae(temp, 22) = sum(ae(temp,9));
%         ae(temp, 23) = sum(ae(temp,10));
%         ae(temp, 24) = sum(ae(temp,11));
%     end
% end
% 
% 
% temp = ae(:, 13:24);
% j = 1;
% for i = 1:size(temp,1)
%     if temp(i, 3) ~= 0
%         ae_sj(j, :) = temp(i, :);
%         j = j+1;
%     end
% end
% 
% [ae_sj_cluster, ia_aesj, ic_aesj] = unique(ae_sj(:,1), 'rows', 'stable');
% ae_sj_cluster = ae_sj(ia_aesj,:);

%
% MS

% filename = fopen(strcat('ms',fnumber,'-1.txt'));
% 
% %id cyc0 cyclast x y z     rad  mag M11  M12  M22  M31  M32  M33  iso  dev  aenum
% C = textscan(filename, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
% ms = cell2mat(C);
% 
% clear temp
% clear temp2
% for i = 1:length(sjcrk(:,1))
%     temp = find(ms(:,2) == sjcrk(i,1));
%     temp2 = find(ms(:,3) == sjcrk(i,1));
%     ms(temp, 18) = 1;
%     ms(temp2, 18) = 1;
% end
% ms_sj = ms(find(ms(:,18)==1),:);
% 
% % MS filtering
% ms_trans = [1 0 0; 0 1 0; sqrt(3) 0 1];
% ms_sj_temp = ms_trans * [ms_sj(:,4), ms_sj(:,5), ms_sj(:,6)]';
% sjcrk_transcoord = ms_trans * [sjcrk(:,2),sjcrk(:,3),sjcrk(:,4)]';
% ms_sj((abs(ms_sj_temp(3,:))>0.01),:) = [];
% % figure(f); f = f+1;
% % scatter3(sjcrk(:,2),sjcrk(:,3),sjcrk(:,4), '*');
% % hold on
% % scatter3(ms_sj(:,4), ms_sj(:,5), ms_sj(:,6));
% % % hold on
% % % scatter3(ms_sj_temp(1,:), ms_sj_temp(2,:), ms_sj_temp(3,:));
% % % hold on
% % % scatter3(sjcrk_transcoord(1,:), sjcrk_transcoord(2,:), sjcrk_transcoord(3,:), '*');
% % hold off
% % view([0,0]);


% % AE unclustered
% % unclustered AE, moment tensor from every single crack
% if strcmp(fnumber, '201') == 1
%     filename = fopen(strcat('ae-unclustered',fnumber,'-5.txt'));
% else
%     filename = fopen(strcat('ae-unclustered',fnumber,'-1.txt'));
% end
% %id   cyc0    x      y      z     M11  M12  M22  M31  M32  M33  cluster_id
% C = textscan(filename, '%f %f %f %f %f %f %f %f %f %f %f %f');
% ae_un = cell2mat(C);
% 
% clear temp
% clear temp2
% 
% ae_un_transcoord = ( ms_trans * [ae_un(:,3), ae_un(:,4), ae_un(:,5)]' )';
% %{
% figure(f); f= f+1;
% scatter3(ae_un_transcoord(:,1), ae_un_transcoord(:,2), ae_un_transcoord(:,3));
% hold on
% scatter3(sjcrk_transcoord(1,:), sjcrk_transcoord(2,:), sjcrk_transcoord(3,:) ,'k*' );
% hold off
% view([0,0]);
% %}
% % ae_un_sj = ae_un;
% % ae_un_sj((abs(ae_un_transcoord(3,:))>0.01),:) = [];
% 
% %for i = 1:length(ae_un(:,1))
%     %temp = find(ae_un(:,2) == sjcrk(i,1));
%     %temp = find(ae_un(:,5)<=0.03);
%     temp = find(abs(ae_un_transcoord(:,3))<0.01);
%     ae_un(temp, 12) = 1;
% %end
% ae_un_sj = ae_un(find(ae_un(:,12)==1),:);
% %{
% figure(f); f = f+1;
% scatter3(sjcrk(:,2),sjcrk(:,3),sjcrk(:,4), '*');
% hold on
% scatter3(ae_un_sj(:,3), ae_un_sj(:,4), ae_un_sj(:,5));
% hold off
% view([0,0]);
% %}
% % calculate seismic moment based on moment tensor
% M0_un = 1/sqrt(2) * sqrt( ae_un_sj(:,6).^2 + ae_un_sj(:,8).^2 + ae_un_sj(:,11).^2 ...
%     +  2*ae_un_sj(:,7).^2 +  2*ae_un_sj(:,9).^2 +  2*ae_un_sj(:,10).^2 ...
%     );
% Mw_un = 2/3*log10(M0_un) - 6.0;
% ae_un_sj(:,13) = M0_un;
% ae_un_sj(:,14) = Mw_un; % magnitudes of unclustered AE are stored in the 13 and 14th columns

%
% segment
%{
% locs obtained from the above command is based on sjcrk array, whose cycles
% are recorded in 1's instead of 20's as all other histories are recorded.
% Need to find the corresponding cycles recorded based on SH
temp = zeros(length(locs), 1);
locs2 = zeros(length(locs), 1);
for i = 1:length(locs)
    temp(i) = sjcrk(locs(i), 1);
    temp(i) = round(temp(i)/sample_rate) * sample_rate;
    locs2(i) = find(cyc == temp(i));
end
seg_bounds_ss = zeros(length(locs)+1, 1);
%seg_bound_offset = max(find(cyc<=sjcrk(1,1)));
for i = 1:length(locs)-1
    seg_bounds_ss(i+1) = ceil( (locs2(i) + locs2(i+1))/2 ); 
end
seg_bounds_ss(1) = 1;
seg_bounds_ss(end) = num_record;
sbs = seg_bounds_ss; %create a simpler variable just so that it can be called more easily
%}


% identify foreshock sequences
%is_sjcrk = zeros(size(is_sj_crk));
%is_sjcrk(nonzeros(is_sj_crk)) = 1;

sjcrk_cyc_clip = floor(sjcrk(:, 1) / 20 ) * 20;

cyc_sjclip = zeros(length(cyc),1);
for i = 1:length(sjcrk_cyc_clip)
    tempi = find(cyc == sjcrk_cyc_clip(i));
    %cyc_sjclip(tempi) = i;
    cyc_sjclip(tempi) = cyc(tempi);
end

pdSix = fitdist(sjcrk_cyc_clip,'Kernel','BandWidth',2000);
x_kernel = (1:1:length(cyc_sjclip))*sample_rate + cyc_offset;
ySix = pdf(pdSix,x_kernel); %ySix is the kernel density of the indices of sj cracks
[kernelmins, kernelminslocs] = findpeaks(-ySix);
[kernelpeaks, kernelpeakslocs] = findpeaks(ySix);
cyc_kernelpeakslocs = cyc(kernelpeakslocs);
%{
figure(f); f = f+1;
% plot(x_kernel, cyc_sjclip*1e-8, 'color',[0,0,0]+0.5);
% hold on
plot(x_kernel,ySix,'r-', 'linewidth', 1);
%hold off
ylim([0, 50e-5]);
%}
%{
figure(f); f = f+1;
plot(x_kernel,ySix,'r-', 'linewidth', 1);
hold on
scatter(x_kernel(kernelpeakslocs), kernelpeaks);
hold off
%}

seg_bounds_ss = zeros(length(kernelminslocs)+2, 1);
seg_bounds_ss(1) = 1;
seg_bounds_ss(end) = num_record;
seg_bounds_ss(2:end-1) = kernelminslocs';
sbs = seg_bounds_ss;
locs2 = kernelpeakslocs;

%segment sjcrk first to see if there is any empty rows. Combine the empty
%rows with the previus rows
seg_sjcrk = cell(length(sbs)-1, 1);
for i = 1:length(sbs)-1
    temp = find( (sjcrk(:,1))>cyc(sbs(i)) & (sjcrk(:,1))<cyc(sbs(i+1)) );
    temp_array = zeros(length(temp), 7);
    temp_array = sjcrk(temp,:);
    seg_sjcrk{i} = temp_array;
end

for i = 1:length(sbs)-1
    if isempty(seg_sjcrk{i}) == 1
        token = 1;
        break;
    else
        token = 0;
        continue;
    end
end

if token == 1
    clear temp
    clear temp2
    clear temp3 
    temp = zeros(length(sbs),1);
    temp(end) = sbs(end);
    temp2 = zeros(2, length(sbs)-1);
    for i = 1:length(sbs)-1
        if isempty(seg_sjcrk{i}) ~= 1
            temp(i) = sbs(i);
            temp2(1,i) = kernelpeaks(i);
            temp2(2,i) = kernelpeakslocs(i);

        else
            continue;
        end
    end
    temp(temp==0) = [];
    clear sbs
    sbs = temp;
    temp2(:, temp2(1,:)==0) = [];
    clear kernelpeaks
    kernelpeaks = temp2(1,:);
    clear kernelpeakslocs
    kernelpeakslocs = temp2(2,:);
    
    clear seg_sjcrk
    seg_sjcrk = cell(length(sbs)-1, 1);
    for i = 1:length(sbs)-1
        temp = find( (sjcrk(:,1))>cyc(sbs(i)) & (sjcrk(:,1))<cyc(sbs(i+1)) );
        temp_array = zeros(length(temp), 7);
        temp_array = sjcrk(temp,:);
        seg_sjcrk{i} = temp_array;
    end
    
    token = 0;
end

locs2 = kernelpeakslocs;
kernelpks_leverarm = zeros(length(kernelpeaks),1);
for i = 1:length(kernelpeaks)
    kernelpks_leverarm(i) = x_kernel(kernelpeakslocs(kernelpeaks == (max(kernelpeaks)))) - x_kernel(kernelpeakslocs(i)); 
end
%{
figure1 = figure(f); f= f+1;
plot(x_kernel, ySix, 'r-', 'linewidth', 1);
hold on
scatter(x_kernel(kernelpeakslocs), kernelpeaks, 'b*');
hold on
%scatter(x_kernel(kernelminslocs), kernelmins);
scatter(x_kernel(sbs), ySix(sbs));
hold off
xlabel('PFC cycles');
ylabel('Kernel density function');
figurename = strcat(fnumber, 'sj_kdf.png');
saveas(figure1,figurename);


figure1 = figure(f); f=f+1;
plot(t, shh);
hold on 
scatter(t(seg_bounds_ss), shh(seg_bounds_ss), 'go');
hold on 
scatter(t(locs2), shh(locs2), 'r*');
hold off;
xlabel('time');
ylabel('Stress (Pa)');
figurename = strcat(fnumber, 'seg_SH.png');
saveas(figure1, figurename);

figure1 = figure(f); f=f+1;
plot(ss(:,1), ss(:,2));
hold on 
scatter(ss(seg_bounds_ss,1), ss(seg_bounds_ss,2), 'go');
hold on 
scatter(ss(locs2,1), ss(locs2,2), 'r*');
hold off;
xlabel('stress');
ylabel('Stress (Pa)');
figurename = strcat(fnumber, 'seg_SS.png');
saveas(figure1, figurename);
%}
% Segment the contact forces based on seg_bounds_ss(sbs)
%close all
clear temp
clear temp_array

seg_cn_array = cell(length(sbs)-1, 1);
seg_t = cell(length(sbs)-1 , 1);
seg_cyc = cell(length(sbs)-1 , 1);
seg_strain = cell(length(sbs)-1, 1);
seg_stress = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = cn(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_cn_array{j, 1} = temp_array;
    seg_t{j, 1} = t(sbs(j): sbs(j+1));
    seg_cyc{j, 1} = cyc(sbs(j): sbs(j+1));
    seg_strain{j,1} = ss(sbs(j):sbs(j+1), 1);
    seg_stress{j,1} = ss(sbs(j):sbs(j+1), 2);
end
    
seg_cs_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = cs(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_cs_array{j, 1} = temp_array;
end
   
temp_t = seg_t{j};
temp_mean = mean(temp_array);
seg_mean =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean{j} = mean(  seg_cn_array{j}(:,:)  );
end
    
seg_mean_cs =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_cs{j} = mean(  seg_cs_array{j}(:,:)  );
end
    
seg_bdz_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bdz(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bdz_array{j, 1} = temp_array;
end

seg_bvz_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bvz(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bvz_array{j, 1} = temp_array;
end

seg_bdy_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bdy(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bdy_array{j, 1} = temp_array;
end

seg_bvy_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bvy(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bvy_array{j, 1} = temp_array;
end

seg_bdx_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bdx(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bdx_array{j, 1} = temp_array;
end

seg_bvx_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bvx(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bvx_array{j, 1} = temp_array;
end

seg_bdxoe_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bdxoe(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bdxoe_array{j, 1} = temp_array;
end

seg_bvxoe_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bvxoe(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bvxoe_array{j, 1} = temp_array;
end

seg_bdyoe_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bdyoe(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bdyoe_array{j, 1} = temp_array;
end

seg_bvyoe_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bvyoe(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bvyoe_array{j, 1} = temp_array;
end

seg_bdzoe_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bdzoe(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bdzoe_array{j, 1} = temp_array;
end

seg_bvzoe_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = bvzoe(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_bvzoe_array{j, 1} = temp_array;
end

seg_cdx_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = cdx(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_cdx_array{j, 1} = temp_array;
end

seg_cdy_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = cdy(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_cdy_array{j, 1} = temp_array;
end

seg_cdz_array = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        temp = cdz(sbs(j): sbs(j+1), i);
        temp_array(:, i ) = temp;
    end
    seg_cdz_array{j, 1} = temp_array;
end
    
seg_mean_bdz =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bdz{j} = mean(  seg_bdz_array{j}(:,:)  );
end

seg_mean_bvz =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bvz{j} = mean(  seg_bvz_array{j}(:,:)  );
end

seg_mean_bdy =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bdy{j} = mean(  seg_bdy_array{j}(:,:)  );
end

seg_mean_bvy =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bvy{j} = mean(  seg_bvy_array{j}(:,:)  );
end

seg_mean_bdx =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bdx{j} = mean(  seg_bdx_array{j}(:,:)  );
end

seg_mean_bvx =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bvx{j} = mean(  seg_bvx_array{j}(:,:)  );
end

seg_mean_bdxoe =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bdxoe{j} = mean(  seg_bdxoe_array{j}(:,:)  );
end

seg_mean_bvxoe =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bvxoe{j} = mean(  seg_bvxoe_array{j}(:,:)  );
end

seg_mean_bdyoe =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bdyoe{j} = mean(  seg_bdyoe_array{j}(:,:)  );
end

seg_mean_bvyoe =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bvyoe{j} = mean(  seg_bvyoe_array{j}(:,:)  );
end

seg_mean_bdzoe =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bdzoe{j} = mean(  seg_bdzoe_array{j}(:,:)  );
end

seg_mean_bvzoe =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_bvzoe{j} = mean(  seg_bvzoe_array{j}(:,:)  );
end

seg_mean_cdx =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_cdx{j} = mean(  seg_cdx_array{j}(:,:)  );
end

seg_mean_cdy =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_cdy{j} = mean(  seg_cdy_array{j}(:,:)  );
end

seg_mean_cdz =  cell(length(sbs)-1, 1);
for j = 1:length(sbs) - 1
    seg_mean_cdz{j} = mean(  seg_cdz_array{j}(:,:)  );
end

seg_cn_coef_array = cell(length(sbs)-1, 1);
seg_model = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    temp_coef_array = zeros(3, num_usefulsensor);
    temp_length = length([ sbs(j): sbs(j+1) ]);
    temp_model_array = zeros(temp_length, num_usefulsensor);
    for i = 1:num_usefulsensor
        [temp_coef_array(:, i), temp_model_array(:, i)] = onestep_linfit(seg_t{j}', seg_cn_array{j}(:, i));
    end
    seg_cn_coef_array{j} = temp_coef_array;
    seg_model{j} = temp_model_array;
end

    seg_cs_coef_array = cell(length(sbs)-1, 1);
    for j = 1:length(sbs)-1 % j window slides in time
        temp_coef_array = zeros(3, num_usefulsensor);
        for i = 1:num_usefulsensor
            [temp_coef, S] = polyfit(seg_t{j}', seg_cs_array{j, 1}(:, i), 1);
            temp_r2 = 1 - (S.normr/norm( seg_cs_array{j, 1}(:, i) - mean( seg_cs_array{j, 1}(:, i) )))^2;
            if isnan(temp_r2) == 1
                temp_r2 = 0; %some r2 are not NaN since the fitting data and the model coefficients are zero
            end
            temp_coef_array(:, i) = [temp_coef, temp_r2]; %The third number is r2
        end
        seg_cs_coef_array{j , 1} = temp_coef_array;
    end
    
    seg_model = cell(length(sbs)-1, 1);
    for j = 1:length(sbs)-1
        temp_length = length([ sbs(j): sbs(j+1) ]);
        temp_model_array = zeros(temp_length, num_usefulsensor);
        for i = 1:num_usefulsensor
            temp_model = (seg_cn_coef_array{j}(1, i)*seg_t{j} + seg_cn_coef_array{j}(2, i))';
            temp_model_array(:, i) = temp_model;
        end
        seg_model{j} = temp_model_array;
    end

    seg_cs_model = cell(length(sbs)-1, 1);
    for j = 1:length(sbs)-1
        temp_length = length([ sbs(j): sbs(j+1) ]);
        temp_model_array = zeros(temp_length, num_usefulsensor);
        for i = 1:num_usefulsensor
            temp_model = (seg_cs_coef_array{j}(1, i)*seg_t{j} + seg_cs_coef_array{j}(2, i))';
            temp_model_array(:, i) = temp_model;
        end
        seg_cs_model{j} = temp_model_array;
    end
  
    seg_cn_sd = cell(length(sbs)-1, 1);
    for j = 1:length(sbs)-1
        %temp_length = length([ sbs(j): sbs(j+1) ]);
        temp_model_array = zeros(1, num_usefulsensor);
        for i = 1:num_usefulsensor
            temp_sd = sum(...
                ( seg_cn_array{j}(:,i) - seg_cn_coef_array{j}(1,i)*seg_t{j}' - seg_cn_coef_array{j}(2,i) ).^2 ...
                );        
            temp_model_array(1, i) = temp_sd;
        end
        seg_cn_sd{j} = temp_model_array;
    end
    
    seg_cs_sd = cell(length(sbs)-1, 1);
    for j = 1:length(sbs)-1
        %temp_length = length([ sbs(j): sbs(j+1) ]);
        temp_model_array = zeros(1, num_usefulsensor);
        for i = 1:num_usefulsensor
            temp_sd = sum(...
                ( seg_cs_array{j}(:,i) - seg_cs_coef_array{j}(1,i)*seg_t{j}' - seg_cs_coef_array{j}(2,i) ).^2 ...
                );
            temp_model_array(1, i) = temp_sd;
        end
        seg_cs_sd{j} = temp_model_array;
    end
    
seg_is_shear = cell(length(sbs)-1, 1);
for i = 1:length(sbs)-1
    temp = find( (sjcrk(:,1))>cyc(sbs(i)) & (sjcrk(:,1))<cyc(sbs(i+1)) );
    temp_array = zeros(length(temp), 7);
    temp_array = is_shear(temp,:);
    seg_is_shear{i} = temp_array;
end

seg_sj_strength = cell(length(sbs)-1, 1);
for i = 1:length(sbs)-1
    temp = find( (sjcrk(:,1))>cyc(sbs(i)) & (sjcrk(:,1))<cyc(sbs(i+1)) );
    temp_array = zeros(length(temp), 7);
    temp_array = sj_strength(temp,:);
    seg_sj_strength{i} = temp_array;
end
%{
seg_ae_sj = cell(length(sbs)-1, 1);
for i = 1:length(sbs)-1
    temp = find( (ae_sj_cluster(:,3))>cyc(sbs(i)) & (ae_sj_cluster(:,3))<cyc(sbs(i+1)) );
    temp_array = zeros(length(temp), 12);
    temp_array = ae_sj_cluster(temp,:);
    seg_ae_sj{i} = temp_array;
end

seg_ms_sj = cell(length(sbs)-1, 1);
for i = 1:length(sbs)-1
    temp = find( (ms_sj(:,3))>cyc(sbs(i)) & (ms_sj(:,3))<cyc(sbs(i+1)) );
    temp_array = zeros(length(temp), 17);
    temp_array = ms_sj(temp,:);
    seg_ms_sj{i} = temp_array;
end

seg_ae_un_sj = cell(length(sbs)-1, 1);
for i = 1:length(sbs)-1
    temp = find( (ae_un_sj(:,2))>cyc(sbs(i)) & (ae_un_sj(:,2))<cyc(sbs(i+1)) );
    temp_array = zeros(length(temp), 14);
    temp_array = ae_un_sj(temp,:);
    seg_ae_un_sj{i} = temp_array;
end
%}
seg_sjstrength = cell(length(sbs)-1, 1);
for i = 1:length(sbs)-1
    temp = find( (sjcrk(:,2))>cyc(sbs(i)) & (sjcrk(:,2))<cyc(sbs(i+1)) );
    temp_array = zeros(length(temp), 1);
    temp_array = sj_strength(temp,:);
    seg_sjstrength{i} = temp_array;
end

ipt = zeros(1, length(sbs)-1);
ladder_idx_cs  = zeros(1, length(sbs)-1);
ladder_prejump_cs = cell(length(sbs)-1, 1);
ladder_postjump_cs = cell(length(sbs)-1, 1);
ladder_jump_cs = zeros( length(sbs)-1 , num_usefulsensor);
coef_prejump_cs = zeros(length(sbs)-1, num_usefulsensor, 3);
coef_postjump_cs = zeros(length(sbs)-1, num_usefulsensor, 3);
for j = 1:length(sbs)-1
    clear temp
    clear temp1
    clear temp2
    clear temp_length
    temp = mean(seg_cs_array{j},2);
    ladder_idx_cs(j) = find(seg_cyc{j}==cyc(locs2(j)));
    %ipt(j) = findchangepts(seg_cs_array{j}(:, 2));
    %ipt(j) = findchangepts(temp, 'Statistic','rms')

    for i = 1:num_usefulsensor        
%         temp = seg_cs_array{j}(:, i);      
%         ipt = findchangepts(temp);      
        temp_length = length(seg_cs_array{j}(1:ladder_idx_cs(j), i));
        temp1 = zeros(temp_length, 1);
        
        temp_length = length(seg_cs_array{j}(ladder_idx_cs(j):end, i));
        temp2 = zeros(temp_length, 1);
        
        [coef_prejump_cs(j,i,:),temp1] = onestep_linfit( seg_cyc{j}(1:ladder_idx_cs(j))', seg_cs_array{j}(1:ladder_idx_cs(j), i) );
        [coef_postjump_cs(j,i,:),temp2] = onestep_linfit(seg_cyc{j}(ladder_idx_cs(j):end)', seg_cs_array{j}(ladder_idx_cs(j):end, i));
        
        ladder_prejump_cs{j}(:, i) = temp1;
        ladder_postjump_cs{j}(:, i) = temp2;
        ladder_jump_cs(j,i) = temp2(1) - temp1(end); 
    end
end

% figure(f); f = f+1;
% sp1 = ceil(sqrt(num_usefulsensor)); sp2 = ceil(sqrt(num_usefulsensor));
% for j = 16:16%1:length(sbs)-1
%     for i = 1:num_usefulsensor
%         subplot(sp1, sp2, i);
%         plot(seg_cyc{j}, seg_cs_array{j}(:,i));
%         hold on
%         plot(seg_cyc{j}(1:ladder_idx_cs(j)) , ladder_prejump_cs{j}(:, i));
%         hold on
%         plot(seg_cyc{j}(ladder_idx_cs(j):end) , ladder_postjump_cs{j}(:, i));
%         hold off
%     end
% end

ladder_relative_cs = zeros(length(sbs)-1 , num_usefulsensor);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        ladder_relative_cs(j,i) = ladder_jump_cs(j, i)/seg_mean_cs{j}(i);
    end
end


% Multi segment linear fit to CN

ladder_idx_cn  = zeros(1, length(sbs)-1);
ladder_prejump_cn = cell(length(sbs)-1, 1);
ladder_postjump_cn = cell(length(sbs)-1, 1);
ladder_jump_cn = zeros( length(sbs)-1 , num_usefulsensor);
coef_prejump_cn = zeros(length(sbs)-1, num_usefulsensor, 3);
coef_postjump_cn = zeros(length(sbs)-1, num_usefulsensor, 3);
for j = 1:length(sbs)-1
    clear temp
    clear temp1
    clear temp2
    clear temp_length
    temp = mean(seg_cn_array{j},2);
    ladder_idx_cn(j) = find(seg_cyc{j}==cyc(locs2(j)));

    for i = 1:num_usefulsensor             
        temp_length = length(seg_cn_array{j}(1:ladder_idx_cn(j), i));
        temp1 = zeros(temp_length, 1);
        
        temp_length = length(seg_cn_array{j}(ladder_idx_cn(j):end, i));
        temp2 = zeros(temp_length, 1);
        
        [coef_prejump_cn(j,i,:),temp1] = onestep_linfit( seg_cyc{j}(1:ladder_idx_cn(j))', seg_cn_array{j}(1:ladder_idx_cn(j), i) );
        [coef_postjump_cn(j,i,:),temp2] = onestep_linfit(seg_cyc{j}(ladder_idx_cn(j):end)', seg_cn_array{j}(ladder_idx_cn(j):end, i));
        
        ladder_prejump_cn{j}(:, i) = temp1;
        ladder_postjump_cn{j}(:, i) = temp2;
        ladder_jump_cn(j,i) = temp2(1) - temp1(end); 
    end
end

% figure(f); f = f+1;
% sp1 = ceil(sqrt(num_usefulsensor)); sp2 = ceil(sqrt(num_usefulsensor));
% for j = 16:16%1:length(sbs)-1
%     for i = 1:num_usefulsensor
%         subplot(sp1, sp2, i);
%         plot(seg_cyc{j}, seg_cn_array{j}(:,i));
%         hold on
%         plot(seg_cyc{j}(1:ladder_idx_cn(j)) , ladder_prejump_cn{j}(:, i));
%         hold on
%         plot(seg_cyc{j}(ladder_idx_cn(j):end) , ladder_postjump_cn{j}(:, i));
%         hold off
%     end
% end
ladder_relative_cn = zeros(length(sbs)-1 , num_usefulsensor);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        ladder_relative_cn(j,i) = ladder_jump_cn(j, i)/seg_mean{j}(i);
    end
end


% Multi segment linear fit to bdz

ladder_idx_bdz  = zeros(1, length(sbs)-1);
ladder_prejump_bdz = cell(length(sbs)-1, 1);
ladder_postjump_bdz = cell(length(sbs)-1, 1);
ladder_jump_bdz = zeros( length(sbs)-1 , num_usefulsensor);
coef_prejump_bdz = zeros(length(sbs)-1, num_usefulsensor, 3);
coef_postjump_bdz = zeros(length(sbs)-1, num_usefulsensor, 3);
for j = 1:length(sbs)-1
    clear temp
    clear temp1
    clear temp2
    clear temp_length
    temp = mean(seg_bdz_array{j},2);
    ladder_idx_bdz(j) = find(seg_cyc{j}==cyc(locs2(j)));

    for i = 1:num_usefulsensor             
        temp_length = length(seg_bdz_array{j}(1:ladder_idx_bdz(j), i));
        temp1 = zeros(temp_length, 1);
        
        temp_length = length(seg_bdz_array{j}(ladder_idx_bdz(j):end, i));
        temp2 = zeros(temp_length, 1);
        
        [coef_prejump_bdz(j,i,:),temp1] = onestep_linfit( seg_cyc{j}(1:ladder_idx_bdz(j))', seg_bdz_array{j}(1:ladder_idx_bdz(j), i) );
        [coef_postjump_bdz(j,i,:),temp2] = onestep_linfit(seg_cyc{j}(ladder_idx_bdz(j):end)', seg_bdz_array{j}(ladder_idx_bdz(j):end, i));
        
        ladder_prejump_bdz{j}(:, i) = temp1;
        ladder_postjump_bdz{j}(:, i) = temp2;
        ladder_jump_bdz(j,i) = temp2(1) - temp1(end); 
    end
end

ladder_relative_bdz = zeros(length(sbs)-1 , num_usefulsensor);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        ladder_relative_bdz(j,i) = ladder_jump_bdz(j, i)/seg_mean_bdz{j}(i);
    end
end

% Multi segment linear fit to bdy

ladder_idx_bdy  = zeros(1, length(sbs)-1);
ladder_prejump_bdy = cell(length(sbs)-1, 1);
ladder_postjump_bdy = cell(length(sbs)-1, 1);
ladder_jump_bdy = zeros( length(sbs)-1 , num_usefulsensor);
coef_prejump_bdy = zeros(length(sbs)-1, num_usefulsensor, 3);
coef_postjump_bdy = zeros(length(sbs)-1, num_usefulsensor, 3);
for j = 1:length(sbs)-1
    clear temp
    clear temp1
    clear temp2
    clear temp_length
    temp = mean(seg_bdy_array{j},2);
    ladder_idx_bdy(j) = find(seg_cyc{j}==cyc(locs2(j)));

    for i = 1:num_usefulsensor             
        temp_length = length(seg_bdy_array{j}(1:ladder_idx_bdy(j), i));
        temp1 = zeros(temp_length, 1);
        
        temp_length = length(seg_bdy_array{j}(ladder_idx_bdy(j):end, i));
        temp2 = zeros(temp_length, 1);
        
        [coef_prejump_bdy(j,i,:),temp1] = onestep_linfit( seg_cyc{j}(1:ladder_idx_bdy(j))', seg_bdy_array{j}(1:ladder_idx_bdy(j), i) );
        [coef_postjump_bdy(j,i,:),temp2] = onestep_linfit(seg_cyc{j}(ladder_idx_bdy(j):end)', seg_bdy_array{j}(ladder_idx_bdy(j):end, i));
        
        ladder_prejump_bdy{j}(:, i) = temp1;
        ladder_postjump_bdy{j}(:, i) = temp2;
        ladder_jump_bdy(j,i) = temp2(1) - temp1(end); 
    end
end

% figure(f); f = f+1;
% sp1 = ceil(sqrt(num_usefulsensor)); sp2 = ceil(sqrt(num_usefulsensor));
% for j = 16:16%1:length(sbs)-1
%     for i = 1:num_usefulsensor
%         subplot(sp1, sp2, i);
%         plot(seg_cyc{j}, seg_bdy_array{j}(:,i));
%         hold on
%         plot(seg_cyc{j}(1:ladder_idx_bdy(j)) , ladder_prejump_bdy{j}(:, i));
%         hold on
%         plot(seg_cyc{j}(ladder_idx_bdy(j):end) , ladder_postjump_bdy{j}(:, i));
%         hold off
%     end
% end

ladder_relative_bdy = zeros(length(sbs)-1 , num_usefulsensor);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        ladder_relative_bdy(j,i) = ladder_jump_bdy(j, i)/seg_mean_bdy{j}(i);
    end
end

% Multi segment linear fit to bdx

ladder_idx_bdx  = zeros(1, length(sbs)-1);
ladder_prejump_bdx = cell(length(sbs)-1, 1);
ladder_postjump_bdx = cell(length(sbs)-1, 1);
ladder_jump_bdx = zeros( length(sbs)-1 , num_usefulsensor);
coef_prejump_bdx = zeros(length(sbs)-1, num_usefulsensor, 3);
coef_postjump_bdx = zeros(length(sbs)-1, num_usefulsensor, 3);
for j = 1:length(sbs)-1
    clear temp
    clear temp1
    clear temp2
    clear temp_length
    temp = mean(seg_bdx_array{j},2);
    ladder_idx_bdx(j) = find(seg_cyc{j}==cyc(locs2(j)));

    for i = 1:num_usefulsensor             
        temp_length = length(seg_bdx_array{j}(1:ladder_idx_bdx(j), i));
        temp1 = zeros(temp_length, 1);
        
        temp_length = length(seg_bdx_array{j}(ladder_idx_bdx(j):end, i));
        temp2 = zeros(temp_length, 1);
        
        [coef_prejump_bdx(j,i,:),temp1] = onestep_linfit( seg_cyc{j}(1:ladder_idx_bdx(j))', seg_bdx_array{j}(1:ladder_idx_bdx(j), i) );
        [coef_postjump_bdx(j,i,:),temp2] = onestep_linfit(seg_cyc{j}(ladder_idx_bdx(j):end)', seg_bdx_array{j}(ladder_idx_bdx(j):end, i));
        
        ladder_prejump_bdx{j}(:, i) = temp1;
        ladder_postjump_bdx{j}(:, i) = temp2;
        ladder_jump_bdx(j,i) = temp2(1) - temp1(end); 
    end
end

% figure(f); f = f+1;
% sp1 = ceil(sqrt(num_usefulsensor)); sp2 = ceil(sqrt(num_usefulsensor));
% for j = 16:16%1:length(sbs)-1
%     for i = 1:num_usefulsensor
%         subplot(sp1, sp2, i);
%         plot(seg_cyc{j}, seg_bdx_array{j}(:,i));
%         hold on
%         plot(seg_cyc{j}(1:ladder_idx_bdx(j)) , ladder_prejump_bdx{j}(:, i));
%         hold on
%         plot(seg_cyc{j}(ladder_idx_bdx(j):end) , ladder_postjump_bdx{j}(:, i));
%         hold off
%     end
% end

ladder_relative_bdx = zeros(length(sbs)-1 , num_usefulsensor);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        ladder_relative_bdx(j,i) = ladder_jump_bdx(j, i)/seg_mean_bdx{j}(i);
    end
end

% Multi segment linear fit to bdxoe

ladder_idx_bdxoe  = zeros(1, length(sbs)-1);
ladder_prejump_bdxoe = cell(length(sbs)-1, 1);
ladder_postjump_bdxoe = cell(length(sbs)-1, 1);
ladder_jump_bdxoe = zeros( length(sbs)-1 , num_usefulsensor);
coef_prejump_bdxoe = zeros(length(sbs)-1, num_usefulsensor, 3);
coef_postjump_bdxoe = zeros(length(sbs)-1, num_usefulsensor, 3);
for j = 1:length(sbs)-1
    clear temp
    clear temp1
    clear temp2
    clear temp_length
    temp = mean(seg_bdxoe_array{j},2);
    ladder_idx_bdxoe(j) = find(seg_cyc{j}==cyc(locs2(j)));

    for i = 1:num_usefulsensor             
        temp_length = length(seg_bdxoe_array{j}(1:ladder_idx_bdxoe(j), i));
        temp1 = zeros(temp_length, 1);
        
        temp_length = length(seg_bdxoe_array{j}(ladder_idx_bdxoe(j):end, i));
        temp2 = zeros(temp_length, 1);
        
        [coef_prejump_bdxoe(j,i,:),temp1] = onestep_linfit( seg_cyc{j}(1:ladder_idx_bdxoe(j))', seg_bdxoe_array{j}(1:ladder_idx_bdxoe(j), i) );
        [coef_postjump_bdxoe(j,i,:),temp2] = onestep_linfit(seg_cyc{j}(ladder_idx_bdxoe(j):end)', seg_bdxoe_array{j}(ladder_idx_bdxoe(j):end, i));
        
        ladder_prejump_bdxoe{j}(:, i) = temp1;
        ladder_postjump_bdxoe{j}(:, i) = temp2;
        ladder_jump_bdxoe(j,i) = temp2(1) - temp1(end); 
    end
end

ladder_relative_bdxoe = zeros(length(sbs)-1 , num_usefulsensor);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        ladder_relative_bdxoe(j,i) = ladder_jump_bdxoe(j, i)/seg_mean_bdxoe{j}(i);
    end
end

% Multi segment linear fit to bdyoe

ladder_idx_bdyoe  = zeros(1, length(sbs)-1);
ladder_prejump_bdyoe = cell(length(sbs)-1, 1);
ladder_postjump_bdyoe = cell(length(sbs)-1, 1);
ladder_jump_bdyoe = zeros( length(sbs)-1 , num_usefulsensor);
coef_prejump_bdyoe = zeros(length(sbs)-1, num_usefulsensor, 3);
coef_postjump_bdyoe = zeros(length(sbs)-1, num_usefulsensor, 3);
for j = 1:length(sbs)-1
    clear temp
    clear temp1
    clear temp2
    clear temp_length
    temp = mean(seg_bdyoe_array{j},2);
    ladder_idx_bdyoe(j) = find(seg_cyc{j}==cyc(locs2(j)));

    for i = 1:num_usefulsensor             
        temp_length = length(seg_bdyoe_array{j}(1:ladder_idx_bdyoe(j), i));
        temp1 = zeros(temp_length, 1);
        
        temp_length = length(seg_bdyoe_array{j}(ladder_idx_bdyoe(j):end, i));
        temp2 = zeros(temp_length, 1);
        
        [coef_prejump_bdyoe(j,i,:),temp1] = onestep_linfit( seg_cyc{j}(1:ladder_idx_bdyoe(j))', seg_bdyoe_array{j}(1:ladder_idx_bdyoe(j), i) );
        [coef_postjump_bdyoe(j,i,:),temp2] = onestep_linfit(seg_cyc{j}(ladder_idx_bdyoe(j):end)', seg_bdyoe_array{j}(ladder_idx_bdyoe(j):end, i));
        
        ladder_prejump_bdyoe{j}(:, i) = temp1;
        ladder_postjump_bdyoe{j}(:, i) = temp2;
        ladder_jump_bdyoe(j,i) = temp2(1) - temp1(end); 
    end
end

ladder_relative_bdyoe = zeros(length(sbs)-1 , num_usefulsensor);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        ladder_relative_bdyoe(j,i) = ladder_jump_bdyoe(j, i)/seg_mean_bdyoe{j}(i);
    end
end

% Multi segment linear fit to bdzoe

ladder_idx_bdzoe  = zeros(1, length(sbs)-1);
ladder_prejump_bdzoe = cell(length(sbs)-1, 1);
ladder_postjump_bdzoe = cell(length(sbs)-1, 1);
ladder_jump_bdzoe = zeros( length(sbs)-1 , num_usefulsensor);
coef_prejump_bdzoe = zeros(length(sbs)-1, num_usefulsensor, 3);
coef_postjump_bdzoe = zeros(length(sbs)-1, num_usefulsensor, 3);
for j = 1:length(sbs)-1
    clear temp
    clear temp1
    clear temp2
    clear temp_length
    temp = mean(seg_bdzoe_array{j},2);
    ladder_idx_bdzoe(j) = find(seg_cyc{j}==cyc(locs2(j)));

    for i = 1:num_usefulsensor             
        temp_length = length(seg_bdzoe_array{j}(1:ladder_idx_bdzoe(j), i));
        temp1 = zeros(temp_length, 1);
        
        temp_length = length(seg_bdzoe_array{j}(ladder_idx_bdzoe(j):end, i));
        temp2 = zeros(temp_length, 1);
        
        [coef_prejump_bdzoe(j,i,:),temp1] = onestep_linfit( seg_cyc{j}(1:ladder_idx_bdzoe(j))', seg_bdzoe_array{j}(1:ladder_idx_bdzoe(j), i) );
        [coef_postjump_bdzoe(j,i,:),temp2] = onestep_linfit(seg_cyc{j}(ladder_idx_bdzoe(j):end)', seg_bdzoe_array{j}(ladder_idx_bdzoe(j):end, i));
        
        ladder_prejump_bdzoe{j}(:, i) = temp1;
        ladder_postjump_bdzoe{j}(:, i) = temp2;
        ladder_jump_bdzoe(j,i) = temp2(1) - temp1(end); 
    end
end

ladder_relative_bdzoe = zeros(length(sbs)-1 , num_usefulsensor);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        ladder_relative_bdzoe(j,i) = ladder_jump_bdzoe(j, i)/seg_mean_bdzoe{j}(i);
    end
end

% Multi segment linear fit to cdx

ladder_idx_cdx  = zeros(1, length(sbs)-1);
ladder_prejump_cdx = cell(length(sbs)-1, 1);
ladder_postjump_cdx = cell(length(sbs)-1, 1);
ladder_jump_cdx = zeros( length(sbs)-1 , num_usefulsensor);
coef_prejump_cdx = zeros(length(sbs)-1, num_usefulsensor, 3);
coef_postjump_cdx = zeros(length(sbs)-1, num_usefulsensor, 3);
for j = 1:length(sbs)-1
    clear temp
    clear temp1
    clear temp2
    clear temp_length
    temp = mean(seg_cdx_array{j},2);
    ladder_idx_cdx(j) = find(seg_cyc{j}==cyc(locs2(j)));

    for i = 1:num_usefulsensor             
        temp_length = length(seg_cdx_array{j}(1:ladder_idx_cdx(j), i));
        temp1 = zeros(temp_length, 1);
        
        temp_length = length(seg_cdx_array{j}(ladder_idx_cdx(j):end, i));
        temp2 = zeros(temp_length, 1);
        
        [coef_prejump_cdx(j,i,:),temp1] = onestep_linfit( seg_cyc{j}(1:ladder_idx_cdx(j))', seg_cdx_array{j}(1:ladder_idx_cdx(j), i) );
        [coef_postjump_cdx(j,i,:),temp2] = onestep_linfit(seg_cyc{j}(ladder_idx_cdx(j):end)', seg_cdx_array{j}(ladder_idx_cdx(j):end, i));
        
        ladder_prejump_cdx{j}(:, i) = temp1;
        ladder_postjump_cdx{j}(:, i) = temp2;
        ladder_jump_cdx(j,i) = temp2(1) - temp1(end); 
    end
end

ladder_relative_cdx = zeros(length(sbs)-1 , num_usefulsensor);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        ladder_relative_cdx(j,i) = ladder_jump_cdx(j, i)/seg_mean_cdx{j}(i);
    end
end

% Multi segment linear fit to cdy

ladder_idx_cdy  = zeros(1, length(sbs)-1);
ladder_prejump_cdy = cell(length(sbs)-1, 1);
ladder_postjump_cdy = cell(length(sbs)-1, 1);
ladder_jump_cdy = zeros( length(sbs)-1 , num_usefulsensor);
coef_prejump_cdy = zeros(length(sbs)-1, num_usefulsensor, 3);
coef_postjump_cdy = zeros(length(sbs)-1, num_usefulsensor, 3);
for j = 1:length(sbs)-1
    clear temp
    clear temp1
    clear temp2
    clear temp_length
    temp = mean(seg_cdy_array{j},2);
    ladder_idx_cdy(j) = find(seg_cyc{j}==cyc(locs2(j)));

    for i = 1:num_usefulsensor             
        temp_length = length(seg_cdy_array{j}(1:ladder_idx_cdy(j), i));
        temp1 = zeros(temp_length, 1);
        
        temp_length = length(seg_cdy_array{j}(ladder_idx_cdy(j):end, i));
        temp2 = zeros(temp_length, 1);
        
        [coef_prejump_cdy(j,i,:),temp1] = onestep_linfit( seg_cyc{j}(1:ladder_idx_cdy(j))', seg_cdy_array{j}(1:ladder_idx_cdy(j), i) );
        [coef_postjump_cdy(j,i,:),temp2] = onestep_linfit(seg_cyc{j}(ladder_idx_cdy(j):end)', seg_cdy_array{j}(ladder_idx_cdy(j):end, i));
        
        ladder_prejump_cdy{j}(:, i) = temp1;
        ladder_postjump_cdy{j}(:, i) = temp2;
        ladder_jump_cdy(j,i) = temp2(1) - temp1(end); 
    end
end

ladder_relative_cdy = zeros(length(sbs)-1 , num_usefulsensor);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        ladder_relative_cdy(j,i) = ladder_jump_cdy(j, i)/seg_mean_cdy{j}(i);
    end
end

% Multi segment linear fit to cdz

ladder_idx_cdz  = zeros(1, length(sbs)-1);
ladder_prejump_cdz = cell(length(sbs)-1, 1);
ladder_postjump_cdz = cell(length(sbs)-1, 1);
ladder_jump_cdz = zeros( length(sbs)-1 , num_usefulsensor);
coef_prejump_cdz = zeros(length(sbs)-1, num_usefulsensor, 3);
coef_postjump_cdz = zeros(length(sbs)-1, num_usefulsensor, 3);
for j = 1:length(sbs)-1
    clear temp
    clear temp1
    clear temp2
    clear temp_length
    temp = mean(seg_cdz_array{j},2);
    ladder_idx_cdz(j) = find(seg_cyc{j}==cyc(locs2(j)));

    for i = 1:num_usefulsensor             
        temp_length = length(seg_cdz_array{j}(1:ladder_idx_cdz(j), i));
        temp1 = zeros(temp_length, 1);
        
        temp_length = length(seg_cdz_array{j}(ladder_idx_cdz(j):end, i));
        temp2 = zeros(temp_length, 1);
        
        [coef_prejump_cdz(j,i,:),temp1] = onestep_linfit( seg_cyc{j}(1:ladder_idx_cdz(j))', seg_cdz_array{j}(1:ladder_idx_cdz(j), i) );
        [coef_postjump_cdz(j,i,:),temp2] = onestep_linfit(seg_cyc{j}(ladder_idx_cdz(j):end)', seg_cdz_array{j}(ladder_idx_cdz(j):end, i));
        
        ladder_prejump_cdz{j}(:, i) = temp1;
        ladder_postjump_cdz{j}(:, i) = temp2;
        ladder_jump_cdz(j,i) = temp2(1) - temp1(end); 
    end
end

ladder_relative_cdz = zeros(length(sbs)-1 , num_usefulsensor);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        ladder_relative_cdz(j,i) = ladder_jump_cdz(j, i)/seg_mean_cdz{j}(i);
    end
end

% Find if a contact is broken due to stress transfer
% find all the dead contacts for each segment
%broken_threshold = 2; %if the jump in force is twice the size of force moment prior to broken, it is considered broken due to triggering
broken_threshold = 0.05:0.05:1;
broken_stressthreshold = 0:5e6:90e6;
%nb = length(broken_threshold);
nb = length(broken_stressthreshold);
broken_intensitythreshold = linspace(0,70e7,nb);
broken_contact = cell(length(sbs)-1, 1);
broken_preratio_cn = cell(length(sbs)-1, 1);
broken_preratio_cs = cell(length(sbs)-1, 1);
broken_sjstrength = cell(length(sbs)-1, 1);
broken_stressspike_cs = cell(length(sbs)-1, 1);
broken_trigger = cell(length(sbs)-1, nb);%indices of the broken contacts due to triggering
broken_stresstrigger = cell(length(sbs)-1, nb);%indices of the broken contacts due to triggering
broken_vx = cell(length(sbs)-1, 1); %velocity of contact (average of bv and bv_oe)
broken_vy = cell(length(sbs)-1, 1);
broken_vz = cell(length(sbs)-1, 1);
broken_intensityx = cell(length(sbs)-1, 1);
broken_intensityy = cell(length(sbs)-1, 1);
broken_intensityz = cell(length(sbs)-1, 1);
broken_inttriggerx = cell(length(sbs)-1, 1);
broken_inttriggery = cell(length(sbs)-1, 1);
broken_inttriggerz = cell(length(sbs)-1, 1);
seg_cn_prejump = zeros(num_usefulsensor,length(sbs)-1);%cn forces prior to broken
seg_cs_prejump = zeros(num_usefulsensor,length(sbs)-1);%cs forces prior to broken

for i = 1:length(sbs)-1 %number of episodes
    broken_contact{i} = find(seg_cn_array{i}(1,:)~=0 ...
        & seg_cn_array{i}(end,:)==0); %indices of all the broken contacts (among those sampled) within i-th episode. It's a row vector
    tempi = broken_contact{i}(1,:); %indices of the broken contacts. temporarily store them in a separate variable for easier coding
    broken_preratio_cn{i} = (max(seg_cn_array{i}(:,tempi)) - ladder_prejump_cn{i}(end, tempi))...
        ./ladder_prejump_cn{i}(end, tempi);
    broken_preratio_cs{i} = (max(seg_cs_array{i}(:,tempi)) - ladder_prejump_cs{i}(end, tempi))...
        ./ladder_prejump_cs{i}(end, tempi);
    % reverse calculating stress pre-spike. First estimate roughly the failure strength
    % average area of contact is 2e-5; cohesion + tan(phi) * (average Fn_spike / average contact area) 
    broken_sjstrength{i} = ...
        (mean(max(seg_cn_array{i}(:,tempi)))/(2e-5))*tand(30) + 55e6; 
    broken_stressspike_cs{i} = 1./(1+1./broken_preratio_cs{i}) .* broken_sjstrength{i};
%     for j = 1:nb
%         temp1 = broken_contact{i}(broken_preratio_cs{i}>=broken_threshold(j));
%         temp2 = broken_contact{i}(broken_preratio_cs{i}>=broken_threshold(j));
%         temp3 = unique([temp1, temp2]);
%         broken_trigger{i,j} = temp3;
%     end
   
    % Intensity = pressure * particle velocity
    broken_vx{i} = mean( [max(seg_bvxoe_array{i}(:,tempi))...
        ;max(seg_bvx_array{i}(:,tempi))] );
    broken_vy{i} = mean( [max(seg_bvyoe_array{i}(:,tempi))...
        ;max(seg_bvy_array{i}(:,tempi))] );
    broken_vz{i} = mean( [max(seg_bvzoe_array{i}(:,tempi))...
        ;max(seg_bvz_array{i}(:,tempi))] );
    broken_intensityx{i} = broken_stressspike_cs{i} .* broken_vx{i};
    broken_intensityy{i} = broken_stressspike_cs{i} .* broken_vy{i};
    broken_intensityz{i} = broken_stressspike_cs{i} .* broken_vz{i};
    for j = 1:nb
        temp1 = broken_contact{i}( ...
            broken_intensityx{i}>=broken_intensitythreshold(j) );
        temp2 = broken_contact{i}( ...
            broken_intensityy{i}>=broken_intensitythreshold(j) );
        temp3 = broken_contact{i}( ...
            broken_intensityz{i}>=broken_intensitythreshold(j) );
        broken_inttriggerx{i,j} = temp1;
        broken_inttriggery{i,j} = temp2;
        broken_inttriggerz{i,j} = temp3;
        temp0 = broken_contact{i}( ...
            broken_stressspike_cs{i}>=broken_stressthreshold(j) );
        broken_stresstrigger{i,j} = temp0;
    end
    
    seg_cn_prejump(:,i) = ladder_prejump_cn{i}(end,:)';
    seg_cs_prejump(:,i) = ladder_prejump_cs{i}(end,:)';
end
save(strcat('S',fnumber,'broken_inttriggerx'),'broken_inttriggerx');
save(strcat('S',fnumber,'broken_inttriggery'),'broken_inttriggery');
save(strcat('S',fnumber,'broken_inttriggerz'),'broken_inttriggerz');
save(strcat('S',fnumber,'broken_stresstrigger'),'broken_stresstrigger');
fclose all
end