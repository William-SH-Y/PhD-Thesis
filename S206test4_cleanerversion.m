
% S67test4, code based on S67test1, refined based on S101test2
% For stacking data: S201, S208, S214, S215, S216, S220, S221, S222, S223,
% S224, S225, S226, S227, S228, S229, S230, S231, S232, S233, S234, S235,
% S236, S237, s238, s239
close all
clear
clc
farray = {'201','208','214','215','216','220','221','222','223','224','225', '226','227','228',...
    '229','230','231','232','233','234','235', '236','237','238','239','240','241','242', '243', '244',...
    '245','246','247','248','249','250','251','254','255',...
    '449','451','452','453','454','455','456','457','458','459','460',...
    '461','462','463','464','465','466','467','468','469','470','471','472','473','474','487',...
    '488',...
    '489',...
    '490',...
    '486',...
    '485',...
    '484',...
    '476',...
    '475',...
    '474',...
    '473',...
    '472',...
    '464',...
    '483',...
    '471',...
    '449',...
    '470',...
    '455',...
    '469',...
    '458',...
    '468',...
    '467',...
    '479',...
    '482',...
    '481',...
    '480',...
    '463'};
n_simulation = length(farray);
sim_workingonthis = '487';
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
elseif strcmp(sim_workingonthis(1),'4') == 1
    num_sensor = 100;
else
    num_sensor = 30;
end
num_rec = 6e4;
num_record = length(shh); %technically it is the total nonzero elements of bv, but here for simplicity sake (since bv is not read), we temporarily set it to total number of recordings per trial
%
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

read_step = textscan(fopen(strcat('global_step',fnumber,'-1.txt')), '%f');
if n_workingonthis  >= 40
    cyc_offset = read_step{1}(1)-1;
else
    cyc_offset = cyc_offset_array(n_workingonthis); 
    cyc_offset_array(n_workingonthis) = cyc_offset;
end

%save('cyc_offset_array','cyc_offset_array');

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
    
    num_usefulsensor = num_sensor - length(isnt_sj_ind)
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
        bvzoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvzoe_temp;
        bdzoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdzoe_temp;
        bvyoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvyoe_temp;
        bdyoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdyoe_temp;
        bvxoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvxoe_temp;
        bdxoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdxoe_temp;
    elseif strcmp(sim_workingonthis,'490')==1
        cn(:, ( popcnind: popcnind+num_usefulsensor-1) ) = cn_temp;
        cs(:, ( popcnind: popcnind+num_usefulsensor-1) ) = cs_temp;
        bvz(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bvz_temp;
        bdz(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bdz_temp;
        bvy(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bvy_temp;
        bdy(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bdy_temp;
        bvx(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bvx_temp;
        bdx(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bdx_temp;
        bvzoe(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bvzoe_temp;
        bdzoe(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bdzoe_temp;
        bvyoe(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bvyoe_temp;
        bdyoe(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bdyoe_temp;
        bvxoe(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bvxoe_temp;
        bdxoe(2:end, ( popcnind: popcnind+num_usefulsensor-1) ) = bdxoe_temp;
    else
    cn(:, ( popcnind: popcnind+num_usefulsensor-1) ) = cn_temp;
    cs(:, ( popcnind: popcnind+num_usefulsensor-1) ) = cs_temp;
    bvz(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvz_temp;
    bdz(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdz_temp;
    bvy(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvy_temp;
    bdy(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdy_temp;
    bvx(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvx_temp;
    bdx(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdx_temp;
    bvzoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvzoe_temp;
    bdzoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdzoe_temp;
    bvyoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvyoe_temp;
    bdyoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdyoe_temp;
    bvxoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bvxoe_temp;
    bdxoe(:, ( popcnind: popcnind+num_usefulsensor-1) ) = bdxoe_temp;
    end
    
    
    bcoord(( popcnind: popcnind+num_usefulsensor-1), :) = bcoord_temp;
    popcnind = popcnind+num_usefulsensor;
    %cn(:, (num_sensor*(k-2)+1):(num_sensor*(k-2)+1)+num_usefulsensor-1 )= cn_temp;
    %cs(:, (num_sensor*(k-2)+1):(num_sensor*(k-2)+1)+num_usefulsensor-1 )= cs_temp;
    
    
end
%
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
if strcmp(fnumber, '201')
    filename = fopen(strcat('ae',fnumber,'-5.txt'));
elseif strcmp(fnumber, '202')
    filename = fopen(strcat('ae',fnumber,'-2.txt'));
else
    filename = fopen(strcat('ae',fnumber,'-1.txt'));
end
%id   cyc0    x      y      z     M11  M12  M22  M31  M32  M33  cluster_id
C = textscan(filename, '%f %f %f %f %f %f %f %f %f %f %f %f');
ae = cell2mat(C);
ae = sortrows(ae, 12);

clear temp
for i = 1:length(sjcrk(:,1))
    temp = find(ae(:,2) == sjcrk(i,1));
    ae(temp, 14) = 1;
end

[ae_unique, ia_ae, ic_ae] = unique(ae(:,12), 'rows', 'stable');
clear ae_unique
ae_unique = ae(ia_ae,:);
clear temp
for i = 1:size(ae_unique, 1)
    %find all the ae whose 12th column is one of the ae_unique(:,12)
    %find if within this group, if there is at least one 1 in the 18th
    %column
    %if there is at least one 1 in 18th column, calculate the average
    %coordinate, calculate total moment-tensor, write that into a new array
    temp = find(ae(:, 12) == ae_unique(i, 12));
    ae(temp, 13) = i;
    if nnz(ae(temp, 14)) > 0
        ae(temp, 15) = mean(ae(temp,2), 1); %cycle
        ae(temp, 16) = mean(ae(temp,3), 1); %average x location
        ae(temp, 17) = mean(ae(temp,4), 1); %average x location
        ae(temp, 18) = mean(ae(temp,5), 1); %average x location
        ae(temp, 19) = sum(ae(temp,6));
        ae(temp, 20) = sum(ae(temp,7));
        ae(temp, 21) = sum(ae(temp,8));
        ae(temp, 22) = sum(ae(temp,9));
        ae(temp, 23) = sum(ae(temp,10));
        ae(temp, 24) = sum(ae(temp,11));
    end
end


temp = ae(:, 13:24);
j = 1;
for i = 1:size(temp,1)
    if temp(i, 3) ~= 0
        ae_sj(j, :) = temp(i, :);
        j = j+1;
    end
end

[ae_sj_cluster, ia_aesj, ic_aesj] = unique(ae_sj(:,1), 'rows', 'stable');
ae_sj_cluster = ae_sj(ia_aesj,:);

%
% MS

filename = fopen(strcat('ms',fnumber,'-1.txt'));

%id cyc0 cyclast x y z     rad  mag M11  M12  M22  M31  M32  M33  iso  dev  aenum
C = textscan(filename, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
ms = cell2mat(C);

clear temp
clear temp2
for i = 1:length(sjcrk(:,1))
    temp = find(ms(:,2) == sjcrk(i,1));
    temp2 = find(ms(:,3) == sjcrk(i,1));
    ms(temp, 18) = 1;
    ms(temp2, 18) = 1;
end
ms_sj = ms(find(ms(:,18)==1),:);

% MS filtering
ms_trans = [1 0 0; 0 1 0; sqrt(3) 0 1];
ms_sj_temp = ms_trans * [ms_sj(:,4), ms_sj(:,5), ms_sj(:,6)]';
sjcrk_transcoord = ms_trans * [sjcrk(:,2),sjcrk(:,3),sjcrk(:,4)]';
ms_sj((abs(ms_sj_temp(3,:))>0.01),:) = [];
% figure(f); f = f+1;
% scatter3(sjcrk(:,2),sjcrk(:,3),sjcrk(:,4), '*');
% hold on
% scatter3(ms_sj(:,4), ms_sj(:,5), ms_sj(:,6));
% % hold on
% % scatter3(ms_sj_temp(1,:), ms_sj_temp(2,:), ms_sj_temp(3,:));
% % hold on
% % scatter3(sjcrk_transcoord(1,:), sjcrk_transcoord(2,:), sjcrk_transcoord(3,:), '*');
% hold off
% view([0,0]);


% AE unclustered
% unclustered AE, moment tensor from every single crack
if strcmp(fnumber, '201') == 1
    filename = fopen(strcat('ae-unclustered',fnumber,'-5.txt'));
else
    filename = fopen(strcat('ae-unclustered',fnumber,'-1.txt'));
end
%id   cyc0    x      y      z     M11  M12  M22  M31  M32  M33  cluster_id
C = textscan(filename, '%f %f %f %f %f %f %f %f %f %f %f %f');
ae_un = cell2mat(C);

clear temp
clear temp2

ae_un_transcoord = ( ms_trans * [ae_un(:,3), ae_un(:,4), ae_un(:,5)]' )';
%{
figure(f); f= f+1;
scatter3(ae_un_transcoord(:,1), ae_un_transcoord(:,2), ae_un_transcoord(:,3));
hold on
scatter3(sjcrk_transcoord(1,:), sjcrk_transcoord(2,:), sjcrk_transcoord(3,:) ,'k*' );
hold off
view([0,0]);
%}
% ae_un_sj = ae_un;
% ae_un_sj((abs(ae_un_transcoord(3,:))>0.01),:) = [];

%for i = 1:length(ae_un(:,1))
    %temp = find(ae_un(:,2) == sjcrk(i,1));
    %temp = find(ae_un(:,5)<=0.03);
    temp = find(abs(ae_un_transcoord(:,3))<0.01);
    ae_un(temp, 12) = 1;
%end
ae_un_sj = ae_un(find(ae_un(:,12)==1),:);
%{
figure(f); f = f+1;
scatter3(sjcrk(:,2),sjcrk(:,3),sjcrk(:,4), '*');
hold on
scatter3(ae_un_sj(:,3), ae_un_sj(:,4), ae_un_sj(:,5));
hold off
view([0,0]);
%}
% calculate seismic moment based on moment tensor
M0_un = 1/sqrt(2) * sqrt( ae_un_sj(:,6).^2 + ae_un_sj(:,8).^2 + ae_un_sj(:,11).^2 ...
    +  2*ae_un_sj(:,7).^2 +  2*ae_un_sj(:,9).^2 +  2*ae_un_sj(:,10).^2 ...
    );
Mw_un = 2/3*log10(M0_un) - 6.0;
ae_un_sj(:,13) = M0_un;
ae_un_sj(:,14) = Mw_un; % magnitudes of unclustered AE are stored in the 13 and 14th columns

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

figure(f); f = f+1;
% plot(x_kernel, cyc_sjclip*1e-8, 'color',[0,0,0]+0.5);
% hold on
plot(x_kernel,ySix,'r-', 'linewidth', 1);
%hold off
ylim([0, 50e-5]);


figure(f); f = f+1;
plot(x_kernel,ySix,'r-', 'linewidth', 1);
hold on
scatter(x_kernel(kernelpeakslocs), kernelpeaks);
hold off
%}
%

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
%
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

seg_sjstrength = cell(length(sbs)-1, 1);
for i = 1:length(sbs)-1
    temp = find( (sjcrk(:,2))>cyc(sbs(i)) & (sjcrk(:,2))<cyc(sbs(i+1)) );
    temp_array = zeros(length(temp), 1);
    temp_array = sj_strength(temp,:);
    seg_sjstrength{i} = temp_array;
end
save(strcat('S',fnumber,'seg_sjstrength'),'seg_sjstrength');

% figure(f); f=f+1;
% for i = 1:num_usefulsensor
%     subplot(5, 8, i);
%     plot( seg_t{8}, seg_cs_array{8}(:, i));
%     ylabel('CN (N)');
%     xlabel('time');
% end
% 
% figure(f); f=f+1;
% for i = 1:num_usefulsensor
%     subplot(5, 8, i);
%     plot( seg_t{8}, seg_cs_array{8}(:, i));
%     hold on
%     plot( seg_t{8}, seg_cs_model{8}(:, i), 'k');
%     hold off
%     ylabel('CN (N)');
%     xlabel('time');
% end

        
% SJ propagation
sp1 = ceil(sqrt(length(sbs)));
sp2 = ceil(sqrt(length(sbs)));
%{
figure(f); f=f+1;
for i = 1:length(sbs)-1
    ah = subplot(sp1, sp2, i);
    plot(seg_sjcrk{i}(:,1), seg_sjcrk{i}(:,2));
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',...
    sprintf('Crack occurance (PFC cycle) vs. crack x coordinates'), 'EdgeColor', 'none');
figure(f); f=f+1;
for i = 1:length(sbs)-1
    ah = subplot(sp1, sp2, i);
    plot(seg_sjcrk{i}(:,1), seg_sjcrk{i}(:,3));
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',...
    sprintf('Crack occurance (PFC cycle) vs. crack y coordinates'), 'EdgeColor', 'none');
figure(f); f=f+1;
for i = 1:length(sbs)-1
    ah = subplot(sp1, sp2, i);
    plot(seg_sjcrk{i}(:,1), seg_sjcrk{i}(:,4));
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',...
    sprintf('Crack occurance (PFC cycle) vs. crack z coordinates'), 'EdgeColor', 'none');

figure(f); f=f+1;
for i = 1:length(sbs)-1
    ah = subplot(sp1, sp2, i);
    plot(seg_sjcrk{i}(:,1), 1:1:length(seg_sjcrk{i}(:,1)));
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',...
    sprintf('Accumulative crack counts'), 'EdgeColor', 'none');
%}



%measurement of how offcentered the cracks are. Add up all the vectors
%(coordinates of the cracks), and then normalize it. If the cracks are
%scattered in all direction with respect to zero, it would add up to close
%to 0, but if there are located on one side, it would add up to be far from
%zero, depends on how clustered they are.
sj_spatdist = zeros(length(sbs)-1, 3); %smoothjoint cracks spatial distriction vector
for i = 1:length(sbs)-1
    %sj_I(i,1) = norm
    sj_spatdist(i,1) = mean(seg_sjcrk{i}(:,2));
    sj_spatdist(i,2) = mean(seg_sjcrk{i}(:,3));
    sj_spatdist(i,3) = mean(seg_sjcrk{i}(:,4));
    %sj_spatdist(i,:) = sj_spatdist(i,:)/norm(sj_spatdist(i,:));
end
%{
figure1 = figure(f);f=f+1;
scatter3(sj_spatdist(:,1), sj_spatdist(:,2), sj_spatdist(:,3));
view([0, 90]);
title('Cracks (in each time window) spacial distribution vectors');
figurename = strcat(fnumber, 'sjcrk_spacial.png');
saveas(figure1, figurename);

sp1 = ceil(sqrt(length(sbs)));
sp2 = ceil(sqrt(length(sbs)));
figure1 = figure(f); f=f+1;
for i = 1:length(sbs)-1
    ah = subplot(sp1, sp2, i);
    a = [1:length(seg_sjcrk{i}(:,2))]'; b = num2str(a); c=cellstr(b);
    scatter3(seg_sjcrk{i}(:,2), seg_sjcrk{i}(:,3), seg_sjcrk{i}(:,4));
    %text(seg_sjcrk{i}(:,2), seg_sjcrk{i}(:,3), seg_sjcrk{i}(:,4), c);
    view([0, 90]); %x-z plane
    %view([0, 0]); %x-y plane
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]); 
    clear a
    clear b
    clear c
    %plot(seg_sjcrk{i}(:,1), seg_sjcrk{i}(:,2));
    hold on
    scatter3(sj_spatdist(i,1), sj_spatdist(i,2), sj_spatdist(i,3),'r');
    hold off
end
figurename = strcat(fnumber, 'sjprop_timewindows.png');
saveas(figure1,figurename);
%}


%Focus on stick slip, finding fracture propagation velocity 
%Identify which time window contains stick slip
temp = zeros(length(sbs)-1, 1);
for i = 1:length(sbs)-1
    temp(i) = size(seg_sjcrk{i},1);
end
window_stickslip = find(temp == max(temp));
%

clear temp
temp = size(seg_sjcrk{window_stickslip}, 1);
temp_n = 25;
ub = floor(temp/temp_n).*(1:1:temp_n);
lb = floor(temp/temp_n).*(1:1:temp_n) - floor(temp/temp_n) + 1;
stickslip_occur = zeros(floor(temp/temp_n), temp_n);
stickslip_x = zeros(floor(temp/temp_n), temp_n);
stickslip_y = zeros(floor(temp/temp_n), temp_n);
stickslip_z = zeros(floor(temp/temp_n), temp_n);
for i = 1:temp_n
    stickslip_occur(:, i) = seg_sjcrk{window_stickslip}(lb(i):ub(i), 1);
    stickslip_x(:, i) = seg_sjcrk{window_stickslip}(lb(i):ub(i), 2);
    stickslip_y(:, i) = seg_sjcrk{window_stickslip}(lb(i):ub(i), 3);
    stickslip_z(:, i) = seg_sjcrk{window_stickslip}(lb(i):ub(i), 4);
end
%{
figure(f); f=f+1;
sp1 = ceil(sqrt(temp_n));
sp2 = ceil(sqrt(temp_n));
for i = 1:temp_n
    ah = subplot(sp1, sp2, i);
    %scatter3(stickslip_x(:,i), stickslip_y(:,i), stickslip_z(:,i));
    scatter3(reshape(stickslip_x(:,1:i), [], 1), reshape(stickslip_y(:,1:i), [], 1), reshape(stickslip_z(:,1:i), [], 1));
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);
    %plot(stickslip_occur(:,i), stickslip_z(:,i));
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',sprintf('Crack propagation during stick slip, splitted into 25 sub-episodes, projected onto x-y plane.'), 'EdgeColor', 'none');

figure(f); f=f+1;
sp1 = ceil(sqrt(temp_n));
sp2 = ceil(sqrt(temp_n));
for i = 1:temp_n
    ah = subplot(sp1, sp2, i);
    scatter3(reshape(stickslip_x(:,1:i), [], 1), reshape(stickslip_y(:,1:i), [], 1), reshape(stickslip_z(:,1:i), [], 1));
    view([0, 0]);
    xlim([-20e-3, 20e-3]);
    zlim([-55e-3, 55e-3]);
    %plot(stickslip_occur(:,i), stickslip_z(:,i));
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',sprintf('Crack propagation during stick slip, splitted into 25 sub-episodes, projected onto x-z plane.'), 'EdgeColor', 'none');
%}
%
%{
figure(f); f = f+1;
clear temp
temp = size(seg_sjcrk{window_stickslip}, 1);
temp_n = 100;
figure(f); f = f+1;
for i = 1:temp_n
    subplot(sqrt(temp_n), sqrt(temp_n), i);
    scatter3(seg_sjcrk{window_stickslip}(1:i, 2), seg_sjcrk{window_stickslip}(1:i, 3), seg_sjcrk{window_stickslip}(1:i, 4));
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);
end

if strcmp(fnumber, '216') == 1
    n_ini = 4;
    n_inix = [seg_sjcrk{window_stickslip}(1,2) ; seg_sjcrk{window_stickslip}(8,2) ; seg_sjcrk{window_stickslip}(20,2) ; seg_sjcrk{window_stickslip}(49,2) ];
    n_iniy = [seg_sjcrk{window_stickslip}(1,3) ; seg_sjcrk{window_stickslip}(8,3) ; seg_sjcrk{window_stickslip}(20,3) ; seg_sjcrk{window_stickslip}(49,3) ];
    n_iniz = [seg_sjcrk{window_stickslip}(1,4) ; seg_sjcrk{window_stickslip}(8,4) ; seg_sjcrk{window_stickslip}(20,4) ; seg_sjcrk{window_stickslip}(49,4) ];
end
%}

%figure(f); f=f+1;
%plot(seg_strain{window_stickslip}, seg_stress{window_stickslip});
clear temp1
clear temp2

ind_stress_ws_min = find(seg_stress{window_stickslip} == min(seg_stress{window_stickslip}) );
ind_stress_ws_max = find(seg_stress{window_stickslip} == max(seg_stress{window_stickslip}(1:ind_stress_ws_min)) );

% temp1 = seg_cyc{window_stickslip}( seg_stress{window_stickslip} == max(seg_stress{window_stickslip}) );
% temp2 = seg_cyc{window_stickslip}( seg_stress{window_stickslip} == min(seg_stress{window_stickslip}) );
temp1 = seg_cyc{window_stickslip}(ind_stress_ws_max);
temp2 = seg_cyc{window_stickslip}(ind_stress_ws_min);

sjcrk_ss = seg_sjcrk{window_stickslip}( (seg_sjcrk{window_stickslip}(:,1)>=temp1&seg_sjcrk{window_stickslip}(:,1)<=temp2) , : );
%figure(f); f = f+1;
%plot(sjcrk_ss(:,1), sjcrk_ss(:,2));


cyc_stress_ws_max = temp1;
cyc_stress_ws_min = temp2;


% apply a low pass filter to sj propagation during stickslip to find its propagation pass
signalt = sjcrk_ss(:,1);
signalx = sjcrk_ss(:,2);
signaly = sjcrk_ss(:,3);
fs = sample_rate;
%{
figure(f); f = f+1;
temp = lowpass(signaly,0.0001,fs);
figure(f); f = f+1;
plot(signalt, temp);
hold on
plot(signalt,signaly);
hold off
%}
% low pass filter doesn't work well. Try polynomial fitting
%p_ss_sj = polyfit(signalt, signalx, 4);
p_order = 4; %polynomial order
p_ss_sj = zeros(3, p_order+1);
pf_ss_sj = zeros(length(sjcrk_ss(:,1)), 3);
for i = 2:4
    p_ss_sj(i-1,:) = polyfit(sjcrk_ss(:,1), sjcrk_ss(:,i), p_order);
    pf_ss_sj(:, i-1) = polyval(p_ss_sj(i-1,:), sjcrk_ss(:,1));
end
%signalxf = polyval(p_ss_sj, signalt);
%{
figure1 = figure(f); f=f+1;
plot(sjcrk_ss(:,1), sjcrk_ss(:,2), '-', sjcrk_ss(:,1), pf_ss_sj(:,1), '-.');
hold on
plot(sjcrk_ss(:,1), sjcrk_ss(:,3), '-', sjcrk_ss(:,1), pf_ss_sj(:,2), '-.');
hold on
plot(sjcrk_ss(:,1), sjcrk_ss(:,4), '-', sjcrk_ss(:,1), pf_ss_sj(:,3), '-.');
hold off
figurename = strcat(fnumber,'sjprop_stickslip.png');
saveas(figure1, figurename);
%}
%
%{
sp = ceil(sqrt(length(sbs)));
    figure(f); f = f+1;
for i = 1:length(sbs) -1
subplot(sp, sp, i);
%figure(f); f = f+1;
    scatter1 = scatter3(sj_loc(:,1), sj_loc(:,2),sj_loc(:,3), 50, sj_area);

    %colormap(jet);
    colorbar;
        alpha(scatter1, 0.2);
    hold on
    scatter3(seg_sjcrk{i}((seg_is_shear{i}==1),2), seg_sjcrk{i}((seg_is_shear{i}==1),3), seg_sjcrk{i}((seg_is_shear{i}==1),4),'*', 'red');
    hold on
    scatter3(seg_sjcrk{i}((seg_is_shear{i}==0),2), seg_sjcrk{i}((seg_is_shear{i}==0),3), seg_sjcrk{i}((seg_is_shear{i}==0),4),'*', 'black');
    hold off
    %view([0, 0]);%
    view([0, 90]);
    xlabel('x(m)');
    ylabel('y(m)');
    
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',sprintf('SJ Area') );
%}
% Fractal dimension based on sj crk

%For each recorded cycle, calculate fractal dimension
%Perhaps too much calculation? Test fractal dimension calculation with only
%the last cycle
%cyc_plt = n0(:, 1) - n0(1, 1);
% cyc_plt = n0(:, 1);
% fact_plt = 1000;
fracDcoef_sj = zeros(length(sbs)-1, 4);
figure(f); f = f+1;
ah1 = ceil(sqrt(length(sbs)));
for j =  1:length(sbs)-1
    %X = [seg_sjcrk{window_stickslip}(:, 2), seg_sjcrk{window_stickslip}(:, 3)]; %x and y coordinates of the sj_cracks.
    clear X
    clear Y
    X = [seg_sjcrk{j}(:, 2), seg_sjcrk{j}(:, 3)]; %x and y coordinates of the sj_cracks. 
    if size(X, 1) == 1 | size(X,1) == 0
        j = j+1;
        continue
    end
    Y = pdist(X);
    Y_sq = squareform(Y); %pair-wise distances between every pair of sj_cracks.
    Z = linkage(Y);
    %dendrogram(Z, 0);
    range = linspace(min(Y)*1.5, max(Y)/2.5, 20); %Tweak the range a bit (*1.2 and /2) such that all log(Cr) and log(range) form a straight line, with no sample interval boundary problem (for instance Cr(min_dist) = 0 and Cr(large_distance) flattens)
    N_pair = zeros(length(range), 1);
    for i = 1:length(range)
         %N_pair(i) = ( sum(Y<range(i)) - length(n0(:, 1)) )/2;
         N_pair(i) =  sum(Y<range(i));
    end
    Cr = N_pair ./ ( length(sjcrk(:, 1)) .* (length(sjcrk(:, 1))-1)   );
    if size(X, 1)> 10 %More than 10 cracks in j-th time window, so Cr is actually meaningful (aka no zero element)
        [fracDcoef_sj(j, 1:3), fracDmodel] = onestep_linfit(log(range'), log(Cr));
    end
%     subplot(ah1, ah1, j);
%     scatter(log(range), log(Cr));
%     xlabel('Range');
%     ylabel('C(r)');
    %title('Fractal dimension is calculated based on the slope');
end
fracDcoef_sj(window_stickslip, 4) = 1;
csvname = strcat('fracD_sj',fnumber,'.csv');
writematrix(fracDcoef_sj, csvname);


% Fractal dimension based on MS

%For each recorded cycle, calculate fractal dimension
%Perhaps too much calculation? Test fractal dimension calculation with only
%the last cycle
%cyc_plt = n0(:, 1) - n0(1, 1);
% cyc_plt = n0(:, 1);
% fact_plt = 1000;
fracDcoef_ms = zeros(length(sbs), 4); %The last entry is based on overall. The fourth column identifies if it's during window_stickslip or overall
figure(f); f = f+1;
ah1 = ceil(sqrt(length(sbs)));
for j =  1:length(sbs)-1
    %X = [seg_sjcrk{window_stickslip}(:, 2), seg_sjcrk{window_stickslip}(:, 3)]; %x and y coordinates of the sj_cracks.
    clear X
    clear Y
    X = [seg_ms_sj{j}(:, 4), seg_ms_sj{j}(:, 5)]; %x and y coordinates of the sj_cracks.
    if size(X, 1) == 1 | size(X,1) == 0
        j = j+1;
        continue
    end
    Y = pdist(X);
    Y_sq = squareform(Y); %pair-wise distances between every pair of sj_cracks.
    Z = linkage(Y);
    %dendrogram(Z, 0);
    range = linspace(min(Y)*1.5, max(Y)/2.5, 20); %Tweak the range a bit (*1.2 and /2) such that all log(Cr) and log(range) form a straight line, with no sample interval boundary problem (for instance Cr(min_dist) = 0 and Cr(large_distance) flattens)
    N_pair = zeros(length(range), 1);
    for i = 1:length(range)
         %N_pair(i) = ( sum(Y<range(i)) - length(n0(:, 1)) )/2;
         N_pair(i) =  sum(Y<range(i));
    end
    Cr = N_pair ./ ( length(sjcrk(:, 1)) .* (length(sjcrk(:, 1))-1)   );
    if size(X, 1)> 10 %More than 10 cracks in j-th time window, so Cr is actually meaningful (aka no zero element)
        [fracDcoef_ms(j, 1:3), fracDmodel] = onestep_linfit(log(range'), log(Cr));
    end
%     subplot(ah1, ah1, j);
%     scatter(log(range), log(Cr));
%     xlabel('Range');
%     ylabel('C(r)');
    %title('Fractal dimension is calculated based on the slope');
end
%
clear X
clear Y
X = [ms_sj(:, 4), ms_sj(:, 5)]; %x and y coordinates of the sj_cracks.
Y = pdist(X);
Y_sq = squareform(Y); %pair-wise distances between every pair of sj_cracks.
Z = linkage(Y);
%dendrogram(Z, 0);
range = linspace(min(Y)*1.5, max(Y)/2.5, 20); %Tweak the range a bit (*1.2 and /2) such that all log(Cr) and log(range) form a straight line, with no sample interval boundary problem (for instance Cr(min_dist) = 0 and Cr(large_distance) flattens)
N_pair = zeros(length(range), 1);
for i = 1:length(range)
    %N_pair(i) = ( sum(Y<range(i)) - length(n0(:, 1)) )/2;
    N_pair(i) =  sum(Y<range(i));
end
Cr = N_pair ./ ( length(sjcrk(:, 1)) .* (length(sjcrk(:, 1))-1)   );
if size(X, 1)> 10 %More than 10 cracks in j-th time window, so Cr is actually meaningful (aka no zero element)
    [fracDcoef_ms(length(sbs), 1:3), fracDmodel, S] = onestep_linfit2(log(range'), log(Cr));
end
[Y,DELTA] = polyconf(fracDcoef_ms(length(sbs), 1:2),log(range'),S,'predopt','observation','simopt','on');

figure1 = figure(f); f = f+1;
scatter(log(range), log(Cr)); hold on
plot(log(range),Y+DELTA,'b--'); hold on
plot(log(range),Y-DELTA,'b--'); hold off
xlabel('Range');
ylabel('C(r)');
figurename = strcat(fnumber, 'MS_fracD.png');
saveas(figure1, figurename);
 %  
fracDcoef_ms(window_stickslip, 4) = 1;
fracDcoef_ms(end, 4) = 2;
csvname = strcat('fracD_ms',fnumber,'.csv');
writematrix(fracDcoef_ms, csvname);

% b value based on ms
sp = ceil(sqrt(length(sbs)));
temp_n = 10;
figure(f); f=f+1;
b_value = zeros(length(sbs), 1);
b_coef = zeros(length(sbs), 4);
mag_occur = zeros(1, temp_n);
for j = 1:length(sbs)-1
    if size(seg_ms_sj{j},1) == 0
        j = j+1;
        continue;
    end
    mag = seg_ms_sj{j}(:, 8);
    mag_range = linspace(min(mag), max(mag), temp_n); %/2 to clip out the nonlinear portion
    for i = 1:temp_n
        mag_occur(i) = length(find(mag>=mag_range(i)));
    end
    mag_range_linfit = mag_range( 1:min(find(mag_occur == 1)) );
    mag_occur_linfit = mag_occur( 1:min(find(mag_occur == 1)) );
    %[b_coef(j,:), b_linfit{j}] = onestep_linfit(mag_range_linfit, log10(mag_occur_linfit));
    [b_coef(j,1:3), b_linfit{j}] = onestep_linfit(mag_range, log10(mag_occur));
%     subplot(sp, sp, j);
%     scatter(mag_range, log10(mag_occur));
%     xlabel('magnitude');
%     ylabel('occurance');
end
b_value = -b_coef(:, 1); 

% overall b value, stored as the last element of b_coef array
mag = ms_sj(:, 8);
mag_range = linspace(min(mag), max(mag), temp_n);
for i = 1:temp_n
    mag_occur(i) = length(find(mag>=mag_range(i)));
end
mag_range_linfit = mag_range( 1:min(find(mag_occur == 1)) );
mag_occur_linfit = mag_occur( 1:min(find(mag_occur == 1)) );
%[b_coef(length(sbs),:), b_linfit{length(sbs)}] = onestep_linfit(mag_range_linfit, log10(mag_occur_linfit));
[b_coef(length(sbs),1:3), b_linfit{length(sbs)}, S] = onestep_linfit2(mag_range, log10(mag_occur));
[Y,DELTA] = polyconf(b_coef(length(sbs),1:2),mag_range,S,'predopt','observation','simopt','on');
figure1 = figure(f); f=f+1;
scatter(mag_range, log10(mag_occur)); hold on
plot(mag_range,Y+DELTA,'b--'); hold on
plot(mag_range,Y-DELTA,'b--'); hold off
figurename = strcat(fnumber, 'b_value.png'); saveas(figure1, figurename);
%
b_coef(end, 4) = 2;
b_coef(window_stickslip, 4) = 1;
writematrix(b_coef, strcat('b_coef_ms', fnumber, '.csv'));


%% Fractal dimension based on ae

%For each recorded cycle, calculate fractal dimension
%Perhaps too much calculation? Test fractal dimension calculation with only
%the last cycle
%cyc_plt = n0(:, 1) - n0(1, 1);
% cyc_plt = n0(:, 1);
% fact_plt = 1000;
fracDcoef_ae = zeros(length(sbs)-1, 4);
figure(f); f = f+1;
ah1 = ceil(sqrt(length(sbs)));
for j =  1:length(sbs)-1
    %X = [seg_sjcrk{window_stickslip}(:, 2), seg_sjcrk{window_stickslip}(:, 3)]; %x and y coordinates of the sj_cracks.
    clear X
    clear Y
    X = [seg_ae_sj{j}(:, 4), seg_ae_sj{j}(:, 5)]; %x and y coordinates of the sj_cracks.
    if size(X, 1) == 1
        j = j+1;
        continue
    end
    Y = pdist(X);
    Y_sq = squareform(Y); %pair-wise distances between every pair of sj_cracks.
    Z = linkage(Y);
    %dendrogram(Z, 0);
    range = linspace(min(Y)*1.5, max(Y)/2.5, 20); %Tweak the range a bit (*1.2 and /2) such that all log(Cr) and log(range) form a straight line, with no sample interval boundary problem (for instance Cr(min_dist) = 0 and Cr(large_distance) flattens)
    N_pair = zeros(length(range), 1);
    for i = 1:length(range)
         %N_pair(i) = ( sum(Y<range(i)) - length(n0(:, 1)) )/2;
         N_pair(i) =  sum(Y<range(i));
    end
    Cr = N_pair ./ ( length(sjcrk(:, 1)) .* (length(sjcrk(:, 1))-1)   );
    if size(X, 1)> 10 %More than 10 cracks in j-th time window, so Cr is actually meaningful (aka no zero element)
        [fracDcoef_ae(j, 1:3), fracDmodel] = onestep_linfit(log(range'), log(Cr));
    end
    subplot(ah1, ah1, j);
    %figure(f); f = f+1;
    scatter(log(range), log(Cr));
    xlabel('Range');
    ylabel('C(r)');
    %title('Fractal dimension is calculated based on the slope');
end

%%
% for i = 1:3
%     sj_loc(:, i) = cell2mat(C(i));
% end
% sj_radius = cell2mat(C(4));
% sj_area = cell2mat(C(5));
% sj_num = length(sj_radius);
% 
% sj_kn = 3.6e12;
% sj_ks = 1.4e12;
% for i = sj_num
%     sj_kn_dist = sj_area * sj_kn;
%     sj_ks_dist = sj_area * sj_ks;
% end
% figure(f); f = f+1;
% scatter3(sj_loc(:,1), sj_loc(:, 2), sj_loc(:,3), 50, sj_kn);
% colormap(jet);
% colorbar;
% hold on
% scatter3(seg_sjcrk{7}(:,2), seg_sjcrk{7}(:,3), seg_sjcrk{7}(:,4),'filled', 'red');
% hold off
%% Check the cs forces of those contacts during time window 6
for n = 4:7
    figure(f); f=f+1;
    for j = 1:size(cs,2)-1
        subplot(10, 10, j);
        plot(seg_cyc{n}, seg_cs_array{n}(:, j));
        hold on
        scatter(seg_sjcrk{n}(:,1), ones(length(seg_sjcrk{n}(:,1)), 1)*seg_mean_cs{n}(j), 'LineWidth', 0.5);
        hold off
    end
end

%% Check CN forces

for n = 4:7
    figure(f); f=f+1;
    for j = 1:size(cs,2)-1
        subplot(10, 10, j);
        plot(seg_cyc{n}, seg_cn_array{n}(:, j));
        hold on
        scatter(seg_sjcrk{n}(:,1), ones(length(seg_sjcrk{n}(:,1)), 1)*seg_mean{n}(j), 'LineWidth', 0.5);
        hold off
    end
end

    %%
figure(f); f = f+1;
a = [1:size(cn,2)-1]'; b = num2str(a); c=cellstr(b);
for i = 1:4
    j = i+3;
    subplot(2, 2, i);
    scatter3(ccoord(:, 1), ccoord(:, 2), ccoord(:, 3));
    text(ccoord(:, 1), ccoord(:, 2), ccoord(:, 3), c);
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
    hold off
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);   
        
end
%% Cross-correlation the segments of CS with a step function
step_func = cell(length(sbs)-1, 1);
step_match = cell(length(sbs)-1, 1);
step_rise = cell(length(sbs)-1, 1);
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        clear temp
        temp = seg_cs_array{j}(:, i);
        %if temp'*temp ~= 0
            %temp = temp/sqrt(temp'*temp); %Need to normalize CS and the resultant step function so that the cross correlation results are comparable. Otherwise, one CS with large values but less step-function like behaviour can still result in higher xcorr value.
            ipt = findchangepts(temp);
            step_prejump = mean(temp(1:ipt));
            step_postjump = mean(temp(ipt:end));

            step_func{j}(:, i) = [ones(ipt, 1)*step_prejump; ones(length(temp)-ipt, 1)*step_postjump ];
            step_match{j}(i) = max( xcorr(temp, step_func{j}(:, i)) );
            step_rise{j}(i) = 1*(step_prejump < step_postjump);
%         else
%             step_func{j}(:, i) = zeros(size(temp, 1), 1);
%             step_match{j}(i) = 0;
%             step_rise{j}(i) = 0;
%         end
    end
end
%{
for j = 1:length(sbs)-1
    for i = 1:num_usefulsensor
        clear temp
        temp = seg_cs_array{j}(:, i);
        if temp'*temp ~= 0
            temp = temp/sqrt(temp'*temp); %Need to normalize CS and the resultant step function so that the cross correlation results are comparable. Otherwise, one CS with large values but less step-function like behaviour can still result in higher xcorr value.
            ipt = findchangepts(temp);
            step_prejump = mean(temp(1:ipt));
            step_postjump = mean(temp(ipt:end));

            step_func{j}(:, i) = [ones(ipt, 1)*step_prejump; ones(length(temp)-ipt, 1)*step_postjump ];
            step_match{j}(i) = max( xcorr(temp, step_func{j}(:, i)) );
            step_rise{j}(i) = 1*(step_prejump < step_postjump);
        else
            step_func{j}(:, i) = zeros(size(temp, 1), 1);
            step_match{j}(i) = 0;
            step_rise{j}(i) = 0;
        end
    end
end
%}

figure(f); f = f+1;
for j = 1:num_usefulsensor
    subplot(10, 10, j);
    plot(seg_cyc{5}, step_func{5}(:, j));
    hold on
    plot(seg_cyc{5}, seg_cs_array{5}(:, j)/sqrt(seg_cs_array{5}(:, j)'*seg_cs_array{5}(:, j)));
    hold off
end

figure(f); f = f+1;
a = [1:size(cn,2)-1]'; b = num2str(a); c=cellstr(b);

for i = 1:4
    j = i+3;
    subplot(2, 2, i);
    scatter3(ccoord(:, 1), ccoord(:, 2), ccoord(:, 3), 50, step_rise{j});
    colorbar;
    text(ccoord(:, 1), ccoord(:, 2), ccoord(:, 3), c);
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
     %text(ccoord(:, 1), ccoord(:, 2), ccoord(:, 3), c1);
    hold off
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);          
end

a1 = [1:16]'; b1 = num2str(a1); c1=cellstr(b1);
figure(f); f = f+1;
scatter3(seg_sjcrk{7}(:,2), seg_sjcrk{7}(:,3), seg_sjcrk{7}(:,4),'*');
text(seg_sjcrk{7}(:,2), seg_sjcrk{7}(:,3), seg_sjcrk{7}(:,4), c1);

figure(f); f = f+1;
a = [1:size(cn,2)-1]'; b = num2str(a); c=cellstr(b);
for i = 1:4
    j = i+3;
    subplot(2, 2, i);
    scatter3(ccoord(:, 1), ccoord(:, 2), ccoord(:, 3), 50, step_match{j});
    %colormap(jet);
    colorbar;
    %a = step_match{j}'; b = num2str(a); c=cellstr(b);
    text(ccoord(:, 1), ccoord(:, 2), ccoord(:, 3), c);
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
    hold off
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);          
end
%% Multi-segment linear fit to the CS
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


%% Check if particle scale stress drop can be representative of grain scale stick slip
% Force drop multiplied by deformation gives energy
figure(f); f = f+1;
n = 10;
sp = ceil(sqrt(length(sbs)));
for j = n:n
    for i = 1:num_usefulsensor
        subplot(sp, sp, i);
        plot(seg_cyc{j}, seg_cdx_array{j}(:,i));
        hold on
        plot(seg_cyc{j}(1:ladder_idx_cdx(j)) , ladder_prejump_cdx{j}(:, i));
        hold on
        plot(seg_cyc{j}(ladder_idx_cdx(j):end) , ladder_postjump_cdx{j}(:, i));
        hold off
    end
end

figure(f); f = f+1;
sp = ceil(sqrt(length(sbs)));
for j = n:n
    for i = 1:num_usefulsensor
        subplot(sp, sp, i);
        plot(seg_cyc{j}, seg_cs_array{j}(:,i));
    end
end
%% k mean analysis
% cn_norm = zeros(size(cn));
% for i = 1:num_usefulsensor
%     cn_norm(:, i) = cn(:, i)/ norm(cn(:,i));
% end
num_cluster = 4;
% idx = kmeans(cn_norm',num_cluster);
% 
% idx_ori = 1:1:num_usefulsensor;
% idx_calc = [idx(1:end-1), idx_ori'];
% for i = 1:num_usefulsensor
%     idx_sort = sortrows(idx_calc, [1, 2]);
% end

% figure(f); f=f+1;
% for i = 1:num_usefulsensor
%     subplot(10, 10, i);
%     plot(t, cn_norm(:, idx_sort(i, 2)));
%      ylabel('CN (N)');
%     xlabel('time');
% end

%     cn_cluster = cell(num_cluster, 1);
%     for j = 1:num_cluster
%         temp = idx_sort(:,1) == j;
%         temp2 = idx_sort(temp, 2);
%         cn_cluster{j, 1} = cn(:, temp2);
%     end
    
    
 [idx_sort, cn_cluster] = cluster_shape(cn, num_cluster, num_usefulsensor);
    
for i = 1:num_cluster
    figure(f); f =f+1;
    for j = 1:size(cn_cluster{i, 1}, 2)
        subplot(6, 6, j);
        plot(t, cn_cluster{i}(:, j));
        ylabel('CN (N)');
        xlabel('time');
    end
end

figure(f); f =f+1;
for i = 1:num_cluster
        subplot(3, 3, i);
        temp = idx_sort(idx_sort(:,1) == i, 2);
        scatter3(ccoord(temp, 1), ccoord(temp, 2), ccoord(temp, 3));
       
        view([0, 90]);
        xlim([-20e-3, 20e-3]);
         ylim([-20e-3, 20e-3]);       
end


%% Cluster the segment of contact force that contains stick sl

[idx_sort7, cn_cluster7] = cluster_shape_nonnormal(seg_cn_array{7}, num_cluster, num_usefulsensor);
for i = 1:num_cluster
    figure(f); f =f+1;
    for j = 1:size(cn_cluster7{i, 1}, 2)
        subplot(7, 10, j);
        plot(seg_t{7}, cn_cluster7{i}(:, j));
    end
end



figure(f); f =f+1;
for i = 1:num_cluster
        subplot(3, 5, i);
        temp = idx_sort7(idx_sort7(:,1) == i, 2);
        scatter3(ccoord(temp, 1), ccoord(temp, 2), ccoord(temp, 3));
        hold on
    scatter3(seg_sjcrk{7}(:,2), seg_sjcrk{7}(:,3), seg_sjcrk{7}(:,4),'*');
    hold off
        view([0, 90]);
        xlim([-20e-3, 20e-3]);
         ylim([-20e-3, 20e-3]);       
end

%%
%%
[idx_cs, cs_cluster] = cluster_shape(cs, num_cluster, num_usefulsensor);
idx_cs_cell = cell(length(sbs), 1);
idx_cs_cell{end} = idx_cs;
cs_cluster_cell = cell(length(sbs), 1);
cs_cluster_cell{end} = cs_cluster;

for i = 1:num_cluster
    figure(f); f =f+1;
    for j = 1:size(cs_cluster{i, 1}, 2)
        subplot(6, 6, j);
        plot(t, cs_cluster{i}(:, j));
                ylabel('CS (N)');
        xlabel('time');
    end
end

figure(f); f =f+1;
for i = 1:num_cluster
        subplot(3, 3, i);
        temp = idx_cs(idx_sort(:,1) == i, 2);
        scatter3(ccoord(temp, 1), ccoord(temp, 2), ccoord(temp, 3));
%         hold on
%     scatter3(seg_sjcrk{7}(:,2), seg_sjcrk{7}(:,3), seg_sjcrk{7}(:,4),'*');
%     hold off
        view([0, 90]);
        xlim([-20e-3, 20e-3]);
         ylim([-20e-3, 20e-3]);       
end
%%
%close all

for i = 1:length(sbs)-1
    [idx_cs_temp, cs_cluster_temp] = cluster_shape(seg_cs_array{i}, num_cluster, num_usefulsensor);
    idx_cs_cell{i} = idx_cs_temp;
    cs_cluster_cell{i} = cs_cluster_temp;
end


%%
n = 16;%window_stickslip; %10;
clear temp
for i = 1:num_cluster
    figure(f); f =f+1;
    temp = idx_cs_cell{n}(idx_cs_cell{n}(:,1)==i, 2);
    if size(temp)~=0, sp1=ceil(sqrt(size(temp,1)));else sp1=1;end, sp2=sp1;
    for j = 1:size(cs_cluster_cell{n}{i, 1}, 2)
        subplot(sp1, sp2, j);
        plot(seg_cyc{n}, cs_cluster_cell{n}{i}(:, j));
        hold on
        
        scatter(seg_sjcrk{n}(:,1), ones(length(seg_sjcrk{n}(:,1)), 1)*seg_mean_cs{n}(temp(j)), 'LineWidth', 0.5);
        hold off
        xlabel('PFC cycles');
        ylabel('CS (N)');
    end
end

sp1 = ceil(sqrt(num_cluster)); sp2 = ceil(sqrt(num_cluster));
figure(f); f =f+1; %Plot cracks and CS during time window 7. Step up, neutral and down in CS are marked as green, yellow and red
for i = 1:num_cluster
    subplot(sp1, sp2, i);
    temp = idx_cs_cell{n}(idx_cs_cell{n}(:,1) == i, 2);
    a = [1:length(temp)]'; b = num2str(a); c=cellstr(b);
    scatter3(ccoord(temp, 1), ccoord(temp, 2), ccoord(temp, 3));
    text(ccoord(temp, 1), ccoord(temp, 2), ccoord(temp, 3), c);
    hold on
    scatter3(seg_sjcrk{n}(:,2), seg_sjcrk{n}(:,3), seg_sjcrk{n}(:,4),'*');
    hold off
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);
    xlabel('x(m)');
    ylabel('y(m)');
    clear a
    clear b
    clear c
end

sp1 = ceil(sqrt(num_cluster)); sp2 = ceil(sqrt(num_cluster));
figure(f); f =f+1; %Plot cracks and CS during time window n. Step up, neutral and down in CS are marked as green, yellow and red
for i = 1:num_cluster
    subplot(sp1, sp2, i);
    temp = idx_cs_cell{n}(idx_cs_cell{n}(:,1) == i, 2);
    temp1 = find(ladder_jump_cs(n,temp)>0);
    scatter3(ccoord(temp(temp1), 1), ccoord(temp(temp1), 2), ccoord(temp(temp1), 3), 'g'); %CS step up
    hold on
    temp1 = find(ladder_jump_cs(n,temp)<0);
    scatter3(ccoord(temp(temp1), 1), ccoord(temp(temp1), 2), ccoord(temp(temp1), 3), 'r'); %CS step down
    hold on
    temp1 = find(ladder_jump_cs(n,temp)==0);
    scatter3(ccoord(temp(temp1), 1), ccoord(temp(temp1), 2), ccoord(temp(temp1), 3), 'y'); %CS no change, likely a dead contact
    hold on
    scatter3(seg_sjcrk{n}(:,2), seg_sjcrk{n}(:,3), seg_sjcrk{n}(:,4),'*');
    hold off
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);
    xlabel('x(m)');
    ylabel('y(m)');
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',sprintf('Cracks and CS during time window %d. Step up, neutral and down in CS are marked as green, yellow and red',n), 'EdgeColor', 'none');

sp1 = ceil(sqrt(num_cluster)); sp2 = ceil(sqrt(num_cluster));
figure(f); f =f+1; %Plot cracks and CS during time window n. positive, 0, negative slope of prejump marked as green, yellow and red
for i = 1:num_cluster
    subplot(sp1, sp2, i);
    temp = idx_cs_cell{n}(idx_cs_cell{n}(:,1) == i, 2);
    temp1 = find(coef_prejump_cs(n, temp, 1)>0);
    scatter3(ccoord(temp(temp1), 1), ccoord(temp(temp1), 2), ccoord(temp(temp1), 3), 'g'); %CS step up
    hold on
    temp1 = find(coef_prejump_cs(n, temp, 1)<0);
    scatter3(ccoord(temp(temp1), 1), ccoord(temp(temp1), 2), ccoord(temp(temp1), 3), 'r'); %CS step down
    hold on
    temp1 = find(coef_prejump_cs(n, temp, 1)==0);
    scatter3(ccoord(temp(temp1), 1), ccoord(temp(temp1), 2), ccoord(temp(temp1), 3), 'y'); %CS no change, likely a dead contact
    hold on
    scatter3(seg_sjcrk{n}(:,2), seg_sjcrk{n}(:,3), seg_sjcrk{n}(:,4),'*');
    hold off
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);
    xlabel('x(m)');
    ylabel('y(m)');
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',sprintf('Cracks and CS during time window %d. Positive, zero, and negative slop in CS pre-jump are marked as green, yellow and red',n), 'EdgeColor', 'none');

sp1 = ceil(sqrt(num_cluster)); sp2 = ceil(sqrt(num_cluster));
figure(f); f =f+1; %Plot cracks and CS during time window n. positive, 0, negative slope of postjump marked as green, yellow and red
for i = 1:num_cluster
    subplot(sp1, sp2, i);
    temp = idx_cs_cell{n}(idx_cs_cell{n}(:,1) == i, 2);
    temp1 = find(coef_postjump_cs(n, temp, 1)>0);
    scatter3(ccoord(temp(temp1), 1), ccoord(temp(temp1), 2), ccoord(temp(temp1), 3), 'g'); %CS step up
    hold on
    temp1 = find(coef_postjump_cs(n, temp, 1)<0);
    scatter3(ccoord(temp(temp1), 1), ccoord(temp(temp1), 2), ccoord(temp(temp1), 3), 'r'); %CS step down
    hold on
    temp1 = find(coef_postjump_cs(n, temp, 1)==0);
    scatter3(ccoord(temp(temp1), 1), ccoord(temp(temp1), 2), ccoord(temp(temp1), 3), 'y'); %CS no change, likely a dead contact
    hold on
    scatter3(seg_sjcrk{n}(:,2), seg_sjcrk{n}(:,3), seg_sjcrk{n}(:,4),'*');
    hold off
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);
    xlabel('x(m)');
    ylabel('y(m)');
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',sprintf('Cracks and CS during time window %d. Positive, zero, and negative slop in CS post-jump are marked as green, yellow and red',n), 'EdgeColor', 'none');


%% Run statistics analysis on jump up or down. First is to check ladder fits in various time windows and see how mnay different levels of goodness of fits can be defined.
%sp1 = 3; sp2 = 3;
n = 10;
figure(f); f =f+1; %Plot cracks and CS during time window n. Step up, neutral and down in CS are marked as green, yellow and red
for i = 1:1
    %subplot(sp1, sp2, i);
    %temp = idx_cs_cell{n}(idx_cs_cell{n}(:,1) == i, 2);
    temp1 = find(ladder_jump_cs(n,:)>0);
    scatter3(ccoord(temp1, 1), ccoord(temp1, 2), ccoord(temp1, 3), 'g'); %CS step up
    hold on
    temp1 = find(ladder_jump_cs(n,:)<0);
    scatter3(ccoord(temp1, 1), ccoord(temp1, 2), ccoord(temp1, 3), 'r'); %CS step down
    hold on
    temp1 = find(ladder_jump_cs(n,:)==0);
    scatter3(ccoord(temp1, 1), ccoord(temp1, 2), ccoord(temp1, 3), 'y'); %CS no change, likely a dead contact
    hold on
    scatter3(seg_sjcrk{n}(:,2), seg_sjcrk{n}(:,3), seg_sjcrk{n}(:,4),'*');
    hold off
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);
    xlabel('x(m)');
    ylabel('y(m)');
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',sprintf('Cracks and CS during time window %d. Step up, neutral and down in CS are marked as green, yellow and red',n), 'EdgeColor', 'none');

figure(f); f = f+1;
sp1 = 10; sp2 = 10;
for j = n%1:length(sbs)-1
    for i = 1:num_usefulsensor
        subplot(sp1, sp2, i);
        plot(seg_cyc{j}, seg_cs_array{j}(:,i));
        hold on
        plot(seg_cyc{j}(1:ladder_idx_cs(j)) , ladder_prejump_cs{j}(:, i));
        hold on
        plot(seg_cyc{j}(ladder_idx_cs(j):end) , ladder_postjump_cs{j}(:, i));
        hold off
    end
end

%% Particle displacement analysis 
n = 10;
sp1 = 10; sp2 = 10;
figure(f); f = f+1;
for j = n%1:length(sbs)-1
    for i = 1:num_usefulsensor
        subplot(sp1, sp2, i);
        plot(seg_cyc{j}, seg_bdz_array{j}(:,i));
        hold on
        plot(seg_cyc{j}(1:ladder_idx_bdz(j)) , ladder_prejump_bdz{j}(:, i));
        hold on
        plot(seg_cyc{j}(ladder_idx_bdz(j):end) , ladder_postjump_bdz{j}(:, i));
        hold off
    end
end
figure(f); f = f+1;
for i = 1:1
    %subplot(sp1, sp2, i);
    %temp = idx_cs_cell{n}(idx_cs_cell{n}(:,1) == i, 2);
    temp1 = find(ladder_jump_bdz(n,:)>0);
    scatter3(ccoord(temp1, 1), ccoord(temp1, 2), ccoord(temp1, 3), 'g'); %CS step up
    hold on
    temp1 = find(ladder_jump_bdz(n,:)<0);
    scatter3(ccoord(temp1, 1), ccoord(temp1, 2), ccoord(temp1, 3), 'r'); %CS step down
    hold on
    temp1 = find(ladder_jump_bdz(n,:)==0);
    scatter3(ccoord(temp1, 1), ccoord(temp1, 2), ccoord(temp1, 3), 'y'); %CS no change, likely a dead contact
    hold on
    scatter3(seg_sjcrk{n}(:,2), seg_sjcrk{n}(:,3), seg_sjcrk{n}(:,4),'*');
    hold off
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);
    xlabel('x(m)');
    ylabel('y(m)');
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',sprintf('BDZ during time window %d. Step up, neutral and down in CS are marked as green, yellow and red',n), 'EdgeColor', 'none');
%% Plot entire force history based on cs_cluster7 (and the respective idx_cs7)
figure(f); f =f+1;
for i = 1:num_cluster
    figure(f); f =f+1;
    temp = cs(:, (idx_cs7(idx_cs7(:,1)==i, 2)));
    for j = 1:size(cs_cluster7{i, 1}, 2)
        subplot(7, 10, j);
        plot(t, temp(:, j) );
        
    end
end

%% Crack cluster analysis
%focus on 7 for now
idx_crk7 = kmeans([seg_sjcrk{7}(:, 2:4)],3);
figure(f); f= f+1;   
for i = 1:3
    subplot(2,2, i);
    scatter3(seg_sjcrk{7}(idx_crk7==i,2), seg_sjcrk{7}(idx_crk7==i,3), seg_sjcrk{7}(idx_crk7==i,4),'*');
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);    
end

%%

temp1 = seg_sjcrk{7}(idx_crk7==1,2);
temp2 = seg_sjcrk{7}(idx_crk7==1,3);
temp3 = seg_sjcrk{7}(idx_crk7==1,4);
boundu = [max(temp1); max(temp2); max(temp3)];
boundl = [min(temp1); min(temp2); min(temp3)];

clear temp1
clear temp2
clear temp3
temp1 = find(ccoord(:,1)<boundu(1)*1.2 & ccoord(:,1)>boundl(1)*1.2);
temp2 = find(ccoord(:,2)<boundu(2)*1.2 & ccoord(:,2)>boundl(2)*1.2);
temp3 = find(ccoord(:,3)<boundu(3)*1.2 & ccoord(:,3)>boundl(3)*1.2);
temp4 = intersect(temp1, temp2);
cs_nearcrk7_ind = intersect(temp4, temp3);
cs_nearcrk7 = sjcrk(cs_nearcrk7_ind, :);
figure(f); f= f+1;
for i = 1:length(cs_nearcrk7)
    subplot(5, 5, i);
    plot(t, cs(:, cs_nearcrk7_ind(i)));
end
figure(f); f= f+1;
for i = 1:length(cs_nearcrk7)
    subplot(5, 5, i);
    plot(seg_t{7}, seg_cs_array{7}(:, cs_nearcrk7_ind(i)));
end
figure(f); f= f+1;
    scatter3(ccoord(cs_nearcrk7_ind,1), ccoord(cs_nearcrk7_ind,2), ccoord(cs_nearcrk7_ind,3));
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);    

    
%% measurement of flatness: std from a flat line
    seg_cs_flat = cell(length(sbs)-1, 1);
%     for j = 1:length(sbs)-1
%         %temp_length = length([ sbs(j): sbs(j+1) ]);
%         temp_model_array = zeros(1, num_usefulsensor);
%         for i = 1:num_usefulsensor
%             temp_sd = sum(...
%                 ( seg_cs_array{j}(:,i) - seg_cs_array{j}(1,i) ).^2 ...
%                 );        
%             temp_model_array(1, i) = temp_sd;
%         end
%         seg_cs_flat{j} = temp_model_array;
%     end
  
 figure(f); f= f+1;
 plot(seg_cs_flat{7});
 hold on
 scatter(cs_nearcrk7_ind, seg_cs_flat{7}(cs_nearcrk7_ind), 'o');
 hold off
%% Create a grid that I can plot the data for all the samples
mv_W = 40e-3;
mv_H = 107e-3;

x = -mv_W/2 : 1e-4: mv_W/2;
y = -mv_W/2 : 1e-4: mv_W/2;
[X,Y] = meshgrid(x, y);
Z = zeros(length(x), length(y));

figure(f); f=f+1;
scatter3(ccoord(:,1),ccoord(:,2),ccoord(:,3));



%% Crack propogation along the fault through out time windows

% Crack propogation along the fault through out time windows
%close all

%[sjcrk, crk_cyc, is_sj_crk] = openfile_crk(cell2mat(fname(17)));

figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)));
for j = 1:length(sbs)-1    
%     for i = 1:num_usefulsensor
%         indx = find(abs(x - ccoord(i, 1)) < 1e-6);
%         indy = find(abs(y - ccoord(i, 2)) < 1e-6);
%         %Z(indx, indy) = seg_cn_coef_array{j}(3,i);
%         %Z(indx, indy) = seg_cs_sd{j}(1,i);
%         Z(indx, indy) = seg_mean_cs{j}(i);
%     end
    subplot(sp, sp, j);
    %figure(f); f=f+1;
    scatter3(ccoord(:, 1), ccoord(:,2), ccoord(:, 3),50, seg_mean_cs{j},'filled');
    %scatter3(X(Z~=0), Y(Z~=0), Z(Z~=0),50, Z(Z~=0),'filled');
    colormap(jet);
    colorbar;
    view([0, 90]); %xy plane
    %view([90, 0]); %yz plane
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    title(sprintf('Average contact shear forces during time window %d', j));
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
    hold off
end
%%
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)));
for j = 1:length(sbs)-1    
%     for i = 1:num_usefulsensor
%         indx = find(abs(x - ccoord(i, 1)) < 1e-6);
%         indy = find(abs(y - ccoord(i, 2)) < 1e-6);
%         %Z(indx, indy) = seg_cn_coef_array{j}(3,i);
%         %Z(indx, indy) = seg_cs_sd{j}(1,i);
%         Z(indx, indy) = seg_mean_cs{j}(i);
%     end
    subplot(sp, sp, j);
    %figure(f); f=f+1;
    scatter3(ccoord(:, 1), ccoord(:,2), ccoord(:, 3),50, seg_mean{j},'filled');
    %scatter3(X(Z~=0), Y(Z~=0), Z(Z~=0),50, Z(Z~=0),'filled');
    colormap(jet);
    colorbar;
    view([0, 90]); %xy plane
    %view([90, 0]); %yz plane
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    title(sprintf('Average contact normal forces during time window %d', j));
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
    hold off
end
%%
for j = 1:length(sbs)-1    
    for i = 1:num_usefulsensor
        indx = find(abs(x - ccoord(i, 1)) < 1e-6);
        indy = find(abs(y - ccoord(i, 2)) < 1e-6);
        %Z(indx, indy) = seg_cn_coef_array{j}(3,i);
        %Z(indx, indy) = seg_cs_sd{j}(1,i);
        Z(indx, indy) = seg_mean{j}(i);
    end
    figure(f); f=f+1;
    scatter3(ccoord(:, 1), ccoord(:,2), ccoord(:, 3),50, seg_mean_cs{j},'filled');
    %scatter3(X(Z~=0), Y(Z~=0), Z(Z~=0),50, Z(Z~=0),'filled');
    colormap(jet);
    colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    title(sprintf('Average contact shear forces during time window %d', j));
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
    hold off
end
%%
for j = 1:length(sbs)-1    
    for i = 1:num_usefulsensor
        indx = find(abs(x - ccoord(i, 1)) < 1e-6);
        indy = find(abs(y - ccoord(i, 2)) < 1e-6);
        %Z(indx, indy) = seg_cn_coef_array{j}(3,i);
        %Z(indx, indy) = seg_cs_sd{j}(1,i);
        Z(indx, indy) = seg_mean{j}(i);
    end
    figure(f); f=f+1;
    scatter3(ccoord(:, 1), ccoord(:,2), ccoord(:, 3),50, seg_mean_cs{j},'filled');
    %scatter3(X(Z~=0), Y(Z~=0), Z(Z~=0),50, Z(Z~=0),'filled');
    colormap(jet);
    colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    title(sprintf('Average contact shear stresses during time window %d', j));
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
    hold off
end
%%
 figure(f); f=f+1;
 sp = ceil(sqrt(length(sbs)));
for j = 1:length(sbs)-1    
    %figure(f); f=f+1;
    subplot(sp, sp, j);
    scatter3(ccoord(:, 1), ccoord(:,2), ccoord(:, 3),50, ladder_relative_cs(j,:),'filled');   
    colormap(jet);
    colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    title(sprintf('Jump in force relative to average force during time window %d', j));
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
    hold off
end
%%
%%
 figure(f); f=f+1;
 sp = ceil(sqrt(length(sbs)));
for j = 1:length(sbs)-1    
    %figure(f); f=f+1;
    subplot(sp, sp, j);
    scatter3(ccoord(:, 1), ccoord(:,2), ccoord(:, 3),50, ladder_relative_cs(j,:),'filled');   
    colormap(jet);
    colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    title(sprintf('Jump in force relative to average force during time window %d', j));
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
    hold off
end
%%
 figure(f); f=f+1;
 sp = ceil(sqrt(length(sbs)));
for j = 1:length(sbs)-1    
    %figure(f); f=f+1;
    subplot(sp, sp, j);
    scatter3(ccoord(:, 1), ccoord(:,2), ccoord(:, 3),50, seg_mean_bdz{j},'filled');
    colormap(jet);
    colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    title(sprintf('averae bdz during time window %d', j));
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
    hold off
end

%%
for j = 1:length(sbs)-1    

    figure(f); f=f+1;
    scatter3(ccoord(:, 1), ccoord(:,2), ccoord(:, 3),50, coef_postjump_cs(j,:, 1),'filled');
    colormap(jet);
    colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    title(sprintf('Post jump slope in force during time window %d', j));
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
    hold off
end

%%
%%
for j = 1:length(sbs)-1    

    figure(f); f=f+1;
    scatter3(ccoord(:, 1), ccoord(:,2), ccoord(:, 3),50, ladder_jump_bdz(j,:),'filled');
    colormap(jet);
    colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    title(sprintf('Jump in z displacement during time window %d', j));
    hold on
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'*');
    hold off
end

%% rise or fall (positive and negative slopes) of bdz and cs prejump
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1
    subplot(sp, sp, j);
    scatter3(ccoord( (coef_prejump_bdz(j,:,1)> 0), 1), ccoord((coef_prejump_bdz(j,:,1)> 0),2), ccoord((coef_prejump_bdz(j,:,1)> 0), 3) , 'b');
    hold on
    scatter3(ccoord( (coef_prejump_bdz(j,:,1)<= 0), 1), ccoord((coef_prejump_bdz(j,:,1)<= 0),2), ccoord((coef_prejump_bdz(j,:,1)<= 0), 3) , 'r');
    hold on
    scatter3(ccoord( (coef_prejump_cs(j,:,1)> 0), 1), ccoord((coef_prejump_cs(j,:,1)> 0),2), ccoord((coef_prejump_cs(j,:,1)> 0), 3) , '*g');
     hold on
     scatter3(ccoord( (coef_prejump_cs(j,:,1)<= 0), 1), ccoord((coef_prejump_cs(j,:,1)<= 0),2), ccoord((coef_prejump_cs(j,:,1)<= 0), 3) , '*y');
    hold off
    view([0, 90]);
    %view([0, 0]);
    xlabel('x(mm)');
    ylabel('y(mm)');
end
%% Jump in bdz (values including both magnitudes and signs)
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1    
    subplot(sp, sp, j);
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'k*');
    hold on
    scatter3(ccoord(ladder_jump_bdz(j,:)>0, 1), ccoord(ladder_jump_bdz(j,:)>0,2), ccoord(ladder_jump_bdz(j,:)>0, 3), 'g');
    hold on
    scatter3(ccoord(ladder_jump_bdz(j,:)==0, 1), ccoord(ladder_jump_bdz(j,:)==0,2), ccoord(ladder_jump_bdz(j,:)==0, 3), 'y');
    hold on
    scatter3(ccoord(ladder_jump_bdz(j,:)<0, 1), ccoord(ladder_jump_bdz(j,:)<0,2), ccoord(ladder_jump_bdz(j,:)<0, 3), 'r');
    %colormap(jet);
    %colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    %title(sprintf('Jump in z displacement during time window %d', j));
    hold off
end

%% Relative jump in bdz (values including both magnitudes and signs). It takes care of how it jumps/drops relative to the average values, so negative drop is equivalent to positive jump
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1    
    subplot(sp, sp, j);
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'k*');
    hold on
    scatter3(ccoord(ladder_relative_bdz(j,:)>0, 1), ccoord(ladder_relative_bdz(j,:)>0,2), ccoord(ladder_relative_bdz(j,:)>0, 3), 'g');
    hold on
    scatter3(ccoord(ladder_relative_bdz(j,:)==0, 1), ccoord(ladder_relative_bdz(j,:)==0,2), ccoord(ladder_relative_bdz(j,:)==0, 3), 'y');
    hold on
    scatter3(ccoord(ladder_relative_bdz(j,:)<0, 1), ccoord(ladder_relative_bdz(j,:)<0,2), ccoord(ladder_relative_bdz(j,:)<0, 3), 'r');
    %colormap(jet);
    %colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    %title(sprintf('Jump in z displacement during time window %d', j));
    hold off
end
%% Jump in bdy (values including both magnitudes and signs)
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1    
    subplot(sp, sp, j);
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'k*');
    hold on
    scatter3(ccoord(ladder_jump_bdy(j,:)>0, 1), ccoord(ladder_jump_bdy(j,:)>0,2), ccoord(ladder_jump_bdy(j,:)>0, 3), 'g');
    hold on
    scatter3(ccoord(ladder_jump_bdy(j,:)==0, 1), ccoord(ladder_jump_bdy(j,:)==0,2), ccoord(ladder_jump_bdy(j,:)==0, 3), 'y');
    hold on
    scatter3(ccoord(ladder_jump_bdy(j,:)<0, 1), ccoord(ladder_jump_bdy(j,:)<0,2), ccoord(ladder_jump_bdy(j,:)<0, 3), 'r');
    %colormap(jet);
    %colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    %title(sprintf('Jump in z displacement during time window %d', j));
    hold off
end
%% Relative jump in bdy (values including both magnitudes and signs)
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1    
    subplot(sp, sp, j);
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'k*');
    hold on
    scatter3(ccoord(ladder_relative_bdy(j,:)>0, 1), ccoord(ladder_relative_bdy(j,:)>0,2), ccoord(ladder_relative_bdy(j,:)>0, 3), 'g');
    hold on
    scatter3(ccoord(ladder_relative_bdy(j,:)==0, 1), ccoord(ladder_relative_bdy(j,:)==0,2), ccoord(ladder_relative_bdy(j,:)==0, 3), 'y');
    hold on
    scatter3(ccoord(ladder_relative_bdy(j,:)<0, 1), ccoord(ladder_relative_bdy(j,:)<0,2), ccoord(ladder_relative_bdy(j,:)<0, 3), 'r');
    %colormap(jet);
    %colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    %title(sprintf('Jump in z displacement during time window %d', j));
    hold off
end
%% Jump in bdx (values including both magnitudes and signs)
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1    
    subplot(sp, sp, j);
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'k*');
    hold on
    scatter3(ccoord(ladder_jump_bdx(j,:)>0, 1), ccoord(ladder_jump_bdx(j,:)>0,2), ccoord(ladder_jump_bdx(j,:)>0, 3), 'g');
    hold on
    scatter3(ccoord(ladder_jump_bdx(j,:)==0, 1), ccoord(ladder_jump_bdx(j,:)==0,2), ccoord(ladder_jump_bdx(j,:)==0, 3), 'y');
    hold on
    scatter3(ccoord(ladder_jump_bdx(j,:)<0, 1), ccoord(ladder_jump_bdx(j,:)<0,2), ccoord(ladder_jump_bdx(j,:)<0, 3), 'r');
    %colormap(jet);
    %colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    %title(sprintf('Jump in z displacement during time window %d', j));
    hold off
end
%% Relative jump in bdx (values including both magnitudes and signs)
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1    
    subplot(sp, sp, j);
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'k*');
    hold on
    scatter3(ccoord(ladder_relative_bdx(j,:)>0, 1), ccoord(ladder_relative_bdx(j,:)>0,2), ccoord(ladder_relative_bdx(j,:)>0, 3), 'g');
    hold on
    scatter3(ccoord(ladder_relative_bdx(j,:)==0, 1), ccoord(ladder_relative_bdx(j,:)==0,2), ccoord(ladder_relative_bdx(j,:)==0, 3), 'y');
    hold on
    scatter3(ccoord(ladder_relative_bdx(j,:)<0, 1), ccoord(ladder_relative_bdx(j,:)<0,2), ccoord(ladder_relative_bdx(j,:)<0, 3), 'r');
    %colormap(jet);
    %colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    %title(sprintf('Jump in z displacement during time window %d', j));
    hold off
end
%% Jump in cs (values including both magnitudes and signs)
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1    
    subplot(sp, sp, j);
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'k*');
    hold on
    scatter3(ccoord(ladder_jump_cs(j,:)>0, 1), ccoord(ladder_jump_cs(j,:)>0,2), ccoord(ladder_jump_cs(j,:)>0, 3), 'g');
    hold on
    scatter3(ccoord(ladder_jump_cs(j,:)==0, 1), ccoord(ladder_jump_cs(j,:)==0,2), ccoord(ladder_jump_cs(j,:)==0, 3), 'y');
    hold on
    scatter3(ccoord(ladder_jump_cs(j,:)<0, 1), ccoord(ladder_jump_cs(j,:)<0,2), ccoord(ladder_jump_cs(j,:)<0, 3), 'r');
    %colormap(jet);
    %colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    %title(sprintf('Jump in z displacement during time window %d', j));
    hold off
end
%% Relative jump in cs 
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1    
    subplot(sp, sp, j);
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'k*');
    hold on
    scatter3(ccoord(ladder_relative_cs(j,:)>0, 1), ccoord(ladder_relative_cs(j,:)>0,2), ccoord(ladder_relative_cs(j,:)>0, 3), 'g');
    hold on
    scatter3(ccoord(isnan(ladder_relative_cs(j,:)), 1), ccoord((isnan(ladder_relative_cs(j,:))),2), ccoord((isnan(ladder_relative_cs(j,:))), 3), 'y');
    hold on
    scatter3(ccoord(ladder_relative_cs(j,:)<0, 1), ccoord(ladder_relative_cs(j,:)<0,2), ccoord(ladder_relative_cs(j,:)<0, 3), 'r');
    %colormap(jet);
    %colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    %title(sprintf('Jump in z displacement during time window %d', j));
    hold off
end
%% Jump in cn (values including both magnitudes and signs)
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1    
    subplot(sp, sp, j);
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'k*');
    hold on
    scatter3(ccoord(ladder_jump_cn(j,:)>0, 1), ccoord(ladder_jump_cn(j,:)>0,2), ccoord(ladder_jump_cn(j,:)>0, 3), 'g');
    hold on
    scatter3(ccoord(ladder_jump_cn(j,:)==0, 1), ccoord(ladder_jump_cn(j,:)==0,2), ccoord(ladder_jump_cn(j,:)==0, 3), 'y');
    hold on
    scatter3(ccoord(ladder_jump_cn(j,:)<0, 1), ccoord(ladder_jump_cn(j,:)<0,2), ccoord(ladder_jump_cn(j,:)<0, 3), 'r');
    %colormap(jet);
    %colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    %title(sprintf('Jump in z displacement during time window %d', j));
    hold off
end
%% Relative jump in cn 
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1    
    subplot(sp, sp, j);
    scatter3(seg_sjcrk{j}(:,2), seg_sjcrk{j}(:,3), seg_sjcrk{j}(:,4),'k*');
    hold on
    scatter3(ccoord(ladder_relative_cn(j,:)>0, 1), ccoord(ladder_relative_cn(j,:)>0,2), ccoord(ladder_relative_cn(j,:)>0, 3), 'g');
    hold on
    scatter3(ccoord(isnan(ladder_relative_cn(j,:)), 1), ccoord((isnan(ladder_relative_cn(j,:))),2), ccoord((isnan(ladder_relative_cn(j,:))), 3), 'y');
    hold on
    scatter3(ccoord(ladder_relative_cn(j,:)<0, 1), ccoord(ladder_relative_cn(j,:)<0,2), ccoord(ladder_relative_cn(j,:)<0, 3), 'r');
    %colormap(jet);
    %colorbar;
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
    %zlabel('R^2');
    %title(sprintf('Jump in z displacement during time window %d', j));
    hold off
end
%% Compare pre/postjump slope with jump 
ratio_prejump_bdz = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_prejump_bdz(:, j) = transpose( ladder_jump_bdz(j, :) ./ coef_prejump_bdz(j,:,1) ); %positive means jump in the same direction as prejump increase (or increase jump or decrease drop)
end

ratio_postjump_bdz = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_postjump_bdz(:, j) = transpose( ladder_jump_bdz(j, :) ./ coef_postjump_bdz(j,:,1) ); 
end

ratio_prejump_bdx = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_prejump_bdx(:, j) = transpose( ladder_jump_bdx(j, :) ./ coef_prejump_bdx(j,:,1) ); %positive means jump in the same direction as prejump increase (or increase jump or decrease drop)
end

ratio_postjump_bdx = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_postjump_bdx(:, j) = transpose( ladder_jump_bdx(j, :) ./ coef_postjump_bdx(j,:,1) ); 
end

ratio_prejump_bdy = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_prejump_bdy(:, j) = transpose( ladder_jump_bdy(j, :) ./ coef_prejump_bdy(j,:,1) ); %positive means jump in the same direction as prejump increase (or increase jump or decrease drop)
end

ratio_postjump_bdy = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_postjump_bdy(:, j) = transpose( ladder_jump_bdy(j, :) ./ coef_postjump_bdy(j,:,1) ); 
end

ratio_prejump_bdxoe = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_prejump_bdxoe(:, j) = transpose( ladder_jump_bdxoe(j, :) ./ coef_prejump_bdxoe(j,:,1) ); %positive means jump in the same direction as prejump increase (or increase jump or decrease drop)
end

ratio_postjump_bdxoe = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_postjump_bdxoe(:, j) = transpose( ladder_jump_bdxoe(j, :) ./ coef_postjump_bdxoe(j,:,1) ); 
end

ratio_prejump_bdyoe = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_prejump_bdyoe(:, j) = transpose( ladder_jump_bdyoe(j, :) ./ coef_prejump_bdyoe(j,:,1) ); %positive means jump in the same direction as prejump increase (or increase jump or decrease drop)
end

ratio_postjump_bdyoe = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_postjump_bdyoe(:, j) = transpose( ladder_jump_bdyoe(j, :) ./ coef_postjump_bdyoe(j,:,1) ); 
end

ratio_prejump_bdzoe = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_prejump_bdzoe(:, j) = transpose( ladder_jump_bdzoe(j, :) ./ coef_prejump_bdzoe(j,:,1) ); %positive means jump in the same direction as prejump increase (or increase jump or decrease drop)
end

ratio_postjump_bdzoe = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_postjump_bdzoe(:, j) = transpose( ladder_jump_bdzoe(j, :) ./ coef_postjump_bdzoe(j,:,1) ); 
end

ratio_prejump_cdx = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_prejump_cdx(:, j) = transpose( ladder_jump_cdx(j, :) ./ coef_prejump_cdx(j,:,1) ); %positive means jump in the same direction as prejump increase (or increase jump or decrease drop)
end

ratio_postjump_cdx = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_postjump_cdx(:, j) = transpose( ladder_jump_cdx(j, :) ./ coef_postjump_cdx(j,:,1) ); 
end

ratio_prejump_cdy = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_prejump_cdy(:, j) = transpose( ladder_jump_cdy(j, :) ./ coef_prejump_cdy(j,:,1) ); %positive means jump in the same direction as prejump increase (or increase jump or decrease drop)
end

ratio_postjump_cdy = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_postjump_cdy(:, j) = transpose( ladder_jump_cdy(j, :) ./ coef_postjump_cdy(j,:,1) ); 
end

ratio_prejump_cdz = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_prejump_cdz(:, j) = transpose( ladder_jump_cdz(j, :) ./ coef_prejump_cdz(j,:,1) ); %positive means jump in the same direction as prejump increase (or increase jump or decrease drop)
end

ratio_postjump_cdz = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_postjump_cdz(:, j) = transpose( ladder_jump_cdz(j, :) ./ coef_postjump_cdz(j,:,1) ); 
end

ratio_prejump_cs = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_prejump_cs(:, j) = transpose( ladder_jump_cs(j, :) ./ coef_prejump_cs(j,:,1) ); 
end

ratio_postjump_cs = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_postjump_cs(:, j) = transpose( ladder_jump_cs(j, :) ./ coef_postjump_cs(j,:,1) ); 
end

ratio_prejump_cn = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_prejump_cn(:, j) = transpose( ladder_jump_cn(j, :) ./ coef_prejump_cn(j,:,1) ); 
end

ratio_postjump_cn = zeros(num_usefulsensor, length(sbs)-1);
for j = 1:length(sbs)-1
    ratio_postjump_cn(:, j) = transpose( ladder_jump_cn(j, :) ./ coef_postjump_cn(j,:,1) ); 
end

% ratio_jump data list
%{
ratio_jump = zeros(num_usefulsensor, length(sbs)-1, 22);
ratio_jump(:, :, 1) = ratio_prejump_bdx;
ratio_jump(:, :, 2) = ratio_postjump_bdx;
ratio_jump(:, :, 3) = ratio_prejump_bdy;
ratio_jump(:, :, 4) = ratio_postjump_bdy;
ratio_jump(:, :, 5) = ratio_prejump_bdz;
ratio_jump(:, :, 6) = ratio_postjump_bdz;
ratio_jump(:, :, 7) = ratio_prejump_bdxoe;
ratio_jump(:, :, 8) = ratio_postjump_bdxoe;
ratio_jump(:, :, 9) = ratio_prejump_bdyoe;
ratio_jump(:, :, 10) = ratio_postjump_bdyoe;
ratio_jump(:, :, 11) = ratio_prejump_bdzoe;
ratio_jump(:, :, 12) = ratio_postjump_bdzoe;
ratio_jump(:, :, 13) = ratio_prejump_cdx;
ratio_jump(:, :, 14) = ratio_postjump_cdx;
ratio_jump(:, :, 15) = ratio_prejump_cdy;
ratio_jump(:, :, 16) = ratio_postjump_cdy;
ratio_jump(:, :, 17) = ratio_prejump_cdz;
ratio_jump(:, :, 18) = ratio_postjump_cdz;
ratio_jump(:, :, 19) = ratio_prejump_cn;
ratio_jump(:, :, 20) = ratio_postjump_cn;
ratio_jump(:, :, 21) = ratio_prejump_cs;
ratio_jump(:, :, 22) = ratio_postjump_cs;
%}
%{
ratio_jump.prejump_bdx = ratio_prejump_bdx;
ratio_jump.postjump_bdx = ratio_postjump_bdx;
ratio_jump.prejump_bdy = ratio_prejump_bdy;
ratio_jump.postjump_bdy = ratio_postjump_bdy;
ratio_jump.prejump_bdz = ratio_prejump_bdz;
ratio_jump.postjump_bdz = ratio_postjump_bdz;
ratio_jump.prejump_bdxoe = ratio_prejump_bdxoe;
ratio_jump.postjump_bdxoe = ratio_postjump_bdxoe;
ratio_jump.prejump_bdyoe = ratio_prejump_bdyoe;
ratio_jump.postjump_bdyoe = ratio_postjump_bdyoe;
ratio_jump.prejump_bdzoe = ratio_prejump_bdzoe;
ratio_jump.postjump_bdzoe = ratio_postjump_bdzoe;
ratio_jump.prejump_cdx = ratio_prejump_cdx;
ratio_jump.postjump_cdx = ratio_postjump_cdx;
ratio_jump.prejump_cdy = ratio_prejump_cdy;
ratio_jump.postjump_cdy = ratio_postjump_cdy;
ratio_jump.prejump_cdz = ratio_prejump_cdz;
ratio_jump.postjump_cdz = ratio_postjump_cdz;
ratio_jump.prejump_cn = ratio_prejump_cn;
ratio_jump.postjump_cn = ratio_postjump_cn;
ratio_jump.prejump_cs = ratio_prejump_cs;
ratio_jump.postjump_cs = ratio_postjump_cs;
%}
ratio_jump.prejump_bdx = ratio_prejump_bdx;
ratio_jump.prejump_bdxoe = ratio_prejump_bdxoe;
ratio_jump.prejump_bdy = ratio_prejump_bdy;
ratio_jump.prejump_bdyoe = ratio_prejump_bdyoe;
ratio_jump.prejump_bdz = ratio_prejump_bdz;
ratio_jump.prejump_bdzoe = ratio_prejump_bdzoe;
ratio_jump.prejump_cdx = ratio_prejump_cdx;
ratio_jump.prejump_cdy = ratio_prejump_cdy;
ratio_jump.prejump_cdz = ratio_prejump_cdz;
ratio_jump.prejump_cn = ratio_prejump_cn;
ratio_jump.prejump_cs = ratio_prejump_cs;
ratio_jump.postjump_bdx = ratio_postjump_bdx;
ratio_jump.postjump_bdxoe = ratio_postjump_bdxoe;
ratio_jump.postjump_bdy = ratio_postjump_bdy;
ratio_jump.postjump_bdyoe = ratio_postjump_bdyoe;
ratio_jump.postjump_bdz = ratio_postjump_bdz;
ratio_jump.postjump_bdzoe = ratio_postjump_bdzoe;
ratio_jump.postjump_cdx = ratio_postjump_cdx;
ratio_jump.postjump_cdy = ratio_postjump_cdy;
ratio_jump.postjump_cdz = ratio_postjump_cdz;
ratio_jump.postjump_cn = ratio_postjump_cn;
ratio_jump.postjump_cs = ratio_postjump_cs;

% coef_jump data list
%coef_jump = zeros(length(sbs)-1, num_usefulsensor, 3, 22);
%store coef_jump as an array such that it's easier to export and import
%{
coef_jump.prejump_bdx = coef_prejump_bdx;
coef_jump.postjump_bdx = coef_postjump_bdx;
coef_jump.prejump_bdy = coef_prejump_bdy;
coef_jump.postjump_bdy = coef_postjump_bdy;
coef_jump.prejump_bdz = coef_prejump_bdz;
coef_jump.postjump_bdz = coef_postjump_bdz;
coef_jump.prejump_bdxoe = coef_prejump_bdxoe;
coef_jump.postjump_bdxoe = coef_postjump_bdxoe;
coef_jump.prejump_bdyoe = coef_prejump_bdyoe;
coef_jump.postjump_bdyoe = coef_postjump_bdyoe;
coef_jump.prejump_bdzoe = coef_prejump_bdzoe;
coef_jump.postjump_bdzoe = coef_postjump_bdzoe;
coef_jump.prejump_cdx = coef_prejump_cdx;
coef_jump.postjump_cdx = coef_postjump_cdx;
coef_jump.prejump_cdy = coef_prejump_cdy;
coef_jump.postjump_cdy = coef_postjump_cdy;
coef_jump.prejump_cdz = coef_prejump_cdz;
coef_jump.postjump_cdz = coef_postjump_cdz;
coef_jump.prejump_cn = coef_prejump_cn;
coef_jump.postjump_cn = coef_postjump_cn;
coef_jump.prejump_cs = coef_prejump_cs;
coef_jump.postjump_cs = coef_postjump_cs;
%}
coef_jump.prejump_bdx = coef_prejump_bdx;
coef_jump.prejump_bdxoe = coef_prejump_bdxoe;
coef_jump.prejump_bdy = coef_prejump_bdy;
coef_jump.prejump_bdyoe = coef_prejump_bdyoe;
coef_jump.prejump_bdz = coef_prejump_bdz;
coef_jump.prejump_bdzoe = coef_prejump_bdzoe;
coef_jump.prejump_cdx = coef_prejump_cdx;
coef_jump.prejump_cdy = coef_prejump_cdy;
coef_jump.prejump_cdz = coef_prejump_cdz;
coef_jump.prejump_cn = coef_prejump_cn;
coef_jump.prejump_cs = coef_prejump_cs;
coef_jump.postjump_bdx = coef_postjump_bdx;
coef_jump.postjump_bdxoe = coef_postjump_bdxoe;
coef_jump.postjump_bdy = coef_postjump_bdy;
coef_jump.postjump_bdyoe = coef_postjump_bdyoe;
coef_jump.postjump_bdz = coef_postjump_bdz;
coef_jump.postjump_bdzoe = coef_postjump_bdzoe;
coef_jump.postjump_cdx = coef_postjump_cdx;
coef_jump.postjump_cdy = coef_postjump_cdy;
coef_jump.postjump_cdz = coef_postjump_cdz;
coef_jump.postjump_cn = coef_postjump_cn;
coef_jump.postjump_cs = coef_postjump_cs;

%
%categorize each contact based on ratio_prejump, ratio_postjump and
%coef_prejump. (Check notes)

for i = 1:length(sbs)-1
    for j = 1:num_usefulsensor
        if isnan(coef_postjump_bdx(i,j,1))
            jump_cat.bdx(j,i) = 9;
        elseif isnan(coef_prejump_bdx(i,j,1))
            jump_cat.bdx(j,i) = 10;
        elseif coef_jump.prejump_bdx(i,j,1) > 0 %category 1-4
            if ratio_jump.prejump_bdx(j,i) < 0 && ratio_jump.postjump_bdx(j,i) < 0
                jump_cat.bdx(j,i) = 1;
            elseif ratio_jump.prejump_bdx(j,i) > 0 && ratio_jump.postjump_bdx(j,i) > 0
                jump_cat.bdx(j,i) = 2; 
            elseif ratio_jump.prejump_bdx(j,i) > 0 && ratio_jump.postjump_bdx(j,i) < 0
                jump_cat.bdx(j,i) = 3; 
            elseif ratio_jump.prejump_bdx(j,i) < 0 && ratio_jump.postjump_bdx(j,i) > 0
                jump_cat.bdx(j,i) = 4;
            end
        elseif coef_jump.prejump_bdx(i,j,1) < 0 %category 5-8
            if ratio_jump.prejump_bdx(j,i) < 0 && ratio_jump.postjump_bdx(j,i) < 0
                jump_cat.bdx(j,i) = 5;
            elseif ratio_jump.prejump_bdx(j,i) > 0 && ratio_jump.postjump_bdx(j,i) > 0
                jump_cat.bdx(j,i) = 6; 
            elseif ratio_jump.prejump_bdx(j,i) > 0 && ratio_jump.postjump_bdx(j,i) < 0
                jump_cat.bdx(j,i) = 7; 
            elseif ratio_jump.prejump_bdx(j,i) < 0 && ratio_jump.postjump_bdx(j,i) > 0
                jump_cat.bdx(j,i) = 8;
            end
        end

        if isnan(coef_postjump_bdxoe(i,j,1))
            jump_cat.bdxoe(j,i) = 9;
        elseif isnan(coef_prejump_bdxoe(i,j,1))
            jump_cat.bdxoe(j,i) = 10;
        elseif coef_jump.prejump_bdxoe(i,j,1) > 0 %category 1-4
            if ratio_jump.prejump_bdxoe(j,i) < 0 && ratio_jump.postjump_bdxoe(j,i) < 0
                jump_cat.bdxoe(j,i) = 1;
            elseif ratio_jump.prejump_bdxoe(j,i) > 0 && ratio_jump.postjump_bdxoe(j,i) > 0
                jump_cat.bdxoe(j,i) = 2; 
            elseif ratio_jump.prejump_bdxoe(j,i) > 0 && ratio_jump.postjump_bdxoe(j,i) < 0
                jump_cat.bdxoe(j,i) = 3; 
            elseif ratio_jump.prejump_bdxoe(j,i) < 0 && ratio_jump.postjump_bdxoe(j,i) > 0
                jump_cat.bdxoe(j,i) = 4;
            end
        elseif coef_jump.prejump_bdxoe(i,j,1) < 0 %category 5-8
            if ratio_jump.prejump_bdxoe(j,i) < 0 && ratio_jump.postjump_bdxoe(j,i) < 0
                jump_cat.bdxoe(j,i) = 5;
            elseif ratio_jump.prejump_bdxoe(j,i) > 0 && ratio_jump.postjump_bdxoe(j,i) > 0
                jump_cat.bdxoe(j,i) = 6; 
            elseif ratio_jump.prejump_bdxoe(j,i) > 0 && ratio_jump.postjump_bdxoe(j,i) < 0
                jump_cat.bdxoe(j,i) = 7; 
            elseif ratio_jump.prejump_bdxoe(j,i) < 0 && ratio_jump.postjump_bdxoe(j,i) > 0
                jump_cat.bdxoe(j,i) = 8;
            end
        end
        
        if isnan(coef_postjump_bdy(i,j,1))
            jump_cat.bdy(j,i) = 9;
        elseif isnan(coef_prejump_bdy(i,j,1))
            jump_cat.bdy(j,i) = 10;
        elseif coef_jump.prejump_bdy(i,j,1) > 0 %category 1-4
            if ratio_jump.prejump_bdy(j,i) < 0 && ratio_jump.postjump_bdy(j,i) < 0
                jump_cat.bdy(j,i) = 1;
            elseif ratio_jump.prejump_bdy(j,i) > 0 && ratio_jump.postjump_bdy(j,i) > 0
                jump_cat.bdy(j,i) = 2; 
            elseif ratio_jump.prejump_bdy(j,i) > 0 && ratio_jump.postjump_bdy(j,i) < 0
                jump_cat.bdy(j,i) = 3; 
            elseif ratio_jump.prejump_bdy(j,i) < 0 && ratio_jump.postjump_bdy(j,i) > 0
                jump_cat.bdy(j,i) = 4;
            end
        elseif coef_jump.prejump_bdy(i,j,1) < 0 %category 5-8
            if ratio_jump.prejump_bdy(j,i) < 0 && ratio_jump.postjump_bdy(j,i) < 0
                jump_cat.bdy(j,i) = 5;
            elseif ratio_jump.prejump_bdy(j,i) > 0 && ratio_jump.postjump_bdy(j,i) > 0
                jump_cat.bdy(j,i) = 6; 
            elseif ratio_jump.prejump_bdy(j,i) > 0 && ratio_jump.postjump_bdy(j,i) < 0
                jump_cat.bdy(j,i) = 7; 
            elseif ratio_jump.prejump_bdy(j,i) < 0 && ratio_jump.postjump_bdy(j,i) > 0
                jump_cat.bdy(j,i) = 8;
            end
        end
        
        if isnan(coef_postjump_bdyoe(i,j,1)) %dead contact after jump
            jump_cat.bdyoe(j,i) = 9;
        elseif isnan(coef_prejump_bdyoe(i,j,1)) %already dead contact before jump
            jump_cat.bdyoe(j,i) = 10;
        elseif coef_jump.prejump_bdyoe(i,j,1) > 0 %category 1-4
            if ratio_jump.prejump_bdyoe(j,i) < 0 && ratio_jump.postjump_bdyoe(j,i) < 0
                jump_cat.bdyoe(j,i) = 1;
            elseif ratio_jump.prejump_bdyoe(j,i) > 0 && ratio_jump.postjump_bdyoe(j,i) > 0
                jump_cat.bdyoe(j,i) = 2; 
            elseif ratio_jump.prejump_bdyoe(j,i) > 0 && ratio_jump.postjump_bdyoe(j,i) < 0
                jump_cat.bdyoe(j,i) = 3; 
            elseif ratio_jump.prejump_bdyoe(j,i) < 0 && ratio_jump.postjump_bdyoe(j,i) > 0
                jump_cat.bdyoe(j,i) = 4;
            end
        elseif coef_jump.prejump_bdyoe(i,j,1) < 0 %category 5-8
            if ratio_jump.prejump_bdyoe(j,i) < 0 && ratio_jump.postjump_bdyoe(j,i) < 0
                jump_cat.bdyoe(j,i) = 5;
            elseif ratio_jump.prejump_bdyoe(j,i) > 0 && ratio_jump.postjump_bdyoe(j,i) > 0
                jump_cat.bdyoe(j,i) = 6; 
            elseif ratio_jump.prejump_bdyoe(j,i) > 0 && ratio_jump.postjump_bdyoe(j,i) < 0
                jump_cat.bdyoe(j,i) = 7; 
            elseif ratio_jump.prejump_bdyoe(j,i) < 0 && ratio_jump.postjump_bdyoe(j,i) > 0
                jump_cat.bdyoe(j,i) = 8;
            end
        end
        
        
        if isnan(coef_postjump_bdz(i,j,1))
            jump_cat.bdz(j,i) = 9;
        elseif isnan(coef_prejump_bdz(i,j,1))
            jump_cat.bdz(j,i) = 10;
        elseif coef_jump.prejump_bdz(i,j,1) > 0 %category 1-4
            if ratio_jump.prejump_bdz(j,i) < 0 && ratio_jump.postjump_bdz(j,i) < 0
                jump_cat.bdz(j,i) = 1;
            elseif ratio_jump.prejump_bdz(j,i) > 0 && ratio_jump.postjump_bdz(j,i) > 0
                jump_cat.bdz(j,i) = 2; 
            elseif ratio_jump.prejump_bdz(j,i) > 0 && ratio_jump.postjump_bdz(j,i) < 0
                jump_cat.bdz(j,i) = 3; 
            elseif ratio_jump.prejump_bdz(j,i) < 0 && ratio_jump.postjump_bdz(j,i) > 0
                jump_cat.bdz(j,i) = 4;
            end
        elseif coef_jump.prejump_bdz(i,j,1) < 0 %category 5-8
            if ratio_jump.prejump_bdz(j,i) < 0 && ratio_jump.postjump_bdz(j,i) < 0
                jump_cat.bdz(j,i) = 5;
            elseif ratio_jump.prejump_bdz(j,i) > 0 && ratio_jump.postjump_bdz(j,i) > 0
                jump_cat.bdz(j,i) = 6; 
            elseif ratio_jump.prejump_bdz(j,i) > 0 && ratio_jump.postjump_bdz(j,i) < 0
                jump_cat.bdz(j,i) = 7; 
            elseif ratio_jump.prejump_bdz(j,i) < 0 && ratio_jump.postjump_bdz(j,i) > 0
                jump_cat.bdz(j,i) = 8;
            end
        end
        
        if isnan(coef_postjump_bdzoe(i,j,1)) %dead contact after jump
            jump_cat.bdzoe(j,i) = 9;
        elseif isnan(coef_prejump_bdzoe(i,j,1)) %already dead contact before jump
            jump_cat.bdzoe(j,i) = 10;
        elseif coef_jump.prejump_bdzoe(i,j,1) > 0 %category 1-4
            if ratio_jump.prejump_bdzoe(j,i) < 0 && ratio_jump.postjump_bdzoe(j,i) < 0
                jump_cat.bdzoe(j,i) = 1;
            elseif ratio_jump.prejump_bdzoe(j,i) > 0 && ratio_jump.postjump_bdzoe(j,i) > 0
                jump_cat.bdzoe(j,i) = 2; 
            elseif ratio_jump.prejump_bdzoe(j,i) > 0 && ratio_jump.postjump_bdzoe(j,i) < 0
                jump_cat.bdzoe(j,i) = 3; 
            elseif ratio_jump.prejump_bdzoe(j,i) < 0 && ratio_jump.postjump_bdzoe(j,i) > 0
                jump_cat.bdzoe(j,i) = 4;
            end
        elseif coef_jump.prejump_bdzoe(i,j,1) < 0 %category 5-8
            if ratio_jump.prejump_bdzoe(j,i) < 0 && ratio_jump.postjump_bdzoe(j,i) < 0
                jump_cat.bdzoe(j,i) = 5;
            elseif ratio_jump.prejump_bdzoe(j,i) > 0 && ratio_jump.postjump_bdzoe(j,i) > 0
                jump_cat.bdzoe(j,i) = 6; 
            elseif ratio_jump.prejump_bdzoe(j,i) > 0 && ratio_jump.postjump_bdzoe(j,i) < 0
                jump_cat.bdzoe(j,i) = 7; 
            elseif ratio_jump.prejump_bdzoe(j,i) < 0 && ratio_jump.postjump_bdzoe(j,i) > 0
                jump_cat.bdzoe(j,i) = 8;
            end
        end
        
        
        if isnan(coef_postjump_cdz(i,j,1))
            jump_cat.cdz(j,i) = 9;
        elseif isnan(coef_prejump_cdz(i,j,1))
            jump_cat.cdz(j,i) = 10;
        elseif coef_jump.prejump_cdz(i,j,1) > 0 %category 1-4
            if ratio_jump.prejump_cdz(j,i) < 0 && ratio_jump.postjump_cdz(j,i) < 0
                jump_cat.cdz(j,i) = 1;
            elseif ratio_jump.prejump_cdz(j,i) > 0 && ratio_jump.postjump_cdz(j,i) > 0
                jump_cat.cdz(j,i) = 2; 
            elseif ratio_jump.prejump_cdz(j,i) > 0 && ratio_jump.postjump_cdz(j,i) < 0
                jump_cat.cdz(j,i) = 3; 
            elseif ratio_jump.prejump_cdz(j,i) < 0 && ratio_jump.postjump_cdz(j,i) > 0
                jump_cat.cdz(j,i) = 4;
            end
        elseif coef_jump.prejump_cdz(i,j,1) < 0 %category 5-8
            if ratio_jump.prejump_cdz(j,i) < 0 && ratio_jump.postjump_cdz(j,i) < 0
                jump_cat.cdz(j,i) = 5;
            elseif ratio_jump.prejump_cdz(j,i) > 0 && ratio_jump.postjump_cdz(j,i) > 0
                jump_cat.cdz(j,i) = 6; 
            elseif ratio_jump.prejump_cdz(j,i) > 0 && ratio_jump.postjump_cdz(j,i) < 0
                jump_cat.cdz(j,i) = 7; 
            elseif ratio_jump.prejump_cdz(j,i) < 0 && ratio_jump.postjump_cdz(j,i) > 0
                jump_cat.cdz(j,i) = 8;
            end
        end
        
        if isnan(coef_postjump_cdy(i,j,1))
            jump_cat.cdy(j,i) = 9;
        elseif isnan(coef_prejump_cdy(i,j,1))
            jump_cat.cdy(j,i) = 10;
        elseif coef_jump.prejump_cdy(i,j,1) > 0 %category 1-4
            if ratio_jump.prejump_cdy(j,i) < 0 && ratio_jump.postjump_cdy(j,i) < 0
                jump_cat.cdy(j,i) = 1;
            elseif ratio_jump.prejump_cdy(j,i) > 0 && ratio_jump.postjump_cdy(j,i) > 0
                jump_cat.cdy(j,i) = 2; 
            elseif ratio_jump.prejump_cdy(j,i) > 0 && ratio_jump.postjump_cdy(j,i) < 0
                jump_cat.cdy(j,i) = 3; 
            elseif ratio_jump.prejump_cdy(j,i) < 0 && ratio_jump.postjump_cdy(j,i) > 0
                jump_cat.cdy(j,i) = 4;
            end
        elseif coef_jump.prejump_cdy(i,j,1) < 0 %category 5-8
            if ratio_jump.prejump_cdy(j,i) < 0 && ratio_jump.postjump_cdy(j,i) < 0
                jump_cat.cdy(j,i) = 5;
            elseif ratio_jump.prejump_cdy(j,i) > 0 && ratio_jump.postjump_cdy(j,i) > 0
                jump_cat.cdy(j,i) = 6; 
            elseif ratio_jump.prejump_cdy(j,i) > 0 && ratio_jump.postjump_cdy(j,i) < 0
                jump_cat.cdy(j,i) = 7; 
            elseif ratio_jump.prejump_cdy(j,i) < 0 && ratio_jump.postjump_cdy(j,i) > 0
                jump_cat.cdy(j,i) = 8;
            end
        end
        
        if isnan(coef_postjump_cdx(i,j,1))
            jump_cat.cdx(j,i) = 9;
        elseif isnan(coef_prejump_cdx(i,j,1))
            jump_cat.cdx(j,i) = 10;
        elseif coef_jump.prejump_cdx(i,j,1) > 0 %category 1-4
            if ratio_jump.prejump_cdx(j,i) < 0 && ratio_jump.postjump_cdx(j,i) < 0
                jump_cat.cdx(j,i) = 1;
            elseif ratio_jump.prejump_cdx(j,i) > 0 && ratio_jump.postjump_cdx(j,i) > 0
                jump_cat.cdx(j,i) = 2; 
            elseif ratio_jump.prejump_cdx(j,i) > 0 && ratio_jump.postjump_cdx(j,i) < 0
                jump_cat.cdx(j,i) = 3; 
            elseif ratio_jump.prejump_cdx(j,i) < 0 && ratio_jump.postjump_cdx(j,i) > 0
                jump_cat.cdx(j,i) = 4;
            end
        elseif coef_jump.prejump_cdx(i,j,1) < 0 %category 5-8
            if ratio_jump.prejump_cdx(j,i) < 0 && ratio_jump.postjump_cdx(j,i) < 0
                jump_cat.cdx(j,i) = 5;
            elseif ratio_jump.prejump_cdx(j,i) > 0 && ratio_jump.postjump_cdx(j,i) > 0
                jump_cat.cdx(j,i) = 6; 
            elseif ratio_jump.prejump_cdx(j,i) > 0 && ratio_jump.postjump_cdx(j,i) < 0
                jump_cat.cdx(j,i) = 7; 
            elseif ratio_jump.prejump_cdx(j,i) < 0 && ratio_jump.postjump_cdx(j,i) > 0
                jump_cat.cdx(j,i) = 8;
            end
        end
        
        if isnan(coef_postjump_cn(i,j,1))
            jump_cat.cn(j,i) = 9;
        elseif isnan(coef_prejump_cn(i,j,1))
            jump_cat.cn(j,i) = 10;
        elseif coef_jump.prejump_cn(i,j,1) > 0 %category 1-4
            if ratio_jump.prejump_cn(j,i) < 0 && ratio_jump.postjump_cn(j,i) < 0
                jump_cat.cn(j,i) = 1;
            elseif ratio_jump.prejump_cn(j,i) > 0 && ratio_jump.postjump_cn(j,i) > 0
                jump_cat.cn(j,i) = 2; 
            elseif ratio_jump.prejump_cn(j,i) > 0 && ratio_jump.postjump_cn(j,i) < 0
                jump_cat.cn(j,i) = 3; 
            elseif ratio_jump.prejump_cn(j,i) < 0 && ratio_jump.postjump_cn(j,i) > 0
                jump_cat.cn(j,i) = 4;
            end
        elseif coef_jump.prejump_cn(i,j,1) < 0 %category 5-8
            if ratio_jump.prejump_cn(j,i) < 0 && ratio_jump.postjump_cn(j,i) < 0
                jump_cat.cn(j,i) = 5;
            elseif ratio_jump.prejump_cn(j,i) > 0 && ratio_jump.postjump_cn(j,i) > 0
                jump_cat.cn(j,i) = 6; 
            elseif ratio_jump.prejump_cn(j,i) > 0 && ratio_jump.postjump_cn(j,i) < 0
                jump_cat.cn(j,i) = 7; 
            elseif ratio_jump.prejump_cn(j,i) < 0 && ratio_jump.postjump_cn(j,i) > 0
                jump_cat.cn(j,i) = 8;
            end
        end
        
        if isnan(coef_postjump_cs(i,j,1))
            jump_cat.cs(j,i) = 9;
        elseif isnan(coef_prejump_cs(i,j,1))
            jump_cat.cs(j,i) = 10;
        elseif coef_jump.prejump_cs(i,j,1) > 0 %category 1-4
            if ratio_jump.prejump_cs(j,i) < 0 && ratio_jump.postjump_cs(j,i) < 0
                jump_cat.cs(j,i) = 1;
            elseif ratio_jump.prejump_cs(j,i) > 0 && ratio_jump.postjump_cs(j,i) > 0
                jump_cat.cs(j,i) = 2; 
            elseif ratio_jump.prejump_cs(j,i) > 0 && ratio_jump.postjump_cs(j,i) < 0
                jump_cat.cs(j,i) = 3; 
            elseif ratio_jump.prejump_cs(j,i) < 0 && ratio_jump.postjump_cs(j,i) > 0
                jump_cat.cs(j,i) = 4;
            end
        elseif coef_jump.prejump_cs(i,j,1) < 0 %category 5-8
            if ratio_jump.prejump_cs(j,i) < 0 && ratio_jump.postjump_cs(j,i) < 0
                jump_cat.cs(j,i) = 5;
            elseif ratio_jump.prejump_cs(j,i) > 0 && ratio_jump.postjump_cs(j,i) > 0
                jump_cat.cs(j,i) = 6; 
            elseif ratio_jump.prejump_cs(j,i) > 0 && ratio_jump.postjump_cs(j,i) < 0
                jump_cat.cs(j,i) = 7; 
            elseif ratio_jump.prejump_cs(j,i) < 0 && ratio_jump.postjump_cs(j,i) > 0
                jump_cat.cs(j,i) = 8;
            end
        end
        
    end
end

% ladder jump for export (stacking purpose)
ladder_jump.bdx = ladder_jump_bdx;
ladder_jump.bdxoe = ladder_jump_bdxoe;
ladder_jump.bdy = ladder_jump_bdy;
ladder_jump.bdyoe = ladder_jump_bdyoe;
ladder_jump.bdz = ladder_jump_bdz;
ladder_jump.bdzoe = ladder_jump_bdzoe;
ladder_jump.cdx = ladder_jump_cdx;
ladder_jump.cdy = ladder_jump_cdy;
ladder_jump.cdz = ladder_jump_cdz;
ladder_jump.cn = ladder_jump_cn;
ladder_jump.cs = ladder_jump_cs;
%%
figure(f); f=f+1;
sp = ceil(sqrt(length(sbs)-1));
for j = 1:length(sbs)-1
    subplot(sp, sp, j);
    %increase and jump
    scatter3(ccoord( ((ratio_prejump_bdz(:, j)> 0)&ladder_jump_bdz(j, :)'>0), 1), ccoord(((ratio_prejump_bdz(:, j)>0)&ladder_jump_bdz(j, :)'>0),2), ccoord(((ratio_prejump_bdz(:, j)> 0)&ladder_jump_bdz(j, :)'>0), 3) , 'b');
    hold on
    %decrease and drop
    scatter3(ccoord( (ratio_prejump_bdz(:, j)<=0&ladder_jump_bdz(j, :)'>0), 1), ccoord((ratio_prejump_bdz(:, j)<= 0&ladder_jump_bdz(j, :)'>0),2), ccoord((ratio_prejump_bdz(:, j)<= 0&ladder_jump_bdz(j, :)'>0), 3) , 'r');
    hold on
    %decrease and jump
    scatter3(ccoord( (ratio_prejump_bdz(:, j)<=0&ladder_jump_bdz(j, :)'<0), 1), ccoord((ratio_prejump_bdz(:, j)<= 0&ladder_jump_bdz(j, :)'<0),2), ccoord((ratio_prejump_bdz(:, j)<= 0&ladder_jump_bdz(j, :)'<0), 3) , 'r');
    hold on
    %increase and jump   or    decrease and drop
    scatter3(ccoord( (ratio_prejump_cs(:, j)>0), 1), ccoord((ratio_prejump_cs(:, j)> 0),2), ccoord((ratio_prejump_cs(:, j)> 0), 3) , '*g');
    hold on
    % increase and drop  or   decrease and jump
    scatter3(ccoord( (ratio_prejump_cs(:, j)<= 0), 1), ccoord((ratio_prejump_cs(:, j)<= 0),2), ccoord((ratio_prejump_cs(:, j)<= 0), 3) , '*y');
    hold off
    view([0, 90]);
    xlabel('x(mm)');
    ylabel('y(mm)');
end



%% b value of unclustered AE
temp_n = 20;
Mw_un_range = linspace(min(Mw_un), max(Mw_un), temp_n);
for i = 1:temp_n
    mag_occur(i) = length(find(Mw_un>=Mw_un_range(i)));
end
mag_range_linfit = Mw_un_range( 1:min(find(mag_occur == 1)) );
mag_occur_linfit = Mw_un_range( 1:min(find(mag_occur == 1)) );
[b_coef(length(sbs),:), b_linfit{length(sbs)}] = onestep_linfit(mag_range_linfit, log10(mag_occur_linfit));
[b_coef(1,:), b_linfit{1}] = onestep_linfit( Mw_un_range, log10(mag_occur));
figure(f); f=f+1;
scatter(Mw_un_range, log10(mag_occur));
%% moment radius
figure(f); f= f+1;
scatter(ms_sj(:, 7), log10(ms_sj(:, 8)));
%% bvz
bvz_average = mean(bvz, 2);
figure(f); f = f+1;
plot(bvz_average);
%%
bdz_mean = mean(bdz, 1);
bdz_average_pos = mean(bdz(:,bdz_mean>0), 2);
bdz_average_neg = mean(bdz(:,bdz_mean<0), 2);
plot(bdz_average_pos);
hold on
plot(bdz_average_neg);
hold off

%% nominal friction coefficient
cn_mean = mean(cn, 2);
cs_mean = mean(cs, 2);
%mu = normalize(cs_mean) ./ normalize(cn_mean);
mu = cs_mean ./ cn_mean;
ss = csvread('S201-longer-ss.csv',1);
figure(f); f= f+1;
plot(ss(2:end, 1), mu);
figure(f); f= f+1;
plot(ss(2:end, 1), cn_mean);
figure(f); f= f+1;
plot(ss(2:end, 1), cs_mean);

%% Mainshock initiation site
temp = seg_ms_sj{window_stickslip};
sj_stickslip_ini = zeros(size(temp,1), size(seg_sjcrk{window_stickslip}, 2));
for i = 1:size(temp,1)
    tempi = find(seg_sjcrk{window_stickslip}(:,1) == temp(i,2) );
    if isempty(tempi) == 1
        tempi = find((seg_sjcrk{window_stickslip}(:,1) - temp(i,2))>=5 , 1 ); % 1 stands for finding the first
    end
    
    if length(tempi) > 1 %if there are more than one sjcrk happening coincident with an ms_sj, find the sjcrk that locates closest to the centre of the corresponding ms_sj, and identify that particular sjcrk as the initation site
        tempdist = zeros(length(tempi), 1);
        for k = 1:length(tempi)
            temp1 = [seg_sjcrk{window_stickslip}(tempi(k), 2), seg_sjcrk{window_stickslip}(tempi(k), 3), seg_sjcrk{window_stickslip}(tempi(k), 4) ] ;  %coordinate of each sjcrk
            temp2 = [seg_ms_sj{window_stickslip}(i, 4), seg_ms_sj{window_stickslip}(i, 5), seg_ms_sj{window_stickslip}(i, 6) ]; %coordinate of the i-th ms_sj
            tempdist(k) = norm(temp1 - temp2); %distance between each sjcrk and ms_sj
        end
        tempi2 = tempi( find(tempdist == min(tempdist)) ); %the index of tempi corresponding to the particular sjcrk that is closest to the ms_sj
        clear tempi
        tempi = tempi2;
        clear tempi2
    end
    
    sj_stickslip_ini(i, :) = seg_sjcrk{window_stickslip}(tempi, :);
end

figure1 = figure(f); f= f+1;
temp = sj_stickslip_ini(:,2);
a = [1:length(temp)]'; b = num2str(a); c=cellstr(b);
scatter3(sj_stickslip_ini(:,2), sj_stickslip_ini(:,3), sj_stickslip_ini(:,4));
text(sj_stickslip_ini(:,2), sj_stickslip_ini(:,3), sj_stickslip_ini(:,4), c);
view([0, 90]);
xlim([-20e-3, 20e-3]);
ylim([-20e-3, 20e-3]);
title('Initiation sites during stickslip');
figurename = strcat(fnumber, 'sjini_stickslip.png');
saveas(figure1, figurename);

ss_stressdrop = [ seg_strain{window_stickslip}(ind_stress_ws_max:ind_stress_ws_min), ...
    seg_stress{window_stickslip}(ind_stress_ws_max:ind_stress_ws_min) ] ;
figure1 = figure(f); f=f+1;
plot(ss_stressdrop(:,1), ss_stressdrop(:,2));
%plot(seg_strain{window_stickslip}(ind_stress_ws_max:ind_stress_ws_min), seg_stress{window_stickslip}(ind_stress_ws_max:ind_stress_ws_min));
title(strcat('stress vs strain during stickslip for ',fnumber ));
figurename = strcat(fnumber,'stressdrop_stickslip.png');
saveas(figure1, figurename);




%% data exporting
save(strcat('S',fnumber,'ladder_jump'),'ladder_jump');
save(strcat('S',fnumber,'jump_cat'),'jump_cat');
save(strcat('S',fnumber,'coef_jump'),'coef_jump');
save(strcat('S',fnumber,'ratio_jump'),'ratio_jump');
save(strcat('S',fnumber,'sbs'),'sbs');
save(strcat('S',fnumber,'bcoord'),'bcoord');
save(strcat('S',fnumber,'seg_sjcrk'),'seg_sjcrk');
save(strcat('S',fnumber,'sj_ss_ini'),'sj_stickslip_ini');
save(strcat('S',fnumber,'sj_kernel2000'),'ySix');
save(strcat('S',fnumber,'kernelpks_leverarm'),'kernelpks_leverarm');
save(strcat('S',fnumber,'kernelpks'),'kernelpeaks');
save(strcat('S',fnumber,'sj_spatdist'),'sj_spatdist');
save(strcat('S',fnumber,'pf_ss_sj'),'pf_ss_sj');
save(strcat('S',fnumber,'ss_stressdrop'), 'ss_stressdrop');
save(strcat('S',fnumber,'fracDcoef_sj'),'fracDcoef_sj');
save(strcat('S',fnumber,'fracDcoef_ms'),'fracDcoef_ms');
save(strcat('S',fnumber,'b_coef'),'b_coef');
save(strcat('S',fnumber,'cyc_kernelpeakslocs'),'cyc_kernelpeakslocs');
save(strcat('S',fnumber,'seg_ms_sj'),'seg_ms_sj');
save(strcat('S',fnumber,'broken_trigger'),'broken_trigger');
save(strcat('S',fnumber,'broken_inttriggerx'),'broken_inttriggerx');
save(strcat('S',fnumber,'broken_inttriggery'),'broken_inttriggery');
save(strcat('S',fnumber,'broken_inttriggerz'),'broken_inttriggerz');
save(strcat('S',fnumber,'broken_stresstrigger'),'broken_stresstrigger');
save(strcat('S',fnumber,'seg_cn_prejump'),'seg_cn_prejump');
save(strcat('S',fnumber,'seg_cs_prejump'),'seg_cs_prejump');
display('save data complete!');

%%
temp = seg_sjcrk{i}(:,1);
for i = 1:28
    temp = [temp; seg_sjcrk{i}(:,i)];
end
%{
ub = floor(temp/temp_n).*(1:1:temp_n);
lb = floor(temp/temp_n).*(1:1:temp_n) - floor(temp/temp_n) + 1;
stickslip_occur = zeros(floor(temp/temp_n), temp_n);
stickslip_x = zeros(floor(temp/temp_n), temp_n);
stickslip_y = zeros(floor(temp/temp_n), temp_n);
stickslip_z = zeros(floor(temp/temp_n), temp_n);
for i = 1:temp_n
    stickslip_occur(:, i) = seg_sjcrk{window_stickslip}(lb(i):ub(i), 1);
    stickslip_x(:, i) = seg_sjcrk{window_stickslip}(lb(i):ub(i), 2);
    stickslip_y(:, i) = seg_sjcrk{window_stickslip}(lb(i):ub(i), 3);
    stickslip_z(:, i) = seg_sjcrk{window_stickslip}(lb(i):ub(i), 4);
end

figure(f); f=f+1;
sp1 = ceil(sqrt(temp_n));
sp2 = ceil(sqrt(temp_n));
for i = 1:temp_n
    ah = subplot(sp1, sp2, i);
    %scatter3(stickslip_x(:,i), stickslip_y(:,i), stickslip_z(:,i));
    scatter3(reshape(stickslip_x(:,1:i), [], 1), reshape(stickslip_y(:,1:i), [], 1), reshape(stickslip_z(:,1:i), [], 1));
    view([0, 90]);
    xlim([-20e-3, 20e-3]);
    ylim([-20e-3, 20e-3]);
    %plot(stickslip_occur(:,i), stickslip_z(:,i));
end
annotation('textbox', [.9 0.5 0.1 0.2], 'String',sprintf('Crack propagation during stick slip, splitted into 25 sub-episodes, projected onto x-y plane.'), 'EdgeColor', 'none');
%}

%foreshock magnitudes, closeness to mainshock in time, location
%mainshock stress drop, rate of stress drop, propagation direction
