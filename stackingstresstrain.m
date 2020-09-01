f = 1;
%farray = {'201'};
% farray = {'201','208','214','215','216','220','221','222','223','224','225', '226','227','228'...
%     '229','230','231','232','233','234','235', '236','237','238','239','240','241'};
farray = {'413'};
% farray = {'487','488','489','490','486','485','484','476','475','474','473','472',...
%     '464','483','471','449','470','455','469','458','468','467','479','482','481','480','463'};
figure1 = figure(f);f=f+1;
for i = 1:length(farray)
    ss = csvread(strcat('S',farray{i},'-ss.csv'), 1);
    ss = ss(1:(end-1), :);
    plot(ss(:,1), ss(:,2));
    xlabel('strain'); ylabel('deviatoric stress (MPa)'); title('Stress-strain');
    saveas(figure1, strcat('stressstrain',farray{i},'.png'));
    %hold on
end
% hold off
% figure(f);f=f+1;
% for i = 1:length(farray)
%     plot(SS{i}(:,1), SS{i}(:,2));hold on
% end
%hold off; 