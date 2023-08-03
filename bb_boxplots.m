function bb_boxplots(group_num, vals, cont_num, asd_num, adhd_num, asdadhd_num, asd_code, asdadhd_code)

hb = boxplot(vals, group_num, 'PlotStyle', 'compact', 'medianstyle', 'line', 'jitter', 0);
hOutliers = findobj(hb,'Tag','Outliers');
yy = get(hOutliers,'YData');

out_val = [cell2mat(yy(1)),cell2mat(yy(2)),cell2mat(yy(3)),cell2mat(yy(4))];
out_val(isnan(out_val)) = [];

out_idx = zeros(2, length(out_val));
out_idx(2, :) = 1;
for i = 1:length(out_val)
    out_idx(1,i) = find(vals == out_val(i));
    
    % Assign group to outlier
    if out_idx(1,i)>(length(cont_num)+length(asd_num)+length(adhd_num))
        out_idx(1,i) = out_idx(1,i)-(length(cont_num)+length(asd_num)+length(adhd_num));
        out_idx(2,i) = 4;
    elseif out_idx(1,i)>(length(cont_num)+length(asd_num)) && out_idx(1,i)<=(length(cont_num)+length(asd_num)+length(adhd_num))
        out_idx(1,i) = out_idx(1,i)-(length(cont_num)+length(asd_num));
        out_idx(2,i) = 3;
    elseif out_idx(1,i)>(length(cont_num)) && out_idx(1,i)<=(length(cont_num)+length(asd_num))
        out_idx(1,i) = out_idx(1,i)-length(cont_num);
        out_idx(2,i) = 2;
    end
    
end

out1 = out_idx(1, find(out_idx(2,:) ==1));
out2 = out_idx(1, find(out_idx(2,:) ==2));
out3 = out_idx(1, find(out_idx(2,:) ==3));
out4 = out_idx(1, find(out_idx(2,:) ==4));

a = get(get(gca, 'children'), 'children');
box_asdadhd = a(1:4:20); set(box_asdadhd, 'color',[0.9 0.9 0.9])
box_adhd = a(2:4:20); set(box_adhd, 'color',[0.65 0.65 0.65])
box_asd = a(3:4:20); set(box_asd, 'color',[0.35 0.35 0.35])
box_cont = a(4:4:20); set(box_cont, 'color',[0.1 0.1 0.1])

set(findobj(gcf, 'tag', 'Outliers'),'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 0.1);
set(findobj(gcf, 'tag', 'Median'),'color','r', 'LineWidth', 3);
set(findobj(gcf, 'tag', 'Whisker'), 'LineWidth', 3);
set(findobj(gcf, 'tag', 'Box'), 'LineWidth', 40);

help_vect_cont = ones(1,length(cont_num))+randn(1,length(cont_num))/20;
help_vect_asd = ones(1,length(asd_num))+randn(1,length(asd_num))/40;
help_vect_adhd = ones(1,length(adhd_num))+randn(1,length(adhd_num))/60;
help_vect_asdadhd = ones(1,length(asdadhd_num))+randn(1,length(asdadhd_num))/80;
hold on
scatter(help_vect_cont(out1),vals(out1), 180, 'Marker', '+', 'MarkerEdgeColor', 'r', 'Linewidth', 2);
scatter(2*help_vect_asd(out2),vals(out2+length(cont_num)), 180, 'Marker', '+', 'MarkerEdgeColor', 'r', 'Linewidth', 2);
scatter(3*help_vect_adhd(out3),vals(out3+length(cont_num)+length(asd_num)), 180, 'Marker', '+', 'MarkerEdgeColor', 'r', 'Linewidth', 2);
scatter(4*help_vect_asdadhd(out4),vals(out4+length(cont_num)+length(asd_num)+length(adhd_num)), 180, 'Marker', '+', 'MarkerEdgeColor', 'r', 'Linewidth', 2);

f1=scatter(help_vect_cont,vals(group_num==1),110, 'k','filled');f1.MarkerFaceAlpha = 0.5;
f2=scatter(2*help_vect_asd,vals(group_num==2),110,'k','filled');f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f3=scatter(3*help_vect_adhd,vals(group_num==3),110,'k','filled');f3.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f4=scatter(4*help_vect_asdadhd,vals(group_num==4),110,'k','filled');f4.MarkerFaceAlpha = f1.MarkerFaceAlpha;

% Mark infants with ASD diag
f5=scatter(2*help_vect_asd(asd_code),vals(asd_code+length(cont_num)),120,'Marker', 'd','MarkerFaceColor','g', 'MarkerEdgeColor','g');f5.MarkerFaceAlpha = 1; 
f6=scatter(4*help_vect_asdadhd(asdadhd_code),vals(asdadhd_code+length(cont_num)+length(asd_num)+length(adhd_num)),120,'Marker', 'd','MarkerFaceColor','g', 'MarkerEdgeColor','g');f6.MarkerFaceAlpha = 1;



xlim([0.5, 4.5])
set(gca,'Xtick', 1:4, 'XTickLabels', {'TL'; 'ASD'; 'ADHD'; 'ASD-ADHD'}, 'FontSize', 32)
box off