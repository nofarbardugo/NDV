function [] = coverage_fig( cloud, consensus, name )
% plots coverage, mutations and gaps

top = max(sum(cloud));

gap_pos = find(consensus=='N');

figure
hold on

%plot coverage
plot(sum(cloud),'k');

%plot number of mutation from consensus
mut_cloud = get_mut_cloud(cloud);
plot(sum(mut_cloud),'b');

%plot gaps
plot(gap_pos,zeros(1,length(gap_pos))+top/2,'*r');

leg = cell(1,3);
leg{1} = 'coverage';
leg{2} = '# mutations';
leg{3} = 'gaps';
legend(leg)

title(sprintf('%s coverage',name),'fontsize',14);
xlabel('position','fontsize',13);

end

