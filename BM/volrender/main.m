fname = '/Users/skyeong/connectome/atlas/atlas_n116.nii';
idvis = 3:5; % index to visualize
col = hot(length(idvis)); % color to map


figure;
visualize_results(fname,idvis,col,3); view(0,90);
camlight(60,-30);
camlight(-60,30);
% print('-dpng','-r300','iaal_n100.png')
