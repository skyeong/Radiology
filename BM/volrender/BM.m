% idvis = 0:30; % index to visualize
% col = hot(length(idvis)); % color to map
idvis=1:116;
col = jet(length(idvis)); % color to map

% Five colors
% col=[234,84,85;  102,167,225;  244,223,70;  69,201,102;  248,133,254; ]/255;



figure;
fname = 'AAL_L_orig.nii';  % resample in 2mm
visualize_rois(fname,idvis,col,0); view(90,0)
camlight(10,-10);
print('-dpng','-r300','AAL_L_medial.png')


figure;
fname = 'AAL_L_orig.nii';
visualize_rois(fname,idvis,col,0); view(-90,0)
camlight(0,10)
print('-dpng','-r300','AAL_L_lateral.png')


figure;
fname = 'AAL_R_orig.nii';
visualize_rois(fname,idvis,col,0); view(-90,0);
camlight(-10,10);
print('-dpng','-r300','AAL_R_medial.png')


figure;
fname = 'AAL_R_orig.nii';
visualize_rois(fname,idvis,col,0); view(90,0);
camlight(0,-10);
print('-dpng','-r300','AAL_R_laterial.png')


figure;
fname = 'AAL_orig.nii';
visualize_rois(fname,idvis,col,0); view(0,90);
camlight(60,-30);
camlight(-60,30);
print('-dpng','-r300','AAL_top.png')
