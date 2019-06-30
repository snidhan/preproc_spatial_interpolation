%% Written by Sheel Nidhan
%  This code interpolates a (nr, ntheta, nx) grid to (nr/2, ntheta/2, nx/2)
%  Date - December 19, 2018

path = './'
%% Reading and postprocessing x1_grid.in

fid = fopen('x1_grid.in');
x1 = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));

count = 1;
for i = 1:2:size(x1,1)
    x1_half(count,1) = x1(i,2);
    count = count+1;
end

%% Reading and postprocessing x3_grid.in

fid = fopen('x3_grid.in');
x3 = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));

count = 1;
for i = 1:2:size(x3,1)
    x3_half(count,1) = x3(i,2);
    count = count+1;
end

%% Writing x1_grid.in

[fid]=fopen([path 'x1_grid_half.in'],'w');
    
fprintf(fid,'%.0f\n',size(x1_half,1));
    
    
for ii=1:size(x1_half,1)
    fprintf(fid, '%.0f %.6f \n',[ii x1_half(ii)]);    
end
    
fclose(fid);

%% Writing x3_grid.in

[fid]=fopen([path 'x3_grid_half.in'],'w');
    
fprintf(fid,'%.0f\n',size(x3_half,1));
    
    
for ii=1:size(x3_half,1)
    fprintf(fid, '%.0f %.6f \n',[ii x3_half(ii)]);    
end
    
fclose(fid);
