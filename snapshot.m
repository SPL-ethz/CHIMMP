function snapshot(source,callbackdata)

[file path] = uiputfile(...
                {'*.pdf';'*.eps';'*.tif  [compressed]';'*.tif  [not compressed]';'*.png'},...
                 'Save as');
[pathstr,name,ext] = fileparts(file);

%save as a vector graphic file
if strcmp(ext,'.pdf')
    orient landscape 
    print('-painters','-noui',[path file],'-dpdf'); 
elseif strcmp(ext,'.eps')
    orient landscape 
    print('-painters','-noui',[path file],'-depsc2');
%save as a raster graphic file   
elseif strcmp(ext,'.png')
    orient landscape 
    print('-noui',[path file],'-dpng', '-r300'); 
elseif strcmp(ext, '.tif  [compressed]')
    file = [name,'.tif'];
    orient landscape
    print('-noui',[path file],'-dtiff','-r300');
elseif strcmp(ext,'*.tif  [not compressed]')
    file = [name,'.tif'];
    orient landscape
    print('-noui',[path file],'-dtiffn','-r300');
end
end


% [~,printers] = findprinters