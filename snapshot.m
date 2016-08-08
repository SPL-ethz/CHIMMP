function snapshot(source,callbackdata)

% This file is part of CHIMMP.
% 
% CHIMMP is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
% 
% CHIMMP is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

[ file, path ] = uiputfile(...
                {'*.pdf';'*.eps';'*.tif  [compressed]';'*.tif  [not compressed]';'*.png'},...
                 'Save as');
[ ~, name, ext ] = fileparts(file);
 
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