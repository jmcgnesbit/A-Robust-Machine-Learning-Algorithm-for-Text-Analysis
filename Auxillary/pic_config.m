classdef pic_config
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pic_config Class for pictures
    %
    % Inputs (for constructor)
    % fontsize
    % height
    % width
    % directory 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   properties 
    fontsize =12
    height = 4
    width = 10
    directory
   end
   % Constructor
   methods
       function obj = pic_config(fontsize, height, width, directory)
           obj.fontsize = fontsize;
           obj.height = height;
           obj.width = width;
           obj.directory = directory;
       end
        
      end
%    events (Attributes) 
%       EventName
%    end
%    enumeration
%       EnumName
%    end
end