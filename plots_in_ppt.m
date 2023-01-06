%% Create ppt for plots and add plots pictures

indir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_DATA_NEW';

import mlreportgen.ppt.*;
ppt = Presentation('\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_DATA_NEW\DVL_030_T10\DVL_030_T10_plots.pptx');
open(ppt);

% Add a slide to the presentation
slide1 = add(ppt,'Title and Content');
slide2 = add(ppt,'Title and Content');
slide3 = add(ppt,'Title and Content');
slide4 = add(ppt,'Title and Content');

% Create a Picture object using an airplane image available in MATLAB.
 % Specify the size of the picture.
pic1 = Picture(which('DVL_030_T10_DEV1_allSTD.jpg'));
pic2 = Picture(which('DVL_030_T10_DEV1.jpg'));
pic3 = Picture(which('DVL_030_T10_DEV2_allSTD.jpg'));
pic4 = Picture(which('DVL_030_T10_DEV2.jpg'));
%    plane.Width = "5in";
%    plane.Height = "2in";

% Add the plane picture to the slide
replace(slide1,'Content',pic1);
replace(slide2,'Content',pic2);
replace(slide3,'Content',pic3);
replace(slide4,'Content',pic4);

% Close and view the presentation
close(ppt);
%    rptview(ppt);