function lesson=ess_live_setup
lesson.home=fileparts(mfilename('fullpath'));lesson.cases=fullfile(fileparts(lesson.home),'cases');lesson.results=fullfile(lesson.home,'results');if ~isfolder(lesson.results),mkdir(lesson.results);end;addpath(lesson.cases);
end
