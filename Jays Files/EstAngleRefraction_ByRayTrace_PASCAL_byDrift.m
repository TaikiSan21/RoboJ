% Propagation of sound from DASBR at 110m at various angles
%   for estimation of angular deviation of arriving signals
%   at an assumed depth of 1 km

% 
% Setup the sound speed profile
files=dir('j:/jay/acoustic/propagation Modeling/PASCAL PropagationModels/SoundSpeedProfileDrift*.csv')
angle(1,1)= 0;      %initialize matrix of output values (true source angle for measured declination angles 10:85)
for (iFile = 1:length(files))

     M= csvread(files(iFile).name,1,1);
     zz= M(:,1);
     cc= M(:,2);
     %angles= [80 70 60 50 40 30 20 10]
     angles= 5:85
     hpDepth= 110
     animalDepth= 1000+1;
    % % Conduct the raytrace and plot result
     [x z t d] = raytrace(0,hpDepth,angles,7,zz,cc,false);

     range= x{1};
     depth= z{1};

     for (iAngle = 1:length(angles))
        angle(iAngle,iFile)= atan((z{iAngle}(animalDepth)-hpDepth)/x{iAngle}(animalDepth))*180/pi();
     end
     angle
end
csvwrite('j:/jay/acoustic/propagation Modeling/PASCAL PropagationModels/RefractedAngleCorrections.csv',angle)
 
 
 
