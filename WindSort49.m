function [MPTF,lfor] = WindSort49(wsfinal,mpwfinal)
% JAH 11-2020
bfl = [0.3 1.6 3.4 5.5 8.0 10.8 13.9 17.2 20.8 24.5 28.5 32.7];
bfh = [1.6 3.4 5.5 8.0 10.8 13.9 17.2 20.8 24.5 28.5 32.7 50];
%sort into speed bins
force = cell(1,12);
MPTF = force;
for i = 1:12
     force{i} = find(wsfinal > bfl(i)  & wsfinal <= bfh(i) );
      MPTF{i} = mpwfinal(:,force{i});
      indx = length(force{i});
     if (indx > 0)
         lfor = i;
     end
end
   
  