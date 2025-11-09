function [pin] = in_wake(xy,n,k,D,pin)
% Which WTs are within wake for different angles?
% Calculated for one selected WT

% Other WTs in wake of selected WT
for d = 1:360
   G = [cosd(d) sind(d);-sind(d) cosd(d)]; 
   xy2 = xy*G; %rotate coordinate system
   dx = D/2+k*abs(xy2(n,2)-xy2(:,2));
   pin1(:,d) = abs(xy2(n,1)-xy2(:,1))<dx & xy2(n,2)>xy2(:,2);
end
pin1(n,:) = 0;

% Selected WT in wake of other WTs
for d = 1:360
   G = [cosd(d) sind(d);-sind(d) cosd(d)]; 
   xy2 = xy*G; %rotate coordinate system
   dx = D/2+k*abs(xy2(n,2)-xy2(:,2));
   pin2(d,:) = abs(xy2(n,1)-xy2(:,1))<dx & xy2(n,2)<xy2(:,2);
end
pin2(:,n) = 0;

pin(:,:,n) = pin1;
pin(n,:,:) = pin2;

end

