clear all
load fcc100a108.txt
firstcol = fcc100a108(:, 1);
secondcol = fcc100a108(:, 2);
thirdcol = fcc100a108(:, 3);
nat = numel(firstcol);
x = zeros(nat,1);
y = zeros(nat,1);
z = zeros(nat,1);
for i=1:nat
    x(i)=firstcol(i);
    y(i)=secondcol(i);
    z(i)=thirdcol(i);
end
