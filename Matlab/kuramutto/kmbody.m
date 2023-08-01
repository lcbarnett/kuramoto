% Kuramutto connectivity spec
%
% To draw, use drawdog (hacked from wgraph2dot (see https://github.com/lcbarnett/gvmat)


% Body parts

% head

ohead = 0;
nhead = 10;
elist{ohead +  1} = ohead + [3 4];
elist{ohead +  2} = ohead + [4 5];
elist{ohead +  3} = ohead + [1 4 6 8];
elist{ohead +  4} = ohead + [1 2 3 5 7];
elist{ohead +  5} = ohead + [2 4 10];
elist{ohead +  6} = ohead + [3 9];
elist{ohead +  7} = ohead + [4 10];
elist{ohead +  8} = ohead + [3 9];
elist{ohead +  9} = ohead + [6 8 10];
elist{ohead + 10} = ohead + [5 7 9];

% body

obodi = ohead+nhead;
nbodi = 11;
elist{obodi +  1} = obodi + [2 3 4];
elist{obodi +  2} = obodi + [1 3 4];
elist{obodi +  3} = obodi + [1 2 4 5 6 8];
elist{obodi +  4} = obodi + [1 2 3 5 6 7];
elist{obodi +  5} = obodi + [3 4 6 7 8];
elist{obodi +  6} = obodi + [3 4 5 7 8];
elist{obodi +  7} = obodi + [4 5 6 8 9 10];
elist{obodi +  8} = obodi + [5 6 7 9 10];
elist{obodi +  9} = obodi + [7 8 10 11];
elist{obodi + 10} = obodi + [7 8 9 11];
elist{obodi + 11} = obodi + [9 10];

% front leg

ofleg = obodi+nbodi;
nfleg = 2;
elist{ofleg +  1} = ofleg + [2];
elist{ofleg +  2} = ofleg + [1];

% hind leg

ohleg = ofleg+nfleg;
nhleg = 2;
elist{ohleg +  1} = ohleg + [2];
elist{ohleg +  2} = ohleg + [1];

% tail

otail = ohleg+nhleg;
ntail = 3;
elist{otail +  1} = otail + [2];
elist{otail +  2} = otail + [1 3];
elist{otail +  3} = otail + [2];

% Joints

% neck

elist{obodi +  1} = [elist{obodi +  1} ohead + [8 9]];
elist{obodi +  2} = [elist{obodi +  2} ohead + [10]];
elist{ohead +  8} = [elist{ohead +  8} obodi + [1]];
elist{ohead +  9} = [elist{ohead +  9} obodi + [1]];
elist{ohead + 10} = [elist{ohead + 10} obodi + [2]];

% front hip

elist{obodi +  4} = [elist{obodi +  4} ofleg + [1]];
elist{ofleg +  1} = [elist{ofleg +  1} obodi + [4]];

% hind hip

elist{obodi +  8} = [elist{obodi +  8} ohleg + [1]];
elist{ohleg +  1} = [elist{ohleg +  1} obodi + [8]];

% coccyx

elist{obodi + 11} = [elist{obodi + 11} otail + [1]];
elist{otail +  1} = [elist{otail +  1} obodi + [11]];

C = elist2conmat(elist);

function C = elist2conmat(elist)

	n = length(elist);
	C = zeros(n);
	for e = 1:n
		C(elist{e},e) = 1;
	end

end
