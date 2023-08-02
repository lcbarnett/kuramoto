% Kuramutto connectivity spec
%
% To draw, use drawdog (hacked from wgraph2dot (see https://github.com/lcbarnett/gvmat)

% Body parts

% head

ohead = 0;
nhead = 10;
ihead = ohead+[1:nhead];

elist{ihead( 1)} = ihead([3 4]);
elist{ihead( 2)} = ihead([4 5]);
elist{ihead( 3)} = ihead([1 4 6 8]);
elist{ihead( 4)} = ihead([1 2 3 5 7]);
elist{ihead( 5)} = ihead([2 4 10]);
elist{ihead( 6)} = ihead([3 9]);
elist{ihead( 7)} = ihead([4 10]);
elist{ihead( 8)} = ihead([3 9]);
elist{ihead( 9)} = ihead([6 8 10]);
elist{ihead(10)} = ihead([5 7 9]);

% body

obodi = ohead+nhead;
nbodi = 11;
ibodi = obodi+[1:nbodi];

elist{ibodi( 1)} = ibodi([2 3 4]);
elist{ibodi( 2)} = ibodi([1 3 4]);
elist{ibodi( 3)} = ibodi([1 2 4 5 6 8]);
elist{ibodi( 4)} = ibodi([1 2 3 5 6 7]);
elist{ibodi( 5)} = ibodi([3 4 6 7 8]);
elist{ibodi( 6)} = ibodi([3 4 5 7 8]);
elist{ibodi( 7)} = ibodi([4 5 6 8 9 10]);
elist{ibodi( 8)} = ibodi([5 6 7 9 10]);
elist{ibodi( 9)} = ibodi([7 8 10 11]);
elist{ibodi(10)} = ibodi([7 8 9 11]);
elist{ibodi(11)} = ibodi([9 10]);

% front leg

ofleg = obodi+nbodi;
nfleg = 2;
ifleg = ofleg+[1:nfleg];

elist{ifleg( 1)} = ifleg([2]);
elist{ifleg( 2)} = ifleg([1]);

% hind leg

ohleg = ofleg+nfleg;
nhleg = 2;
ihleg = ohleg+[1:nhleg];

elist{ihleg( 1)} = ihleg([2]);
elist{ihleg( 2)} = ihleg([1]);

% tail

otail = ohleg+nhleg;
ntail = 3;
itail = otail+[1:ntail];

elist{itail( 1)} = itail([2]);
elist{itail( 2)} = itail([1 3]);
elist{itail( 3)} = itail([2]);

% Joints

% neck

elist{ibodi( 1)} = [elist{ibodi( 1)} ihead([8 9])];
elist{ibodi( 2)} = [elist{ibodi( 2)} ihead([10])];
elist{ihead( 8)} = [elist{ihead( 8)} ibodi([1])];
elist{ihead( 9)} = [elist{ihead( 9)} ibodi([1])];
elist{ihead(10)} = [elist{ihead(10)} ibodi([2])];

% front hip

elist{ibodi( 4)} = [elist{ibodi( 4)} ifleg([1])];
elist{ifleg( 1)} = [elist{ifleg( 1)} ibodi([4])];

% hind hip

elist{ibodi( 8)} = [elist{ibodi( 8)} ihleg([1])];
elist{ihleg( 1)} = [elist{ihleg( 1)} ibodi([8])];

% coccyx

elist{ibodi(11)} = [elist{ibodi(11)} itail([1])];
elist{itail( 1)} = [elist{itail( 1)} ibodi([11])];

[C,nosc,ncon] = elist2conmat(elist);
