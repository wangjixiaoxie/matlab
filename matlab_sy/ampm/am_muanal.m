%!PS-Adobe-2.0 EPSF-1.2
%%Creator: MATLAB, The Mathworks, Inc. Version 7.5.0.338 (R2007b). Operating System: Linux 2.6.28.4-netboot-a #1 SMP Wed Feb 11 14:06:44 PST 2009 x86_64.
%%Title: ./am_muanal.m
%%CreationDate: 01/27/2010  18:20:43
%%DocumentNeededFonts: Helvetica
%%DocumentProcessColors: Cyan Magenta Yellow Black
%%Extensions: CMYK
%%Pages: 1
%%BoundingBox:    76   215   572   589
%%EndComments

%%BeginProlog
% MathWorks dictionary
/MathWorks 160 dict begin
% definition operators
/bdef {bind def} bind def
/ldef {load def} bind def
/xdef {exch def} bdef
/xstore {exch store} bdef
% operator abbreviations
/c  /clip ldef
/cc /concat ldef
/cp /closepath ldef
/gr /grestore ldef
/gs /gsave ldef
/mt /moveto ldef
/np /newpath ldef
/cm /currentmatrix ldef
/sm /setmatrix ldef
/rm /rmoveto ldef
/rl /rlineto ldef
/s {show newpath} bdef
/sc {setcmykcolor} bdef
/sr /setrgbcolor ldef
/sg /setgray ldef
/w /setlinewidth ldef
/j /setlinejoin ldef
/cap /setlinecap ldef
/rc {rectclip} bdef
/rf {rectfill} bdef
% page state control
/pgsv () def
/bpage {/pgsv save def} bdef
/epage {pgsv restore} bdef
/bplot /gsave ldef
/eplot {stroke grestore} bdef
% orientation switch
/portraitMode 0 def /landscapeMode 1 def /rotateMode 2 def
% coordinate system mappings
/dpi2point 0 def
% font control
/FontSize 0 def
/FMS {/FontSize xstore findfont [FontSize 0 0 FontSize neg 0 0]
  makefont setfont} bdef
/reencode {exch dup where {pop load} {pop StandardEncoding} ifelse
  exch dup 3 1 roll findfont dup length dict begin
  { 1 index /FID ne {def}{pop pop} ifelse } forall
  /Encoding exch def currentdict end definefont pop} bdef
/isroman {findfont /CharStrings get /Agrave known} bdef
/FMSR {3 1 roll 1 index dup isroman {reencode} {pop pop} ifelse
  exch FMS} bdef
/csm {1 dpi2point div -1 dpi2point div scale neg translate
 dup landscapeMode eq {pop -90 rotate}
  {rotateMode eq {90 rotate} if} ifelse} bdef
% line types: solid, dotted, dashed, dotdash
/SO { [] 0 setdash } bdef
/DO { [.5 dpi2point mul 4 dpi2point mul] 0 setdash } bdef
/DA { [6 dpi2point mul] 0 setdash } bdef
/DD { [.5 dpi2point mul 4 dpi2point mul 6 dpi2point mul 4
  dpi2point mul] 0 setdash } bdef
% macros for lines and objects
/L {lineto stroke} bdef
/MP {3 1 roll moveto 1 sub {rlineto} repeat} bdef
/AP {{rlineto} repeat} bdef
/PDlw -1 def
/W {/PDlw currentlinewidth def setlinewidth} def
/PP {closepath eofill} bdef
/DP {closepath stroke} bdef
/MR {4 -2 roll moveto dup  0 exch rlineto exch 0 rlineto
  neg 0 exch rlineto closepath} bdef
/FR {MR stroke} bdef
/PR {MR fill} bdef
/L1i {{currentfile picstr readhexstring pop} image} bdef
/tMatrix matrix def
/MakeOval {newpath tMatrix currentmatrix pop translate scale
0 0 1 0 360 arc tMatrix setmatrix} bdef
/FO {MakeOval stroke} bdef
/PO {MakeOval fill} bdef
/PD {currentlinewidth 2 div 0 360 arc fill
   PDlw -1 eq not {PDlw w /PDlw -1 def} if} def
/FA {newpath tMatrix currentmatrix pop translate scale
  0 0 1 5 -2 roll arc tMatrix setmatrix stroke} bdef
/PA {newpath tMatrix currentmatrix pop	translate 0 0 moveto scale
  0 0 1 5 -2 roll arc closepath tMatrix setmatrix fill} bdef
/FAn {newpath tMatrix currentmatrix pop translate scale
  0 0 1 5 -2 roll arcn tMatrix setmatrix stroke} bdef
/PAn {newpath tMatrix currentmatrix pop translate 0 0 moveto scale
  0 0 1 5 -2 roll arcn closepath tMatrix setmatrix fill} bdef
/vradius 0 def /hradius 0 def /lry 0 def
/lrx 0 def /uly 0 def /ulx 0 def /rad 0 def
/MRR {/vradius xdef /hradius xdef /lry xdef /lrx xdef /uly xdef
  /ulx xdef newpath tMatrix currentmatrix pop ulx hradius add uly
  vradius add translate hradius vradius scale 0 0 1 180 270 arc 
  tMatrix setmatrix lrx hradius sub uly vradius add translate
  hradius vradius scale 0 0 1 270 360 arc tMatrix setmatrix
  lrx hradius sub lry vradius sub translate hradius vradius scale
  0 0 1 0 90 arc tMatrix setmatrix ulx hradius add lry vradius sub
  translate hradius vradius scale 0 0 1 90 180 arc tMatrix setmatrix
  closepath} bdef
/FRR {MRR stroke } bdef
/PRR {MRR fill } bdef
/MlrRR {/lry xdef /lrx xdef /uly xdef /ulx xdef /rad lry uly sub 2 div def
  newpath tMatrix currentmatrix pop ulx rad add uly rad add translate
  rad rad scale 0 0 1 90 270 arc tMatrix setmatrix lrx rad sub lry rad
  sub translate rad rad scale 0 0 1 270 90 arc tMatrix setmatrix
  closepath} bdef
/FlrRR {MlrRR stroke } bdef
/PlrRR {MlrRR fill } bdef
/MtbRR {/lry xdef /lrx xdef /uly xdef /ulx xdef /rad lrx ulx sub 2 div def
  newpath tMatrix currentmatrix pop ulx rad add uly rad add translate
  rad rad scale 0 0 1 180 360 arc tMatrix setmatrix lrx rad sub lry rad
  sub translate rad rad scale 0 0 1 0 180 arc tMatrix setmatrix
  closepath} bdef
/FtbRR {MtbRR stroke } bdef
/PtbRR {MtbRR fill } bdef
/stri 6 array def /dtri 6 array def
/smat 6 array def /dmat 6 array def
/tmat1 6 array def /tmat2 6 array def /dif 3 array def
/asub {/ind2 exch def /ind1 exch def dup dup
  ind1 get exch ind2 get sub exch } bdef
/tri_to_matrix {
  2 0 asub 3 1 asub 4 0 asub 5 1 asub
  dup 0 get exch 1 get 7 -1 roll astore } bdef
/compute_transform {
  dmat dtri tri_to_matrix tmat1 invertmatrix 
  smat stri tri_to_matrix tmat2 concatmatrix } bdef
/ds {stri astore pop} bdef
/dt {dtri astore pop} bdef
/db {2 copy /cols xdef /rows xdef mul dup 3 mul string
  currentfile exch readhexstring pop
  dup 0 3 index getinterval /rbmap xdef
  dup 2 index dup getinterval /gbmap xdef
  1 index dup 2 mul exch getinterval /bbmap xdef pop pop}bdef
/it {gs np dtri aload pop moveto lineto lineto cp c
  cols rows 8 compute_transform 
  rbmap gbmap bbmap true 3 colorimage gr}bdef
/il {newpath moveto lineto stroke}bdef
currentdict end def
%%EndProlog

%%BeginSetup
MathWorks begin

0 cap

end
%%EndSetup

%%Page: 1 1
%%BeginPageSetup
%%PageBoundingBox:    76   215   572   589
MathWorks begin
bpage
%%EndPageSetup

%%BeginObject: obj1
bplot

/dpi2point 12 def
portraitMode 0216 7344 csm

  698   274  5960  4487 MR c np
85 dict begin %Colortable dictionary
/c0 { 0.000000 0.000000 0.000000 sr} bdef
/c1 { 1.000000 1.000000 1.000000 sr} bdef
/c2 { 0.900000 0.000000 0.000000 sr} bdef
/c3 { 0.000000 0.820000 0.000000 sr} bdef
/c4 { 0.000000 0.000000 0.800000 sr} bdef
/c5 { 0.910000 0.820000 0.320000 sr} bdef
/c6 { 1.000000 0.260000 0.820000 sr} bdef
/c7 { 0.000000 0.820000 0.820000 sr} bdef
c0
1 j
1 sg
   0    0 6917 5186 PR
6 w
0 1118 5360 0 0 -1118 899 1507 4 MP
PP
-5360 0 0 1118 5360 0 0 -1118 899 1507 5 MP stroke
4 w
DO
SO
6 w
0 sg
 899 1507 mt 6259 1507 L
 899  389 mt 6259  389 L
 899 1507 mt  899  389 L
6259 1507 mt 6259  389 L
 899 1507 mt 6259 1507 L
 899 1507 mt  899  389 L
 899 1507 mt  899 1453 L
 899  389 mt  899  442 L
%%IncludeResource: font Helvetica
/Helvetica /ISOLatin1Encoding 120 FMSR

 796 1652 mt 
(-5) s
1435 1507 mt 1435 1453 L
1435  389 mt 1435  442 L
1332 1652 mt 
(-4) s
1971 1507 mt 1971 1453 L
1971  389 mt 1971  442 L
1868 1652 mt 
(-3) s
2507 1507 mt 2507 1453 L
2507  389 mt 2507  442 L
2404 1652 mt 
(-2) s
3043 1507 mt 3043 1453 L
3043  389 mt 3043  442 L
2940 1652 mt 
(-1) s
3579 1507 mt 3579 1453 L
3579  389 mt 3579  442 L
3546 1652 mt 
(0) s
4115 1507 mt 4115 1453 L
4115  389 mt 4115  442 L
4082 1652 mt 
(1) s
4651 1507 mt 4651 1453 L
4651  389 mt 4651  442 L
4618 1652 mt 
(2) s
5187 1507 mt 5187 1453 L
5187  389 mt 5187  442 L
5154 1652 mt 
(3) s
5723 1507 mt 5723 1453 L
5723  389 mt 5723  442 L
5690 1652 mt 
(4) s
6259 1507 mt 6259 1453 L
6259  389 mt 6259  442 L
6226 1652 mt 
(5) s
 899 1507 mt  952 1507 L
6259 1507 mt 6205 1507 L
 798 1551 mt 
(0) s
 899 1227 mt  952 1227 L
6259 1227 mt 6205 1227 L
 698 1271 mt 
(0.1) s
 899  948 mt  952  948 L
6259  948 mt 6205  948 L
 698  992 mt 
(0.2) s
 899  668 mt  952  668 L
6259  668 mt 6205  668 L
 698  712 mt 
(0.3) s
 899  389 mt  952  389 L
6259  389 mt 6205  389 L
 698  433 mt 
(0.4) s
 899 1507 mt 6259 1507 L
 899  389 mt 6259  389 L
 899 1507 mt  899  389 L
6259 1507 mt 6259  389 L
%%IncludeResource: font Helvetica
/Helvetica /ISOLatin1Encoding 192 FMSR

5187 1155 mt 
(mnout= NaN) s
5187 1211 mt 
(ste= NaN) s
5187 1295 mt 
(nvls= 0) s
gs 899 389 5361 1119 MR c np
36 w
/c8 { 1.000000 0.000000 0.000000 sr} bdef
c8
0 65 187 0 0 -65 188 0 187 0 188 0 188 0 187 0 
0 130 188 0 0 -130 187 0 0 260 188 0 0 130 188 0 
0 -260 187 0 0 260 188 0 0 -130 188 0 0 195 187 0 
0 -65 188 0 0 -195 187 0 0 -65 188 0 0 -130 187 0 
188 0 188 0 187 0 188 0 187 0 188 0 188 0 187 0 
188 0 187 0 899 1507 43 MP stroke
gr

36 w
c8
5187 1155 mt 
(mnout= 0.68558) s
5187 1211 mt 
(ste= 0.16358) s
5187 1295 mt 
(nvls= 43) s
0 sg
5187 1155 mt 
(mnout= NaN) s
5187 1211 mt 
(ste= NaN) s
5187 1295 mt 
(nvls= 0) s
gs 899 389 5361 1119 MR c np
c8
0 65 187 0 0 -65 188 0 187 0 188 0 188 0 187 0 
0 130 188 0 0 -130 187 0 0 260 188 0 0 130 188 0 
0 -260 187 0 0 260 188 0 0 -130 188 0 0 195 187 0 
0 -65 188 0 0 -195 187 0 0 -65 188 0 0 -130 187 0 
188 0 188 0 187 0 188 0 187 0 188 0 188 0 187 0 
188 0 187 0 899 1507 43 MP stroke
gr

c8
5187 1155 mt 
(mnout= 0.68558) s
5187 1211 mt 
(ste= 0.16358) s
5187 1295 mt 
(nvls= 43) s
0 sg
5187 1155 mt 
(mnout= NaN) s
5187 1211 mt 
(ste= NaN) s
5187 1295 mt 
(nvls= 0) s
5187 1155 mt 
(mnout= NaN) s
5187 1211 mt 
(ste= NaN) s
5187 1295 mt 
(nvls= 0) s
gs 899 389 5361 1119 MR c np
c8
0 65 187 0 0 -65 188 0 187 0 188 0 188 0 187 0 
0 130 188 0 0 -130 187 0 0 260 188 0 0 130 188 0 
0 -260 187 0 0 260 188 0 0 -130 188 0 0 195 187 0 
0 -65 188 0 0 -195 187 0 0 -65 188 0 0 -130 187 0 
188 0 188 0 187 0 188 0 187 0 188 0 188 0 187 0 
188 0 187 0 899 1507 43 MP stroke
gr

c8
5187 1155 mt 
(mnout= 0.68558) s
5187 1211 mt 
(ste= 0.16358) s
5187 1295 mt 
(nvls= 43) s
0 sg
5187 1155 mt 
(mnout= NaN) s
5187 1211 mt 
(ste= NaN) s
5187 1295 mt 
(nvls= 0) s
gs 899 389 5361 1119 MR c np
187 0 188 0 187 0 188 0 0 100 188 0 0 -100 187 0 
0 300 188 0 0 100 187 0 0 -200 188 0 0 100 188 0 
187 0 0 399 188 0 0 -399 188 0 0 -100 187 0 0 -200 
188 0 187 0 188 0 187 0 188 0 188 0 187 0 188 0 
187 0 188 0 188 0 187 0 188 0 187 0 899 1507 39 MP stroke
gr

5187 1155 mt 
(mnout= 1.3598) s
5187 1211 mt 
(ste= 0.1584) s
5187 1295 mt 
(nvls= 28) s
gs 899 389 5361 1119 MR c np
c8
0 65 187 0 0 -65 188 0 187 0 188 0 188 0 187 0 
0 130 188 0 0 -130 187 0 0 260 188 0 0 130 188 0 
0 -260 187 0 0 260 188 0 0 -130 188 0 0 195 187 0 
0 -65 188 0 0 -195 187 0 0 -65 188 0 0 -130 187 0 
188 0 188 0 187 0 188 0 187 0 188 0 188 0 187 0 
188 0 187 0 899 1507 43 MP stroke
gr

c8
5187 1155 mt 
(mnout= 0.68558) s
5187 1211 mt 
(ste= 0.16358) s
5187 1295 mt 
(nvls= 43) s
6 w
1 sg
0 1119 5360 0 0 -1119 899 3061 4 MP
PP
-5360 0 0 1119 5360 0 0 -1119 899 3061 5 MP stroke
4 w
DO
SO
6 w
0 sg
 899 3061 mt 6259 3061 L
 899 1942 mt 6259 1942 L
 899 3061 mt  899 1942 L
6259 3061 mt 6259 1942 L
 899 3061 mt 6259 3061 L
 899 3061 mt  899 1942 L
 899 3061 mt  899 3007 L
 899 1942 mt  899 1995 L
%%IncludeResource: font Helvetica
/Helvetica /ISOLatin1Encoding 120 FMSR

 796 3206 mt 
(-5) s
1435 3061 mt 1435 3007 L
1435 1942 mt 1435 1995 L
1332 3206 mt 
(-4) s
1971 3061 mt 1971 3007 L
1971 1942 mt 1971 1995 L
1868 3206 mt 
(-3) s
2507 3061 mt 2507 3007 L
2507 1942 mt 2507 1995 L
2404 3206 mt 
(-2) s
3043 3061 mt 3043 3007 L
3043 1942 mt 3043 1995 L
2940 3206 mt 
(-1) s
3579 3061 mt 3579 3007 L
3579 1942 mt 3579 1995 L
3546 3206 mt 
(0) s
4115 3061 mt 4115 3007 L
4115 1942 mt 4115 1995 L
4082 3206 mt 
(1) s
4651 3061 mt 4651 3007 L
4651 1942 mt 4651 1995 L
4618 3206 mt 
(2) s
5187 3061 mt 5187 3007 L
5187 1942 mt 5187 1995 L
5154 3206 mt 
(3) s
5723 3061 mt 5723 3007 L
5723 1942 mt 5723 1995 L
5690 3206 mt 
(4) s
6259 3061 mt 6259 3007 L
6259 1942 mt 6259 1995 L
6226 3206 mt 
(5) s
 899 3061 mt  952 3061 L
6259 3061 mt 6205 3061 L
 798 3105 mt 
(0) s
 899 2781 mt  952 2781 L
6259 2781 mt 6205 2781 L
 698 2825 mt 
(0.1) s
 899 2501 mt  952 2501 L
6259 2501 mt 6205 2501 L
 698 2545 mt 
(0.2) s
 899 2221 mt  952 2221 L
6259 2221 mt 6205 2221 L
 698 2265 mt 
(0.3) s
 899 1942 mt  952 1942 L
6259 1942 mt 6205 1942 L
 698 1986 mt 
(0.4) s
 899 3061 mt 6259 3061 L
 899 1942 mt 6259 1942 L
 899 3061 mt  899 1942 L
6259 3061 mt 6259 1942 L
gs 899 1942 5361 1120 MR c np
36 w
187 0 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
188 0 188 0 187 0 188 0 0 933 188 0 0 -933 187 0 
188 0 187 0 0 933 188 0 187 0 0 -933 188 0 188 0 
187 0 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
899 3061 33 MP stroke
gr

36 w
%%IncludeResource: font Helvetica
/Helvetica /ISOLatin1Encoding 192 FMSR

5187 2709 mt 
(mnout= -0.6298) s
5187 2765 mt 
(ste= 0.55953) s
5187 2849 mt 
(nvls= 3) s
gs 899 1942 5361 1120 MR c np
c8
187 0 188 0 187 0 188 0 188 0 187 0 188 0 0 165 
187 0 0 -165 188 0 0 330 188 0 187 0 0 -165 188 0 
0 329 188 0 0 -164 187 0 0 -165 188 0 0 329 187 0 
0 -329 188 0 0 -165 187 0 188 0 188 0 0 165 187 0 
0 -165 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
899 3061 41 MP stroke
gr

c8
5187 2709 mt 
(mnout= 0.24895) s
5187 2765 mt 
(ste= 0.2596) s
5187 2849 mt 
(nvls= 17) s
gs 899 1942 5361 1120 MR c np
0 sg
187 0 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
188 0 188 0 187 0 188 0 0 933 188 0 0 -933 187 0 
188 0 187 0 0 933 188 0 187 0 0 -933 188 0 188 0 
187 0 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
899 3061 33 MP stroke
gr

0 sg
5187 2709 mt 
(mnout= -0.6298) s
5187 2765 mt 
(ste= 0.55953) s
5187 2849 mt 
(nvls= 3) s
gs 899 1942 5361 1120 MR c np
c8
187 0 188 0 187 0 188 0 188 0 187 0 188 0 0 165 
187 0 0 -165 188 0 0 330 188 0 187 0 0 -165 188 0 
0 329 188 0 0 -164 187 0 0 -165 188 0 0 329 187 0 
0 -329 188 0 0 -165 187 0 188 0 188 0 0 165 187 0 
0 -165 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
899 3061 41 MP stroke
gr

c8
5187 2709 mt 
(mnout= 0.24895) s
5187 2765 mt 
(ste= 0.2596) s
5187 2849 mt 
(nvls= 17) s
gs 899 1942 5361 1120 MR c np
0 sg
187 0 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
188 0 188 0 187 0 188 0 0 933 188 0 0 -933 187 0 
188 0 187 0 0 933 188 0 187 0 0 -933 188 0 188 0 
187 0 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
899 3061 33 MP stroke
gr

0 sg
5187 2709 mt 
(mnout= -0.6298) s
5187 2765 mt 
(ste= 0.55953) s
5187 2849 mt 
(nvls= 3) s
gs 899 1942 5361 1120 MR c np
c8
187 0 188 0 187 0 188 0 188 0 187 0 188 0 0 165 
187 0 0 -165 188 0 0 330 188 0 187 0 0 -165 188 0 
0 329 188 0 0 -164 187 0 0 -165 188 0 0 329 187 0 
0 -329 188 0 0 -165 187 0 188 0 188 0 0 165 187 0 
0 -165 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
899 3061 41 MP stroke
gr

c8
5187 2709 mt 
(mnout= 0.24895) s
5187 2765 mt 
(ste= 0.2596) s
5187 2849 mt 
(nvls= 17) s
gs 899 1942 5361 1120 MR c np
0 sg
187 0 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
188 0 188 0 187 0 188 0 0 933 188 0 0 -933 187 0 
188 0 187 0 0 933 188 0 187 0 0 -933 188 0 188 0 
187 0 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
899 3061 33 MP stroke
gr

0 sg
5187 2709 mt 
(mnout= -0.6298) s
5187 2765 mt 
(ste= 0.55953) s
5187 2849 mt 
(nvls= 3) s
gs 899 1942 5361 1120 MR c np
c8
187 0 188 0 187 0 188 0 188 0 187 0 188 0 0 165 
187 0 0 -165 188 0 0 330 188 0 187 0 0 -165 188 0 
0 329 188 0 0 -164 187 0 0 -165 188 0 0 329 187 0 
0 -329 188 0 0 -165 187 0 188 0 188 0 0 165 187 0 
0 -165 188 0 187 0 188 0 188 0 187 0 188 0 187 0 
899 3061 41 MP stroke
gr

c8
5187 2709 mt 
(mnout= 0.24895) s
5187 2765 mt 
(ste= 0.2596) s
5187 2849 mt 
(nvls= 17) s
6 w
1 sg
0 1119 5360 0 0 -1119 899 4615 4 MP
PP
-5360 0 0 1119 5360 0 0 -1119 899 4615 5 MP stroke
4 w
DO
SO
6 w
0 sg
 899 4615 mt 6259 4615 L
 899 3496 mt 6259 3496 L
 899 4615 mt  899 3496 L
6259 4615 mt 6259 3496 L
 899 4615 mt 6259 4615 L
 899 4615 mt  899 3496 L
 899 4615 mt  899 4561 L
 899 3496 mt  899 3549 L
%%IncludeResource: font Helvetica
/Helvetica /ISOLatin1Encoding 120 FMSR

 796 4760 mt 
(-5) s
1435 4615 mt 1435 4561 L
1435 3496 mt 1435 3549 L
1332 4760 mt 
(-4) s
1971 4615 mt 1971 4561 L
1971 3496 mt 1971 3549 L
1868 4760 mt 
(-3) s
2507 4615 mt 2507 4561 L
2507 3496 mt 2507 3549 L
2404 4760 mt 
(-2) s
3043 4615 mt 3043 4561 L
3043 3496 mt 3043 3549 L
2940 4760 mt 
(-1) s
3579 4615 mt 3579 4561 L
3579 3496 mt 3579 3549 L
3546 4760 mt 
(0) s
4115 4615 mt 4115 4561 L
4115 3496 mt 4115 3549 L
4082 4760 mt 
(1) s
4651 4615 mt 4651 4561 L
4651 3496 mt 4651 3549 L
4618 4760 mt 
(2) s
5187 4615 mt 5187 4561 L
5187 3496 mt 5187 3549 L
5154 4760 mt 
(3) s
5723 4615 mt 5723 4561 L
5723 3496 mt 5723 3549 L
5690 4760 mt 
(4) s
6259 4615 mt 6259 4561 L
6259 3496 mt 6259 3549 L
6226 4760 mt 
(5) s
 899 4615 mt  952 4615 L
6259 4615 mt 6205 4615 L
 798 4659 mt 
(0) s
 899 4335 mt  952 4335 L
6259 4335 mt 6205 4335 L
 698 4379 mt 
(0.1) s
 899 4055 mt  952 4055 L
6259 4055 mt 6205 4055 L
 698 4099 mt 
(0.2) s
 899 3775 mt  952 3775 L
6259 3775 mt 6205 3775 L
 698 3819 mt 
(0.3) s
 899 3496 mt  952 3496 L
6259 3496 mt 6205 3496 L
 698 3540 mt 
(0.4) s
 899 4615 mt 6259 4615 L
 899 3496 mt 6259 3496 L
 899 4615 mt  899 3496 L
6259 4615 mt 6259 3496 L
gs 899 3496 5361 1120 MR c np
36 w
187 0 188 0 187 0 188 0 0 112 188 0 0 -112 187 0 
0 224 188 0 0 112 187 0 0 -112 188 0 0 112 188 0 
187 0 0 336 188 0 0 -336 188 0 0 -112 187 0 0 -224 
188 0 187 0 188 0 187 0 188 0 188 0 187 0 188 0 
187 0 188 0 188 0 187 0 188 0 187 0 899 4615 39 MP stroke
gr

36 w
%%IncludeResource: font Helvetica
/Helvetica /ISOLatin1Encoding 192 FMSR

5187 4263 mt 
(mnout= 1.297) s
5187 4319 mt 
(ste= 0.16497) s
5187 4403 mt 
(nvls= 25) s
gs 899 3496 5361 1120 MR c np
c8
187 0 188 0 187 0 188 0 188 0 187 0 0 216 188 0 
0 -216 187 0 0 108 188 0 0 323 188 0 0 -323 187 0 
0 108 188 0 0 107 188 0 0 323 187 0 0 -323 188 0 
187 0 0 -215 188 0 0 -108 187 0 188 0 188 0 187 0 
188 0 187 0 188 0 188 0 187 0 188 0 187 0 899 4615 40 MP stroke
gr

c8
5187 4263 mt 
(mnout= 0.51872) s
5187 4319 mt 
(ste= 0.18494) s
5187 4403 mt 
(nvls= 26) s
gs 899 3496 5361 1120 MR c np
0 sg
187 0 188 0 187 0 188 0 0 112 188 0 0 -112 187 0 
0 224 188 0 0 112 187 0 0 -112 188 0 0 112 188 0 
187 0 0 336 188 0 0 -336 188 0 0 -112 187 0 0 -224 
188 0 187 0 188 0 187 0 188 0 188 0 187 0 188 0 
187 0 188 0 188 0 187 0 188 0 187 0 899 4615 39 MP stroke
gr

0 sg
5187 4263 mt 
(mnout= 1.297) s
5187 4319 mt 
(ste= 0.16497) s
5187 4403 mt 
(nvls= 25) s
gs 899 3496 5361 1120 MR c np
c8
187 0 188 0 187 0 188 0 188 0 187 0 0 216 188 0 
0 -216 187 0 0 108 188 0 0 323 188 0 0 -323 187 0 
0 108 188 0 0 107 188 0 0 323 187 0 0 -323 188 0 
187 0 0 -215 188 0 0 -108 187 0 188 0 188 0 187 0 
188 0 187 0 188 0 188 0 187 0 188 0 187 0 899 4615 40 MP stroke
gr

c8
5187 4263 mt 
(mnout= 0.51872) s
5187 4319 mt 
(ste= 0.18494) s
5187 4403 mt 
(nvls= 26) s
gs 899 3496 5361 1120 MR c np
0 sg
187 0 188 0 187 0 188 0 0 112 188 0 0 -112 187 0 
0 224 188 0 0 112 187 0 0 -112 188 0 0 112 188 0 
187 0 0 336 188 0 0 -336 188 0 0 -112 187 0 0 -224 
188 0 187 0 188 0 187 0 188 0 188 0 187 0 188 0 
187 0 188 0 188 0 187 0 188 0 187 0 899 4615 39 MP stroke
gr

0 sg
5187 4263 mt 
(mnout= 1.297) s
5187 4319 mt 
(ste= 0.16497) s
5187 4403 mt 
(nvls= 25) s
gs 899 3496 5361 1120 MR c np
c8
187 0 188 0 187 0 188 0 188 0 187 0 0 216 188 0 
0 -216 187 0 0 108 188 0 0 323 188 0 0 -323 187 0 
0 108 188 0 0 107 188 0 0 323 187 0 0 -323 188 0 
187 0 0 -215 188 0 0 -108 187 0 188 0 188 0 187 0 
188 0 187 0 188 0 188 0 187 0 188 0 187 0 899 4615 40 MP stroke
gr

c8
5187 4263 mt 
(mnout= 0.51872) s
5187 4319 mt 
(ste= 0.18494) s
5187 4403 mt 
(nvls= 26) s
gs 899 3496 5361 1120 MR c np
0 sg
187 0 188 0 187 0 188 0 0 112 188 0 0 -112 187 0 
0 224 188 0 0 112 187 0 0 -112 188 0 0 112 188 0 
187 0 0 336 188 0 0 -336 188 0 0 -112 187 0 0 -224 
188 0 187 0 188 0 187 0 188 0 188 0 187 0 188 0 
187 0 188 0 188 0 187 0 188 0 187 0 899 4615 39 MP stroke
gr

0 sg
5187 4263 mt 
(mnout= 1.297) s
5187 4319 mt 
(ste= 0.16497) s
5187 4403 mt 
(nvls= 25) s
gs 899 3496 5361 1120 MR c np
c8
187 0 188 0 187 0 188 0 188 0 187 0 0 216 188 0 
0 -216 187 0 0 108 188 0 0 323 188 0 0 -323 187 0 
0 108 188 0 0 107 188 0 0 323 187 0 0 -323 188 0 
187 0 0 -215 188 0 0 -108 187 0 188 0 188 0 187 0 
188 0 187 0 188 0 188 0 187 0 188 0 187 0 899 4615 40 MP stroke
gr

c8
5187 4263 mt 
(mnout= 0.51872) s
5187 4319 mt 
(ste= 0.18494) s
5187 4403 mt 
(nvls= 26) s
6 w

end %%Color Dict

eplot
%%EndObject

epage
end

showpage

%%Trailer
%%EOF
