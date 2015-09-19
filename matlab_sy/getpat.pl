#!/usr/bin/perl
# first argument is the pattern being searched for
# second argument is the file with each line a seperate set of
#   labels
#
$pat = shift(@ARGV);
#print "pattern=",$pat,"\n";

$infile = shift(@ARGV);
if (-T $infile){
    #print 'input filename = ',$infile,"\n";
} else {
    print 'Could not open inputfile =',$infile,"\n";
}

#open file with input labels
open(FINP,$infile) || die "Can't open input file :",$infile;

#open file for output
open(FOUT,'> SEARCH_RESULTS.TEMP');

$filenum=0;
$wrtonce=0;
while ($label = <FINP>){
    $filenum++;$wrtonce=1;
    while ($label=~/$pat/g) {
          if ($wrtonce) {
              print FOUT "file number = ",$filenum,"\n";
              $wrtonce=0;
          }
          print FOUT "pos   = ",pos($label),"\n";
          print FOUT "match = ",$&,"\n";
    }
}
close(FINP);close(FOUT);

###old notes
#$label = shift(@ARGV);
#print $label,"\n";
#$label =~ m/$pat/g;
#print '$& =',$&,"\n";
#print 'found at position = ',pos($label),' on line = ',$filenum," \n";
