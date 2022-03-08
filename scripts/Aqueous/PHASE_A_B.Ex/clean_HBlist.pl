#!/usr/bin/perl
#usage:  ./clean_HBlist.pl HBondlist 147
#Script to clean the HBlist

open DF, $ARGV[0];

while (<DF>)
{
	chomp;	
	$line = $_;
	@flds = split" ", $line;
	$a = $flds[0];
	$b = $flds[2];
	$en = $flds[3];
	if(defined $pair{$a}{$b}){
		if($pair{$a}{$b} > $en)  {$pair{$a}{$b}= $en; $pair{$b}{$a}= $en;}
	}
	else {$pair{$a}{$b}= $en; $pair{$b}{$a}= $en; }
};
close DF;
for $x (sort {$a<=>$b} keys %pair ) {
    %tmp = ();
    @yy = ();
    #print "<<<<<<<<< ",$x," >>>>>>>>>>>\n";
    for $y (sort {$a<=>$b} keys %{$pair{$x}}) {
#	print $x,"\t",$y,"\t",$pair{$x}{$y},"\n";
	$tmp{$y} = $pair{$x}{$y};	
    }
    @yy= reverse sort {$tmp{$b} <=> $tmp{$a}} keys %tmp;
    $nk = 1;
    for $yyy (@yy){
	if($nk > 2) {
		if(defined $pair{$x}{$yyy}){undef $pair{$x}{$yyy}; undef $pair{$yyy}{$x};}
	}
#	else {print $yyy,"\t",$tmp{$yyy},"\n";}
	$nk++;
    }
}
 #   print "<<<<<<<<< end  >>>>>>>>>>>\n";
for $x (sort {$a<=>$b} keys %pair) {
    for $y (sort {$a<=>$b} keys %{$pair{$x}}) {
	if(defined $pair{$x}{$y} && $y > $x) {print $x,"\t",$y,"\t",$pair{$x}{$y},"\n";}
    }
}
