#!/usr/bin/perl -w
#Replaces the coordinates for protein only
@new = ();
@shft = ();

if (@ARGV < 2){die "usage: makecrd.pl protein+crowd.crd protein.crd"; }
$pfile = shift @ARGV;
$pcfile = shift @ARGV;


open (PIN, $pfile) or die "no input protein crd";
open (PCIN, $pcfile) or die "no input protein+crowd crd";



my $title = <PIN>;
my $num = <PIN>;
my $k = 0;
while(<PIN>){
	chomp;
	@list= split(" ", $_); 
	for $l (@list){
		$new[$k] = $l;
		$k++;
	}
}


$title = <PCIN>;
$num = <PCIN>;
print $title;
print $num;
$k = 0;
chomp($ll=<PCIN>);	#reading first line of coords
@list= split(" ", $ll);
$shft[0]= $list[0]-$new[0];	#shifts to the new location
$shft[1]= $list[1]-$new[1];
$shft[2]= $list[2]-$new[2];
for $l (@list){
	printf("%12.7f",$new[$k]+$shft[$k%3]);
 	$k++;
} 
print "\n";
while(<PCIN>){
	if($k<$#new) {
		chomp;
		@list= split(" ", $_);
		for $l (@list){
			if($k<=$#new) {printf("%12.7f",$new[$k]+$shft[$k%3]);}
			else {printf("%12.7f",$l);}
                	$k++;
        	}
		print "\n";
	}
	else {print $_;};
}
