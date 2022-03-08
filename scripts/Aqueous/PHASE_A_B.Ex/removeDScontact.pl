#!/usr/bin/perl -w
################################################################################
#     A script to remove the SC contacts that have Disulfid bonds              #
#		#######   Dirar Homouz 7/19/2012   #######                     #
################################################################################

if(@ARGV < 2) {
        print "Usage: removeDScontac.pl dsfile SCmap\n";
        exit;
        }

$dsfile = shift @ARGV;
$scm = shift @ARGV;
%ds = ();
open(DSF,"$dsfile") or die "failed to open $dsfile, ";
open(SCM,"$scm") or die "failed to open $scm, ";

my $nds = <DSF>;
if($nds > 0){
  while(<DSF>){
	@dsbond = split (" ",$_);
	$aa = $dsbond[0];	
	$bb = $dsbond[1];	
	$ds{$aa}{$bb} = 1.;
	#print $aa,"\t",$bb,"\n";
  }
  $count = 1;
  while(<SCM>){
	chomp;
	@l= split (" ",$_);
	#print $l[0],"\t",$l[1],"\t",$l[2],"\t",$l[3],"\n";
	#print $ds{$l[2]}{$l[3]},"\n";
	if(not defined $ds{$l[2]}{$l[3]}){
		print $count," ",$l[1]," ",$l[2]," ",$l[3],"\n";
		#print $l[0]," ",$l[1]," ",$l[2]," ",$l[3],"\n";
		$count++;
	}
  }
}
else{
  while(<SCM>){
	printf $_;
  }
}
