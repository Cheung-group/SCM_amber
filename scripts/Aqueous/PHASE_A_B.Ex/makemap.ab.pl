#!/usr/bin/perl -w

# sidechain-sidechain interactions only. for sidechain-backbone model
# exclude backbone-backbone, sidechain backbone contact
# makemap.ab.pl a (a is the file from resP,makemap.csh) 
open(MYFILE,$ARGV[0]) || die ("open a failed"); 
unless (open (OUTFILE, ">./MapSC")){die("cant open output file\n");}

%map=();

 $p=0;
OUTER: while(<MYFILE>){

	chomp;

		if(/Residues/){
			$firstaatype = substr($_,25,3); 
			$firstaano = substr($_,29,3); 
		#print $firstaatype, $firstaano;
		#exit;
		}

		if(/Table IV/){ $p=1;}

		if($p == 1)
		{
			for($i=0;$i<8;$i++)
			{	
			$_=<MYFILE>;
		#	print $_;
			}
		LABEL: if(/\w/)
			{
			$firstatomtype = substr($_,2,2); 
			$secondatomtype = substr($_,33,2); 

			if($firstatomtype =~ /(N |CA|C |O )/ or $secondatomtype =~ /(N |CA|C |O )/ ) {
					$_=<MYFILE>;
		  	print "$firstaatype $firstaano not qualified ", $firstatomtype, $secondatomtype,"\n";
			
					goto LABEL;
					}
				else{
					$secondaatype = substr($_,19,3); 
					$secondaano = substr($_,24,3); 
					$sep=$secondaano-$firstaano;
		  	print "$firstaatype $firstaano $secondaatype $secondaano qualified ", $firstatomtype, $secondatomtype,"\n";
					if($sep>3)
					{
					$map{$firstaano}{$secondaano}=$firstaatype." ".$secondaatype;	
					#print OUTFILE "$map{$firstaano}{$secondaano} $firstaano $secondaano \n";
					}
					$_=<MYFILE>;
					goto LABEL;
					}	
			}
			else{$p=0;  goto OUTER;}	
		}
}
	
       foreach $first (sort keys %map) {
                foreach $second ( sort keys %{$map{$first}}) {
	print OUTFILE "$map{$first}{$second} $first $second \n";
		}
        }
	

