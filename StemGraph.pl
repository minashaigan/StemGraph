	use warnings;
	use strict;


	open (my $in,"RnaSeq2.txt") or die "Can't read $!";  
	my $Rna = <$in>;
	$Rna = uc($Rna);
	chomp($Rna);
	my @seq = split('',$Rna);
	my @stems = ();
	@stems = Stems();
	Vertices();
	Consistance();
	## calculate energy of each stems
	my %stack_energy = stack_energy();
	calculate_Stem_Energy();
##Functions

# Create Dot Plot
sub DotPlot {
	my @matrix = ();
	for(my $i=0;$i<scalar(@seq);$i++)
	{
		for(my $j=0;$j<scalar(@seq);$j++)
		{
			if(($seq[$i] eq "G" and $seq[$j] eq "C") or ($seq[$i]eq"A" and $seq[$j]eq"U") or ($seq[$i]eq"C" and $seq[$j]eq"G") or ($seq[$i]eq"U" and $seq[$j]eq"A"))
			{
				$matrix[$i][$j]="\\";
			}
			else
			{
				$matrix[$i][$j]=" ";
			}
		}
		
	}
	return @matrix;
}
#End Create Dot Plot


# Declare Stems
sub Stems {
	my @matrix = DotPlot();
	# first declare stems with length 3
	# u1 : initial nucletide position
	# u2 : final nucletide position
	# u3 : length of the stem
	for(my $i=0;$i<scalar(@seq);++$i){
		for(my $j=scalar(@seq)-1;$j>=0;--$j){
			if($j<=$i){
				last;
			}
			if($matrix[$i][$j] eq "\\"){
				my $start=$i;
				my $end=$j;
				my $k;
				for($k=1;$k<scalar(@seq);++$k){
					if($matrix[$i+$k][$j-$k] eq " "){
						if(abs($i+$k-1-($j-$k+1)+1) < 3){
							$k=0;
						}
						last;
					}
					if($j-$k<=$i+$k){
						$k=0;
						last;
					}
				}
				my $length=$k;
				if($length>=3){
					push @stems, {u1=>$start, u2=>$end, u3=>$length };
				}
			}
		}
	}
	my @sorted =  sort { $a->{u3} <=> $b->{u3} } @stems;
	@stems = @sorted;
	return @stems;
}
#End Declare Stems

# Vertices
sub Vertices {
	open(my $fVertex,'>',"vertices.txt") or die "Can't open file for writing: $!";
	print $fVertex length($Rna),"\n";
	for(my $i=0;$i<scalar(@stems);$i++){
		print $fVertex ($stems[$i]{u1},",",$stems[$i]{u2},",",$stems[$i]{u3},"\n");
	}
	close $fVertex or die "Failed to close file: $!";
}
# End Vertices

# Consistance

sub Consistance {
	open(my $fEdge,'>',"edges.txt") or die "Can't open file for writing: $!";
	for(my $i=0;$i<scalar(@stems);$i++){
		for(my $j=0;$j<scalar(@stems);$j++){
			if( (($stems[$i]{u1}+$stems[$i]{u3}<=$stems[$j]{u1} and $stems[$j]{u2}<=$stems[$i]{u2}-$stems[$i]{u3} ) or ( $stems[$i]{u2}<$stems[$j]{u1} ) ) or 
				(($stems[$j]{u1}+$stems[$j]{u3}<=$stems[$i]{u1} and $stems[$i]{u2}<=$stems[$j]{u2}-$stems[$j]{u3} ) or ( $stems[$j]{u2}<$stems[$i]{u1} )) ){
				print $fEdge $i,",",$j,"\n";
			}
		}
	}
	close $fEdge or die "Failed to close file: $!";
}

# End Consistance

# Stack Energy

sub stack_energy {
	my %stack_energy;
	$stack_energy{"AA"}{"UU"} = -0.9;
	$stack_energy{"AC"}{"UG"} = -2.2;
	$stack_energy{"AG"}{"UC"} = -2.1;
	$stack_energy{"AG"}{"UU"} = -0.6;
	$stack_energy{"AU"}{"UA"} = -1.1;
	$stack_energy{"AU"}{"UG"} = -1.4;
	$stack_energy{"CA"}{"GU"} = -2.1;
	$stack_energy{"CC"}{"GG"} = -3.3;
	$stack_energy{"CG"}{"GC"} = -2.4;
	$stack_energy{"CG"}{"GU"} = -1.4;
	$stack_energy{"CU"}{"GA"} = -2.1;
	$stack_energy{"CU"}{"GG"} = -2.1;
	$stack_energy{"GA"}{"CU"} = -2.4;
	$stack_energy{"GA"}{"UU"} = -1.3;
	$stack_energy{"GC"}{"CG"} = -3.4;
	$stack_energy{"GC"}{"UG"} = -2.5;
	$stack_energy{"GG"}{"CC"} = -3.3;
	$stack_energy{"GG"}{"CU"} = -1.5;
	$stack_energy{"GG"}{"UC"} = -2.1;
	$stack_energy{"GG"}{"UU"} = -0.5;
	$stack_energy{"GU"}{"CA"} = -2.2;
	$stack_energy{"GU"}{"CG"} = -2.5;
	$stack_energy{"GU"}{"UA"} = -1.4;
	$stack_energy{"GU"}{"UG"} = 1.3;
	$stack_energy{"UA"}{"AU"} = -1.3;
	$stack_energy{"UA"}{"GU"} = -1;
	$stack_energy{"UC"}{"AG"} = -2.4;
	$stack_energy{"UC"}{"GG"} = -1.5;
	$stack_energy{"UG"}{"AC"} = -2.1;
	$stack_energy{"UG"}{"AU"} = -1;
	$stack_energy{"UG"}{"GC"} = -1.4;
	$stack_energy{"UG"}{"GU"} = 0.3;
	$stack_energy{"UU"}{"AA"} = -0.9;
	$stack_energy{"UU"}{"AG"} = -1.3;
	$stack_energy{"UU"}{"GA"} = -0.6;
	$stack_energy{"UU"}{"GG"} = -0.5;
	return %stack_energy;
}

# End Stack Energy

# Calculate stem energy

sub calculate_Stem_Energy {
	open(my $fEdge,'>>',"edges.txt") or die "Can't open file for writing: $!";
	for(my $i=0;$i<scalar(@stems);++$i){
		my $sum_energy = 0;
		for(my $j=0;$j<$stems[$i]{u3}-1;++$j){
			foreach my $top (sort { $a cmp $b} keys %stack_energy){
				foreach my $bottom ( keys %{$stack_energy{$top}}){
					my $pair1 = join('',$seq[$stems[$i]{u1}+$j],$seq[$stems[$i]{u1}+$j+1]);
					my $pair2 = join('',$seq[$stems[$i]{u2}-$j],$seq[$stems[$i]{u2}-$j-1]);
					if($pair1 eq $top and $pair2 eq $bottom)
					{
						$sum_energy = $sum_energy + $stack_energy{$top}{$bottom};
					}
				}
			}	
		}
		print $fEdge $i,",",$sum_energy,"\n";
	}
	close $fEdge or die "Failed to close file: $!";
}

# End Calculate stem energy

