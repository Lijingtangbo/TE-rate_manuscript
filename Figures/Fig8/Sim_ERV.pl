$pop_size = 1000;
$seed_prop = 1;
$nr_gen = 5000;
$max_transp_rate = 0.05; # effect of one competent erv on transposition prob of each erv (weighted by affinity) in the germline
$max_erv = 250;
$mutation_rate_per_erv = 0.0001; # probability that individual erv undergoes dnm in the germline (which may change comp and affin)
$affin_increase = "OFF";
## Pick ONE!
#$fitness_effect = "NONE";
#$fitness_effect = "LINEAR";
#$fitness_effect = "QUADRATIC";
$fitness_effect = "CUBIC";
$output_file_1 = ">Res_fitness_${fitness_effect}_affinity_${affin_increase}_mut_${mutation_rate_per_erv}.txt";
$output_file_2 = ">>Concat_res.txt";

if ($fitness_effect eq "NONE"){$expon = 0;}
elsif ($fitness_effect eq "LINEAR"){$expon = 1;}
elsif ($fitness_effect eq "QUADRATIC"){$expon = 2;}
elsif ($fitness_effect eq "CUBIC"){$expon = 3;};

print "Power is :",$expon,"\n";

# Functions:
sub comp_transp_rate {
	my ($nr_comp) = @_;
	my $transp_rate = (1-exp(-0.1*$nr_comp))*$max_transp_rate;
	return $transp_rate;
}

sub comp_fitness {
	my ($nr_erv,$power) = @_;
	if ($power == 0){$fitness = 1;}
	else {$fitness = 1 - $nr_erv**$power/$max_erv**$power;};
	return $fitness;
}

sub z_rand {
	my ($u1, $u2);  # uniformly distributed random numbers
	my $w;          # variance, then a weight
	my ($g1, $g2);  # gaussian-distributed numbers

	do {
		$u1 = 2 * rand() - 1;
		$u2 = 2 * rand() - 1;
		$w = $u1*$u1 + $u2*$u2;
	} while ( $w >= 1 );

	$w = sqrt( (-2 * log($w))  / $w );
	$g2 = $u1 * $w;
	$g1 = $u2 * $w;
	return $g1;
}

# Generation 1:
for ($i=1;$i<=$pop_size*$seed_prop;$i++)
{ 
	$nr_erv[$i] = 1;
	$compet_erv[1][$i] = 1;
	$nr_compet_erv[$i] = 1;
	$affin_erv[1][$i] = 1;
	$ident_erv[1][$i] = 1;
	$fitness[$i] = 1;  
};
for ($i=$pop_size*$seed_prop+1;$i<=$pop_size;$i++)
{
	$nr_erv[$i] = 0;
	$fitness[$i] = 1; 
};
$nr_erv_ident = 1;
print "Gen 1:	";
for ($j=1;$j<=$pop_size;$j++)
{
	print "	",$nr_erv[$j];
};
print "\n";
print "Gen 1:	";
for ($j=1;$j<=$pop_size;$j++)
{
	print "	",$nr_compet_erv[$j];
};
print "\n";
print "Gen 1:	";
for ($j=1;$j<=$pop_size;$j++)
{
	print "	",$fitness[$j];
};
print "\n\n";

# Generation 2 to $nr_gen
for ($i=2;$i<=$nr_gen;$i++)
{
	# Generate  parents = copy of present generation
	for ($j=1;$j<=$pop_size;$j++)
	{
		$par_nr_erv[$j] = $nr_erv[$j]; 
		$par_fitness[$j] = $fitness[$j];	
		for ($k=1;$k<=$par_nr_erv[$j];$k++)
		{
			$par_compet_erv[$k][$j] = $compet_erv[$k][$j]; 
			$par_affin_erv[$k][$j] = $affin_erv[$k][$j];
			$par_ident_erv[$k][$j] = $ident_erv[$k][$j];
		};	
	};
	#print "Parents generated\n";
	# Generate new offspring by mating parents 
	for ($j=1;$j<=$pop_size;$j++)
	{
		#print "Working on offspring ",$j,"\n";
		$nr_erv[$j] = 0;
		# Select a father
		$fath_sel = 0;
		while ($fath_sel == 0)
		{
			$fath = int(rand($pop_size)+1);
			$toss_sel = rand();
			if ($toss_sel < $par_fitness[$fath])
			{
				$fath_sel = 1;
			};	
		};
		#print "Father: ",$fath,"\n";
		# Count the number of competent ERVs in the father's genome and the sum of the ERV's affinities
		$nr_compet = 0;
		$sum_affin = 0;
		for ($k=1;$k<=$par_nr_erv[$fath];$k++)
		{
			if ($par_compet_erv[$k][$fath] == 1)
			{
				$nr_compet += 1;
			};
			$sum_affin  += $par_affin_erv[$k][$fath];
		};
		#print "Number of competent ERV: ",$nr_compet,"\n";
		# Compute transposition rate in father
		$transp_rate = comp_transp_rate($nr_compet);
		#print "Father ",$j," - transposition rate ",$transp_rate,"\n";
		#print "Transposition rate: ",$transp_rate,"\n";
		# Mendelian sampling in the father
		for ($k=1;$k<=$par_nr_erv[$fath];$k++)
		{
			$toss_mend = rand();
			if ($toss_mend > 0.5)
			{
				$nr_erv[$j] += 1;
				$compet_erv[$nr_erv[$j]][$j] = $par_compet_erv[$k][$fath];
				$affin_erv[$nr_erv[$j]][$j] = $par_affin_erv[$k][$fath]; 	
				$ident_erv[$nr_erv[$j]][$j] = $par_ident_erv[$k][$fath];
				#print "ERV nr ",$k," transmitted\n";
			};
			
		# Generate dn transposition in father's gamete
			if ($sum_affin > 0)
			{
				$erv_prob_transp = $transp_rate * ($par_affin_erv[$k][$fath]/$sum_affin);
			}
			else
			{
				$erv_prob_transp = 0;
			};
			
			#print "	ERV ",$k," - probability of transposition ",$erv_prob_transp,"\n";
			$toss_dntr = rand();
			if ($toss_dntr < $erv_prob_transp)
			{
				$nr_erv[$j] += 1;
				$compet_erv[$nr_erv[$j]][$j] = $par_compet_erv[$k][$fath];
				$affin_erv[$nr_erv[$j]][$j] = $par_affin_erv[$k][$fath]; 	
				$ident_erv[$nr_erv[$j]][$j] = $par_ident_erv[$k][$fath];
				#print "ERV nr ",$k," transposed\n";
			};	
		};
		# Select a mother
		$moth_sel = 0;
		while ($moth_sel == 0)
		{
			$moth = int(rand($pop_size)+1);
			$toss_sel = rand();
			if ($toss_sel < $par_fitness[$moth])
			{
				$moth_sel = 1;	
			};	
		};
		# Count the number of competent ERVs in the mother's genome and the sum of the ERV's affinities
		$nr_compet = 0;
		$sum_affin = 0;
		for ($k=1;$k<=$par_nr_erv[$moth];$k++)
		{
			if ($par_compet_erv[$k][$moth] == 1)
			{
				$nr_compet += 1;
			};
			$sum_affin += $par_affin_erv[$k][$moth];
		};
		# Compute transposition rate in the mother
		$transp_rate = comp_transp_rate($nr_compet);
		# Mendelian sampling in the mother
		for ($k=1;$k<=$par_nr_erv[$moth];$k++)
		{
			$toss_mend = rand();
			if ($toss_mend > 0.5)
			{
				$nr_erv[$j] += 1;
				$compet_erv[$nr_erv[$j]][$j] = $par_compet_erv[$k][$moth];
				$affin_erv[$nr_erv[$j]][$j] = $par_affin_erv[$k][$moth]; 	
				$ident_erv[$nr_erv[$j]][$j] = $par_ident_erv[$k][$moth];
			};
		# Generate dn transposition in mother's gamete
			if ($sum_affin > 0)
			{
				$erv_prob_transp = $transp_rate * ($par_affin_erv[$k][$moth]/$sum_affin);
			}
			else
			{
				$erv_prob_transp = 0;
			};
			$toss_dntr = rand();
			if ($toss_dntr <= $erv_prob_transp)
			{
				$nr_erv[$j] += 1;
				$compet_erv[$nr_erv[$j]][$j] = $par_compet_erv[$k][$moth];
				$affin_erv[$nr_erv[$j]][$j] = $par_affin_erv[$k][$moth]; 	
				$ident_erv[$nr_erv[$j]][$j] = $par_ident_erv[$k][$moth];
			};		
		};
		# compute fitness of new offspring
		$fitness[$j] = comp_fitness($nr_erv[$j],$expon);
		if ($fitness[$j] < 0)
		{
			$fitness[$j] = 0;
		}; 
		# Mutate the ERVs in the offspring 
		for ($k=1;$k<=$nr_erv[$j];$k++)
		{
			$toss_dnm = rand();
			if ($toss_dnm <= $mutation_rate_per_erv)
			{
				$nr_erv_ident += 1;
				$ident_erv[$k][$j] = $nr_erv_ident;	
				
				$compet_erv[$k][$j] = 0;
				
				if ($affin_increase eq "ON")
				{ 
					if ($affin_erv[$k][$j] != 0)
					{
						$affin_erv[$k][$j] += z_rand();
						if ($affin_erv[$k][$j] < 0)
						{
							$affin_erv[$k][$j] = 0;
						};
					};
				}
				elsif ($affin_increase eq "OFF")
				{ 
					$affin_erv[$k][$j] = 0;
				};
			};	
		};
		# Count number of competent ERVs
		$nr_compet_erv[$j] = 0;
		for ($k=1;$k<=$nr_erv[$j];$k++)
		{
			if ($compet_erv[$k][$j] == 1)
			{
				$nr_compet_erv[$j] += 1;
			};
		};
	};
	$average_nr_erv[$i] = 0;
	$average_nr_compet_erv[$i] = 0;
	$average_nr_high_aff_erv[$i] = 0;
	$average_fitness[$i] = 0;
	for ($j=1;$j<=$pop_size;$j++)
	{
		$average_nr_erv[$i] += $nr_erv[$j];
		$average_nr_compet_erv[$i] += $nr_compet_erv[$j];
		for ($k=1;$k<=$nr_erv[$j];$k++)
		{
			if ($affin_erv[$k][$j] > 1)
			{
				$average_nr_high_aff_erv[$i] += 1;
			};	
		};
		$average_fitness[$i] += $fitness[$j];	
	};
	$average_nr_erv[$i] = $average_nr_erv[$i]/$pop_size;
	$average_nr_compet_erv[$i] = $average_nr_compet_erv[$i]/$pop_size;
	$average_nr_high_aff_erv[$i] = $average_nr_high_aff_erv[$i]/$pop_size;
	$average_nr_def_erv[$i] = $average_nr_erv[$i] - $average_nr_compet_erv[$i] - $average_nr_high_aff_erv[$i];
	$average_fitness[$i] = $average_fitness[$i]/$pop_size;
	print "Average fitness: ",$average_fitness[$i],"\n"; 
	if ($average_fitness[$i] == 0)
	{
		$i = $nr_gen+1;	
	};
	
	print "Gen ",$i,":	";
	for ($j=1;$j<=$pop_size;$j++)
	{
		print "	",$nr_erv[$j];
	};
	print "\n";
	print "Gen ",$i,":	";
	for ($j=1;$j<=$pop_size;$j++)
	{
		print "	",$nr_compet_erv[$j];
	};
	print "\n";
	print "Gen ",$i,":	";
	for ($j=1;$j<=$pop_size;$j++)
	{
		print "	",$fitness[$j];
	};
	print "\n\n";
};
open RES, ">results.txt";
print RES "GENER	AV_NR_ERV	AV_NR_COMP_ERV	AV_NR_HIGH_AFF_ERV	AV_NR_DEF_ERV	AV_FITNESS	MUT_RATE	FITNESS\n";
for ($i=2;$i<=$nr_gen;$i++)
{
	print RES $i,"	",$average_nr_erv[$i],"	",$average_nr_compet_erv[$i],"	",$average_nr_high_aff_erv[$i],"	",$average_nr_def_erv[$i],"	",$average_fitness[$i];
	print RES "	",$mutation_rate_per_erv,"	",$fitness_effect,"\n";
};
close RES;

open RES2, $output_file_1;
print RES2 "GENER	ERV	NR	MUT_RATE	FITNESS\n";
for ($i=2;$i<=$nr_gen;$i++)
{
	print RES2 $i,"	ALL	",$average_nr_erv[$i],"	",$mutation_rate_per_erv,"	",$fitness_effect,"\n";
	print RES2 $i,"	COMPETENT	",$average_nr_compet_erv[$i],"	",$mutation_rate_per_erv,"	",$fitness_effect,"\n";
	print RES2 $i,"	HIGH_AFFINITY	",$average_nr_high_aff_erv[$i],"	",$mutation_rate_per_erv,"	",$fitness_effect,"\n";
	print RES2 $i,"	DEFECTIVE	",$average_nr_def_erv[$i],"	",$mutation_rate_per_erv,"	",$fitness_effect,"\n";
};
close RES2;

open RES3, $output_file_2;
#print RES3 "GENER	ERV	NR	MUT_RATE	FITNESS\n";
for ($i=2;$i<=$nr_gen;$i++)
{
	print RES3 $i,"	ALL	",$average_nr_erv[$i],"	",$mutation_rate_per_erv,"	",$fitness_effect,"\n";
	print RES3 $i,"	COMPETENT	",$average_nr_compet_erv[$i],"	",$mutation_rate_per_erv,"	",$fitness_effect,"\n";
	print RES3 $i,"	HIGH_AFFINITY	",$average_nr_high_aff_erv[$i],"	",$mutation_rate_per_erv,"	",$fitness_effect,"\n";
	print RES3 $i,"	DEFECTIVE	",$average_nr_def_erv[$i],"	",$mutation_rate_per_erv,"	",$fitness_effect,"\n";
};
close RES2;
