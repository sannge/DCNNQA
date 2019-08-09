#################################################################################
# parse the stride output file:
# 1, extract the secondary structure: Amino acid code
# 2, calculate the relative accessible surface area at the 25% exposure threshold
#    Relacc = Acc/Maxacc;
#    Relacc <  25% => "buried: -"
#    Relacc >= 25% => "exposed: e"
#
# Secondary structure Code:
# DSSP has 8 classes(types): 
#   H = alpha-helix
#   B = residue in isolated beta-bridge
#   E = extended strand, participates in beta-ladder
#   G = 3-helix (310 helix)
#   I = 5 helix (pi-helix)
#   T = hydrogen bonded turn
#   S = bend
#   C = ciol
# Whereas Stride has 7 classes:
#   H	    Alpha helix
#   G	    3-10 helix
#   I	    PI-helix
#   E	    Extended conformation
#   B or	b   Isolated bridge
#   T	    Turn
#   C	    Coil (none of the above)
#
# However, eight types of secondary structures are too many for the existing
# methods of secondary structure prediction. Instead usually only three states
# are predicted: helix(H), extended(beta-sheet)(E) and coil(C). There are many
# different methods to translate the eight-letter DSSP alphabet into the 
# three-letter code. The translation used in the CASP experiment is as follows:
#   H,G,I -> H
#   E,B   -> E
#   T,S,C -> C
#
# Stride Output format
# ASG    Detailed secondary structure assignment
#    Format:  6-8  Residue name
#             10-10 Protein chain identifier
#	      12-15 PDB	residue	number
#	      17-20 Ordinal residue number
#	      25-25 One	letter secondary structure code
#	      27-39 Full secondary structure name
#	      43-49 Phi	angle
#	      53-59 Psi	angle
#	      65-69 Residue solvent accessible area
#
# Tong Liu
# 10/10/2015
#################################################################################

sub parse_stride_ss_angle{
  my ($stride_file, $stride) = @_;

  my @aacode = ('A', 'B', 'C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X', 'Y', 'Z');
  my @Maxacc = (106, 160, 135, 163, 194, 197, 84, 184, 169, 205, 164, 188, 157, 136, 198, 248, 130, 142, 142, 227, 180, 222, 196);

  my %Maxacc;
  $Maxacc{'A'} = 106; $Maxacc{'B'} = 160; $Maxacc{'C'} = 135;
  $Maxacc{'D'} = 163; $Maxacc{'E'} = 194; $Maxacc{'F'} = 197;
  $Maxacc{'G'} = 84;  $Maxacc{'H'} = 184; $Maxacc{'I'} = 169;
  $Maxacc{'K'} = 205; $Maxacc{'L'} = 164; $Maxacc{'M'} = 188;
  $Maxacc{'N'} = 157; $Maxacc{'P'} = 136; $Maxacc{'Q'} = 198;
  $Maxacc{'R'} = 248; $Maxacc{'S'} = 130; $Maxacc{'T'} = 142;
  $Maxacc{'V'} = 142; $Maxacc{'W'} = 227; $Maxacc{'X'} = 180;
  $Maxacc{'Y'} = 222; $Maxacc{'Z'} = 196;

  my %amino=();
  $amino{"ALA"} = 'A'; $amino{"CYS"} = 'C';
  $amino{"ASP"} = 'D'; $amino{"GLU"} = 'E';
  $amino{"PHE"} = 'F'; $amino{"GLY"} = 'G';
  $amino{"HIS"} = 'H'; $amino{"ILE"} = 'I';
  $amino{"LYS"} = 'K'; $amino{"LEU"} = 'L';
  $amino{"MET"} = 'M'; $amino{"ASN"} = 'N';
  $amino{"PRO"} = 'P'; $amino{"GLN"} = 'Q';
  $amino{"ARG"} = 'R'; $amino{"SER"} = 'S';
  $amino{"THR"} = 'T'; $amino{"VAL"} = 'V';
  $amino{"TRP"} = 'W'; $amino{"TYR"} = 'Y'; 

  open IN, "$stride_file" or die $!;
  while(my $line = <IN>){
    chomp $line;
    if(substr($line, 0, 4) eq "ASG " ){
        my @items = split(/\s+/, $line);

        my $ordinal_number = $items[4];
        my $residue_name = $items[1];
        my $pdb_residue_number = $items[3];
        my $ss_code_7 = $items[5];
        my $acc = $items[9];
        my $ang1 = $items[7];
        my $ang2 = $items[8];

        $residue_name = uc $residue_name;
        $ss_code_7 = uc $ss_code_7;

        my $ss_code_3;
        if($ss_code_7 eq "H" || $ss_code_7 eq "G" || $ss_code_7 eq "I"){
          $ss_code_3 = "H";
        }
        if($ss_code_7 eq "E" || $ss_code_7 eq "B"){
          $ss_code_3 = "E";
        }
        if($ss_code_7 eq "T" || $ss_code_7 eq "S" || $ss_code_7 eq "C"){
          $ss_code_3 = "C";
        }
        push(@{$stride->{$ordinal_number}}, $pdb_residue_number);
        push(@{$stride->{$ordinal_number}}, $residue_name);
        push(@{$stride->{$ordinal_number}}, $ss_code_3);
        push(@{$stride->{$ordinal_number}}, $ang1);
        push(@{$stride->{$ordinal_number}}, $ang2);
        #$stride->{$pdb_residue_number} = $ss_code_3;
    }
  }
  close IN;
}


1;


# End
