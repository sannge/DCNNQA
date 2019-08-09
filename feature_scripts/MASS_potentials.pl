#!/usr/bin/env perl
use warnings;
use strict;

my $dir_MASS_potentials = $ARGV[0];
my $f_model = $ARGV[1];
my $f_ss = $ARGV[2];
my $f_model2 = $ARGV[3];

my $exe_vol = "$dir_MASS_potentials/calc-volume.exe";
# 1. PAP
my $f_PAP = "$f_model2.PAP";
`perl $dir_MASS_potentials/potentials/PAP/run_PAP_SS.pl $f_model $f_ss $f_PAP $dir_MASS_potentials`;
# 2. TAP
my $f_TAP = "$f_model2.TAP";
`perl $dir_MASS_potentials/potentials/TAP/run_TAP_SS.pl $f_model $f_ss $f_TAP $dir_MASS_potentials`;
# 3. CSP
# CA
my $f_CSP_CA = "$f_model2.CSP_CA";
`perl $dir_MASS_potentials/potentials/CSP/run_CSP_SS_CA.pl $f_model $f_ss $f_CSP_CA $dir_MASS_potentials`;
# CB
my $f_CSP_CB = "$f_model2.CSP_CB";
`perl $dir_MASS_potentials/potentials/CSP/run_CSP_SS_CB.pl $f_model $f_ss $f_CSP_CB $dir_MASS_potentials`;
# 4. ASPR
# CA
my $f_ASPR_CA = "$f_model2.ASPR_CA";
`perl $dir_MASS_potentials/potentials/ASPR/run_ASP_SS_CA.pl $f_model $f_ss $f_ASPR_CA $dir_MASS_potentials`;
# CB
my $f_ASPR_CB = "$f_model2.ASPR_CB";
`perl $dir_MASS_potentials/potentials/ASPR/run_ASP_SS_CB.pl $f_model $f_ss $f_ASPR_CB $dir_MASS_potentials`;
# 5. ASPA
my $f_ASPA = "$f_model2.ASPA";
`perl $dir_MASS_potentials/potentials/ASPA/run_ASPA.pl $f_model $f_ASPA`;
# 6. DDP
# CA
my $f_DDP_CA = "$f_model2.DDP_CA";
`perl $dir_MASS_potentials/potentials/DDP/run_DDP_CA.pl $f_model $f_DDP_CA $dir_MASS_potentials`;
# CB
my $f_DDP_CB = "$f_model2.DDP_CB";
`perl $dir_MASS_potentials/potentials/DDP/run_DDP_CB.pl $f_model $f_DDP_CB $dir_MASS_potentials`;
# 7. SSDP
# CA
my $f_SSDP_CA = "$f_model2.SSDP_CA";
`perl $dir_MASS_potentials/potentials/SSDP/run_SSDP_CA.pl $f_model $f_SSDP_CA $dir_MASS_potentials`;
# CB
my $f_SSDP_CB = "$f_model2.SSDP_CB";
`perl $dir_MASS_potentials/potentials/SSDP/run_SSDP_CB.pl $f_model $f_SSDP_CB $dir_MASS_potentials`;
# 8. CDP
# CA
my $f_CDP_CA = "$f_model2.CDP_CA";
`perl $dir_MASS_potentials/potentials/CDP/run_CDP_CA.pl $f_model $f_ss $f_CDP_CA $dir_MASS_potentials`;
# CB
my $f_CDP_CB = "$f_model2.CDP_CB";
`perl $dir_MASS_potentials/potentials/CDP/run_CDP_CB.pl $f_model $f_ss $f_CDP_CB $dir_MASS_potentials`;
# 9 RSAP
my $f_RSAP = "$f_model2.RSAP";
`perl $dir_MASS_potentials/potentials/RSAP/run_RSAP.pl $f_model $f_ss $f_RSAP $dir_MASS_potentials`;
# 10 VDP
my $f_vol = "$f_model2.volume";
`$exe_vol -i $f_model > $f_vol`;
my $f_VDP = "$f_model2.VDP_CA";
`perl $dir_MASS_potentials/potentials/VDP/run_VDP_CA.pl $f_model $f_ss $f_vol $f_VDP $dir_MASS_potentials`;
	


