use strict;
use Data::Dumper;
use Getopt::Long;

my $G_USAGE = "
$0 <command> -h
<command> [sievegraph]
tablebam    : bam to table of alignments
breakpoints : table of breakpoints
tagbam      : annotate bam with tags

      Copyright (C) 2019-2024 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
      Please contact author for license (cheehongsg at gmail.com)
";

my $command = undef;
if (!defined $ARGV[0] || substr($ARGV[0],0,1) eq '-' || substr($ARGV[0],0,2) eq '--') {
	die("Please specify the command.\n",$G_USAGE);
}
$command = shift @ARGV; $command = lc $command;

# auto-flush for both stderr and stdout
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if ('tablebam' eq $command) {
	tableBam();
} elsif ('breakpoints' eq $command) {
	breakpoints();
} elsif ('tagbam' eq $command) {
	tagBam();
} else {
}

exit 0;

sub getReadLenFromCigar {
  my ($cigar) = @_;
  my $readLen = 0;
  my @ops = split(/([MIDNSHP=X])/, $cigar);
  for(my $i=0; $i<scalar(@ops); ++$i) {
    my $count = int($ops[$i]);
    $i++;
    my $op = $ops[$i];
    if ('M' eq $op || '=' eq $op || 'X' eq $op) {
      $readLen += $count;
    } elsif ('I' eq $op || 'S' eq $op) {
      $readLen += $count;
    } elsif ('D' eq $op || 'N' eq $op) {
    } elsif ('H' eq $op) {
      $readLen += $count;
    }
  }
  return $readLen;
}

sub cigarToCoordinates {
  my ($cigar, $reverse, $countH) = @_;
  # TODO: countP?
  $countH = 0 if (!defined $countH);
  $countH = 1 if (0!=int($countH));

  my $read_pos = 0;
  my $ref_pos = 0;
  my $read_start0 = -1;
  my $read_end1 = -1;
  my $ref_start0 = -1;
  my $ref_end1 = -1;
  my @ops = split(/([MIDNSHP=X])/, $cigar);
  if (0!=$reverse) {
    for(my $i=0; $i<scalar(@ops); $i+=2) {
      my $tmp = $ops[$i];
      $ops[$i] = $ops[$i+1];
      $ops[$i+1] = $tmp;
    }
    @ops = reverse @ops;
  }
  for(my $i=0; $i<scalar(@ops); ++$i) {
    my $count = int($ops[$i]);
    $i++;
    my $op = $ops[$i];
    if ('M' eq $op || '=' eq $op || 'X' eq $op) {
      $read_pos += $count;
      $ref_pos += $count;
    } elsif ('I' eq $op || 'S' eq $op) {
      $read_pos += $count;
    } elsif ('D' eq $op || 'N' eq $op) {
      $ref_pos += $count;
    } elsif ('H' eq $op && 0!=$countH) {
      $read_pos += $count;
    }
    if ('M' eq $op) {
      if ($read_start0 == -1) {
        $read_start0 = $read_pos - $count;
      }
      $read_end1 = $read_pos;
      if ($ref_start0 == -1) {
        $ref_start0 = $ref_pos - $count;
      }
      $ref_end1 = $ref_pos;
    }
  }

  return ($read_start0, $read_end1, $ref_start0, $ref_end1);
}

sub getSourceReadId {
  my ($readId) = @_;
  my $effReadId = $readId;
  $effReadId =~ s/\/rc$//;
  return $effReadId;
}

sub tableBam {
  my $G_USAGE = "
$0 tablebam --bam <bamFile> --easycoord --bed
  --bam STR           bam file
  --exclude STR       query to exclude (multiple)
  --easycoord         pick representative alignment based on coordinates
  --bed               generate alignment bed file
";

  my $bamFile = undef;
  my $easyCoord = 0;
  my $useFirstCG = 0;
  my $bedFormat = 0;
  my @excludeIds = ();
  my $verbose = 0;
	my $help = 0;
  ## my $debugEdge = "";
	
	GetOptions (
	"bam=s"        => \$bamFile,
	"exclude=s"    => \@excludeIds,
	"easycoord!"   => \$easyCoord,
	"usefirstcg!"  => \$useFirstCG,
	"bed!"         => \$bedFormat,
	"verbose!"     => \$verbose,
	"help!"        => \$help)
	or die("Error in command line arguments\n$G_USAGE");
	
	die "$G_USAGE" if ($help);

  open INFILE, "samtools view -h -F 0x700 $bamFile | " || die "Fail to open $bamFile\n$!\n";
  my %excludeIds = ();
  grep { $excludeIds{$_} = 1; } @excludeIds;
  my @readids = ();
  my %reads = ();
  my %refLens = ();
  my $numRefs = 0;

  while (<INFILE>) {
    if (/^@/) {
      if (/^\@SQ\t/) {
        my @bits = split(/\t/);
        my $sqsn = '';
        my $sqln = 0;
        foreach my $bit (@bits) {
          my ($tag, $value) = split(/:/, $bit, 2);
          if ('SN' eq $tag) {
            $sqsn = $value;
          } elsif ('LN' eq $tag) {
            $sqln = int($value);
          }
        }
        if ('' ne $sqsn) {
          $refLens{$sqsn} = $sqln;
          $numRefs++;
        }
      }
      next;
    }
    chomp();
    my @bits = split(/\t/);
    next if (exists $excludeIds{$bits[0]});
    if (!exists $reads{$bits[0]}) {
      push @readids, $bits[0];
      $reads{$bits[0]} = {id=>$bits[0],records=>[]};
    }
    my %record = (read=>$bits[0],flag=>$bits[1],contig=>$bits[2],start=>$bits[3],mapq=>$bits[4],cigar=>$bits[5]);
    push @{$reads{$bits[0]}->{records}}, \%record;

    # calculate the coordinates
    my ($read_start0, $read_end1, $ref_start0, $ref_end1) = cigarToCoordinates($record{cigar}, ($record{flag} & 0x10), 1);
    $record{read_start0} = $read_start0;
    $record{read_end1} = $read_end1;
    $record{ref_start0} = $ref_start0;
    $record{ref_end1} = $ref_end1;
    #printf "%s\t%d\t%s\t%d\tread:%d-%d[%d]\tref:%d-%d[%s%d]\t%s\n",
    #  $record{read}, $record{flag}, $record{contig}, $record{start},
    #  $read_start0+1, $read_end1, $read_end1-$read_start0, $ref_start0+1, $ref_end1, ($record{flag} & 0x10)?'-':'+', $ref_end1-$ref_start0,
    #  $record{cigar};
    # TODO: strand adjusted
  }
  close INFILE;

  @readids = sort @readids;
  if (scalar(@readids)>0) {
    if (0==$bedFormat) {
      my @cols = ('readId','flag','mapq','contig','contigStart','numAlignments','alignId','readLen','refLen','overlapNote','readStartEndLen','refStartEndLen','cigar');
      print join("\t", @cols), "\n";
    }
  }
  my %selectedReadIds = ();
  if (0!=$easyCoord) {
    # representative, prefer non-RC version
    my %selections = ();
    foreach my $readId (@readids) {
      my $effReadId = getSourceReadId($readId);
      if (!exists $selections{$effReadId}) {
        my %item = ();
        $selections{$effReadId} = \%item;
      }
      my $readRef = $reads{$readId};
      if ($effReadId eq $readId) {
        $selections{$effReadId}->{srcRef} = $readRef;
      } else {
        $selections{$effReadId}->{rcRef} = $readRef;
      }
    }
    if (0!=$useFirstCG) {
      # use the first CG's selected strand for all other samples!
      my $readId = $readids[0];
      my $alignmentsRef = $selections{$readId};
      my $selectStrandKey = 'srcRef';
      if (exists $alignmentsRef->{srcRef}) {
        if (exists $alignmentsRef->{rcRef}) {
          my $srcRef = $alignmentsRef->{srcRef};
          my @srcRecords = sort { $a->{read_start0}<=>$b->{read_start0} || $a->{read_end1}<=>$b->{read_end1} } @{$srcRef->{records}};
          my $rcRef = $alignmentsRef->{rcRef};
          my @rcRecords = sort { $a->{read_start0}<=>$b->{read_start0} || $a->{read_end1}<=>$b->{read_end1} } @{$rcRef->{records}};
          if (scalar(@srcRecords)>0) {
            if (scalar(@rcRecords)>0) {
              my $srcRecordRef = $srcRecords[0];
              my $rcRecordRef = $rcRecords[0];
              if (0==($srcRecordRef->{flag} & 0x10)) {
                $selectStrandKey = 'srcRef';
              } elsif (0==($rcRecordRef->{flag} & 0x10)) {
                $selectStrandKey = 'rcRef';
              } else {
                $selectStrandKey = 'srcRef';
              }
            } else {
              $selectStrandKey = 'srcRef';
            }
          } else {
            if (scalar(@rcRecords)>0) {
              $selectStrandKey = 'rcRef';
            } else {
              die "Cannot determine the first CG prefer strand!\n";
            }
          }
        } else {
          # no choice!
          $selectStrandKey = 'srcRef';
        }
      } else {
        if (exists $alignmentsRef->{rcRef}) {
          # no choice!
          $selectStrandKey = 'rcRef';
        } else {
          die "Cannot determine the first CG prefer strand!\n";
        }
      }
      # let's select the same strand for all other CGs
      my $complementStrandKey = ($selectStrandKey eq 'srcRef') ? 'rcRef' : 'srcRef';
      while (my ($readId, $alignmentsRef) = each %selections) {
        if (exists $alignmentsRef->{$selectStrandKey}) {
          my $strandRef = $alignmentsRef->{$selectStrandKey};
          $selectedReadIds{$strandRef->{id}} = $strandRef->{id};
        } elsif (exists $alignmentsRef->{$complementStrandKey}) {
          my $strandRef = $alignmentsRef->{$complementStrandKey};
          $selectedReadIds{$strandRef->{id}} = $strandRef->{id};
        } else {
          # should not happen
          die "No CG?\n";
        }
      }
    } else {
      while (my ($readId, $alignmentsRef) = each %selections) {
        if (exists $alignmentsRef->{srcRef}) {
          if (exists $alignmentsRef->{rcRef}) {
            my $srcRef = $alignmentsRef->{srcRef};
            my @srcRecords = sort { $a->{read_start0}<=>$b->{read_start0} || $a->{read_end1}<=>$b->{read_end1} } @{$srcRef->{records}};
            my $rcRef = $alignmentsRef->{rcRef};
            my @rcRecords = sort { $a->{read_start0}<=>$b->{read_start0} || $a->{read_end1}<=>$b->{read_end1} } @{$rcRef->{records}};
            if (scalar(@srcRecords)>0) {
              if (scalar(@rcRecords)>0) {
                my $srcRecordRef = $srcRecords[0];
                my $rcRecordRef = $rcRecords[0];
                if (0==($srcRecordRef->{flag} & 0x10)) {
                  $selectedReadIds{$srcRef->{id}} = $srcRef->{id};
                } elsif (0==($rcRecordRef->{flag} & 0x10)) {
                  $selectedReadIds{$rcRef->{id}} = $rcRef->{id};
                } else {
                  $selectedReadIds{$srcRef->{id}} = $srcRef->{id};
                }
              } else {
                $selectedReadIds{$srcRef->{id}} = $srcRef->{id};
              }
            } else {
              if (scalar(@rcRecords)>0) {
                $selectedReadIds{$rcRef->{id}} = $rcRef->{id};
              } else {
                # TODO: should not happen!
              }
            }
          } else {
            # no choice!
            $selectedReadIds{$alignmentsRef->{srcRef}->{id}} = $alignmentsRef->{srcRef}->{id};
          }
        } else {
          if (exists $alignmentsRef->{rcRef}) {
            # no choice!
            $selectedReadIds{$alignmentsRef->{rcRef}->{id}} = $alignmentsRef->{rcRef}->{id};
          } else {
            # TODO: should not happen!
          }
        }
      }
    }
  } else {
    # all
    foreach my $readId (@readids) {
      $selectedReadIds{$readId} = $readId;
    }
  }

  foreach my $readId (@readids) {
    my $readRef = $reads{$readId};
    my @records = sort { $a->{read_start0}<=>$b->{read_start0} || $a->{read_end1}<=>$b->{read_end1} } @{$readRef->{records}};
    my $numAlignments = scalar(@records);
    my $alignId = 0;
    my $readLen = 0;
    my $refLen = 0;
    my $numRecords = scalar(@records);
    my $prevRecRef = undef;
    my $nextRecRef = undef;
    for(my $i=0; $i<$numRecords; ++$i) {
      my $recordRef = $records[$i];
      if (0!=$easyCoord && 0==$i) {
        last if (!exists $selectedReadIds{$readId});
        # last if ($recordRef->{flag} & 0x10);
      }
      $prevRecRef = $records[$i-1] if ($i>0);
      $nextRecRef = $records[$i+1] if ($i<($numRecords+1));
      $alignId++;
      $refLen = (exists $refLens{$recordRef->{contig}}) ? $refLens{$recordRef->{contig}} : 0;
      if (0==$readLen) {
        $readLen = getReadLenFromCigar($recordRef->{cigar});
      }
      my $overlapStatus = 'n.a.';
      if (defined $prevRecRef) {
        if ($prevRecRef->{read_end1}>=($recordRef->{read_start0}+1)) {
          my $obp = $prevRecRef->{read_end1} - $recordRef->{read_start0};
          $overlapStatus = 'start('.$obp.')';
        }
      }
      if (defined $nextRecRef) {
        if ($recordRef->{read_end1}>=($nextRecRef->{read_start0}+1)) {
          my $obp = $recordRef->{read_end1} - $nextRecRef->{read_start0};
          if ($overlapStatus ne 'n.a.') {
            $overlapStatus .= ',end('.$obp.')';
          } else {
            $overlapStatus = 'end('.$obp.')';
          }
        }
      }

      if (0==$bedFormat) {
        printf "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\tread:%d-%d[%d]\t%s:%d-%d[%s%d]\t%s\n",
          $recordRef->{read}, $recordRef->{flag}, $recordRef->{mapq}, $recordRef->{contig}, $recordRef->{start},
          $numAlignments, $alignId, $readLen, $refLen,
          $overlapStatus,
          $recordRef->{read_start0}+1, $recordRef->{read_end1}, $recordRef->{read_end1}-$recordRef->{read_start0}, 
          $recordRef->{contig}, $recordRef->{start}+$recordRef->{ref_start0}, $recordRef->{start}+$recordRef->{ref_end1}-1, ($recordRef->{flag} & 0x10)?'-':'+', $recordRef->{ref_end1}-$recordRef->{ref_start0},
          $recordRef->{cigar};
      } else {
        my @cols = ($recordRef->{contig}, $recordRef->{start}+$recordRef->{ref_start0}-1, $recordRef->{start}+$recordRef->{ref_end1}-1);
        push @cols, sprintf("%s_%dof%d:%d-%d[%d]", $recordRef->{read}, $alignId, $numAlignments, $recordRef->{read_start0}+1, $recordRef->{read_end1}, $recordRef->{read_end1}-$recordRef->{read_start0});
        push @cols, $recordRef->{mapq}*10;
        push @cols, ($recordRef->{flag} & 0x10)?'-':'+';
        print join("\t", @cols), "\n";
      }
    }

    ### foreach my $recordRef (@records) {
    ###   $alignId++;
    ###   $refLen = (exists $refLens{$recordRef->{contig}}) ? $refLens{$recordRef->{contig}} : 0;
    ###   if (0==$readLen) {
    ###     $readLen = getReadLenFromCigar($recordRef->{cigar});
    ###   }
    ###   printf "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\tread:%d-%d[%d]\tref:%d-%d[%s%d]\t%s\n",
    ###     $recordRef->{read}, $recordRef->{flag}, $recordRef->{mapq}, $recordRef->{contig}, $recordRef->{start},
    ###     $numAlignments, $alignId, $readLen, $refLen,
    ###     $recordRef->{read_start0}+1, $recordRef->{read_end1}, $recordRef->{read_end1}-$recordRef->{read_start0}, 
    ###     $recordRef->{start}+$recordRef->{ref_start0}, $recordRef->{start}+$recordRef->{ref_end1}-1, ($recordRef->{flag} & 0x10)?'-':'+', $recordRef->{ref_end1}-$recordRef->{ref_start0},
    ###     $recordRef->{cigar};
    ### }
  }
}

sub getIndelsCoordinateFromCigar {
  my ($cigar, $reverse, $countH, $indelsRef, $minDel, $minIns) = @_;
  # TODO: countP?
  $countH = 0 if (!defined $countH);
  $countH = 1 if (0!=int($countH));

  $minDel = 50 if (!defined $minDel);
  $minIns = 50 if (!defined $minIns);

  @{$indelsRef} = ();

  my $read_pos = 0;
  my $ref_pos = 0;
  my $read_start0 = -1;
  my $read_end1 = -1;
  my $ref_start0 = -1;
  my $ref_end1 = -1;
  my @ops = split(/([MIDNSHP=X])/, $cigar);
  if (0!=$reverse) {
    for(my $i=0; $i<scalar(@ops); $i+=2) {
      my $tmp = $ops[$i];
      $ops[$i] = $ops[$i+1];
      $ops[$i+1] = $tmp;
    }
    @ops = reverse @ops;
  }
  for(my $i=0; $i<scalar(@ops); ++$i) {
    my $count = int($ops[$i]);
    $i++;
    my $op = $ops[$i];
    if ('M' eq $op || '=' eq $op || 'X' eq $op) {
      $read_pos += $count;
      $ref_pos += $count;
    } elsif ('I' eq $op) {
      if ($count>=$minIns) {
        $read_pos += $count;
        my %item = (type=>'i', start0=>$ref_pos, end=>$ref_pos);
        push @{$indelsRef}, \%item;
      } else {
        $read_pos += $count;
      }
    } elsif ('S' eq $op) {
      $read_pos += $count;
    } elsif ('D' eq $op) {
      if ($count>=$minDel) {
        my %item = (type=>'d', start0=>$ref_pos);
        $ref_pos += $count;
        $item{end}=$ref_pos;
        push @{$indelsRef}, \%item;
      } else {
        $ref_pos += $count;
      }
    } elsif ('N' eq $op) {
      $ref_pos += $count;
    } elsif ('H' eq $op && 0!=$countH) {
      $read_pos += $count;
    }
    if ('M' eq $op) {
      if ($read_start0 == -1) {
        $read_start0 = $read_pos - $count;
      }
      $read_end1 = $read_pos;
      if ($ref_start0 == -1) {
        $ref_start0 = $ref_pos - $count;
      }
      $ref_end1 = $ref_pos;
    }
  }
}

sub breakpoints {
  my $G_USAGE = "
$0 breakpoints --alignment <alignmentFile> --bed --mindellen <min_del_len> --mininslen <min_ins_len>
  --alignment STR     alignments file
  --mindellen INT     min deletion length
  --mininslen INT     min insertion length
  --bed               generate breakpoints bed file
";

  my $alignFile = undef;
  my $bedFormat = 0;
  my $minDel = 50;
  my $minIns = 50;
  my $verbose = 0;
	my $help = 0;
  ## my $debugEdge = "";
	
	GetOptions (
	"alignment=s"  => \$alignFile,
	"bed!"         => \$bedFormat,
	"mindellen=i"  => \$minDel,
	"mininslen=i"  => \$minIns,
	"verbose!"     => \$verbose,
	"help!"        => \$help)
	or die("Error in command line arguments\n$G_USAGE");
	
	die "$G_USAGE" if ($help);

  open INFILE, $alignFile || die "Fail to open $alignFile\n$!\n";
  my @Hs=();
  my $nH=0;
  my %BPs=();
  my %CGs=();
  while (<INFILE>) {
    if (0==$nH) {
      chomp();
      @Hs = split(/\t/);
      $nH=scalar(@Hs);
    } else {
      chomp();
      my @Cs=split(/\t/);
      my %item = ();
      for(my $i=0; $i<$nH; ++$i) {
        $item{$Hs[$i]} = $Cs[$i];
      }
      my ($c,$s,$e,$strand,$len) = $item{refStartEndLen} =~ /(\w+):(\d+)\-(\d+)\[([\+\-])(\d+)\]/;
      $BPs{$c} = {starts=>{},ends=>{}} if (!exists $BPs{$c});
      $BPs{$c}->{starts}->{$s}++;
      $BPs{$c}->{ends}->{$e}++;

      # cache first and last entry!
      if (!exists $CGs{$item{readId}}) {
        my %CG = (firstRef=>\%item, lastRef=>\%item);
        $CGs{$item{readId}} = \%CG;
      } else {
        $CGs{$item{readId}}->{lastRef} = \%item;
      }

      my @indels = ();
      getIndelsCoordinateFromCigar($item{cigar}, ($item{flag} & 0x10), 1, \@indels, $minDel, $minIns);
      foreach my $indelRef (@indels) {
        if ("d" eq $indelRef->{type}) {
          my $delend = $s+$indelRef->{start0} - 1;
          $BPs{$c}->{ends}->{$delend}++;
          my $delstart = $s+$indelRef->{end};
          $BPs{$c}->{starts}->{$delstart}++;
        }
      }
    }
  }
  close INFILE;

  # let's remove artifical break from linearizing a circular genome
  my $windowBp = 5;
  my %toRemoveStarts = ();
  my %toRemoveEnds = ();
  while (my ($id, $pairRef) = each %CGs) {
      my ($c1,$s1,$e1,$strand1,$len1) = $pairRef->{firstRef}->{refStartEndLen} =~ /(\w+):(\d+)\-(\d+)\[([\+\-])(\d+)\]/;
      my ($c2,$s2,$e2,$strand2,$len2) = $pairRef->{lastRef}->{refStartEndLen} =~ /(\w+):(\d+)\-(\d+)\[([\+\-])(\d+)\]/;
      if ($c1 eq $c2) {
        if ('+' eq $strand1) {
          if ('+' eq $strand2) {
            my $diff = $e2 - $s1;
            ## if ($s1 == ($e2+1)) {
            if (abs($diff)<=$windowBp) {
              $BPs{$c1}->{starts}->{$s1}--;
              $BPs{$c2}->{ends}->{$e2}--;

              if (0==$BPs{$c1}->{starts}->{$s1}) {
                if (!exists $toRemoveStarts{$c1}) {
                  my %empty = ();
                  $toRemoveStarts{$c1} = \%empty;
                }
                $toRemoveStarts{$c1}->{$s1}++;
              }
              if (0==$BPs{$c2}->{ends}->{$e2}) {
                if (!exists $toRemoveEnds{$c2}) {
                  my %empty = ();
                  $toRemoveEnds{$c2} = \%empty;
                }
                $toRemoveEnds{$c2}->{$e2}++;
              }
            }
          } else {
            # FIXME: possible?
          }
        } else {
          if ('+' eq $strand2) {
            # FIXME: possible?
          } else {
            my $diff = $s2 - $e1;
            ## if (($e1+1) == $s2+1) {
            if (abs($diff)<=$windowBp) {
              $BPs{$c1}->{ends}->{$e1}--;
              $BPs{$c2}->{starts}->{$s2}--;

              if (0==$BPs{$c1}->{ends}->{$e1}) {
                if (!exists $toRemoveEnds{$c1}) {
                  my %empty = ();
                  $toRemoveEnds{$c1} = \%empty;
                }
                $toRemoveEnds{$c1}->{$e1}++;
              }
              if (0==$BPs{$c2}->{starts}->{$s2}) {
                if (!exists $toRemoveStarts{$c2}) {
                  my %empty = ();
                  $toRemoveStarts{$c2} = \%empty;
                }
                $toRemoveStarts{$c2}->{$s2}++;
              }
            }
          }
        }
      }
  }
  # let's remove any genomic position that has no more support
  foreach my $chrom (keys %toRemoveStarts) {
    foreach my $pos (keys %{$toRemoveStarts{$chrom}}) {
      delete $BPs{$chrom}->{starts}->{$pos};
    }
  }
  foreach my $chrom (keys %toRemoveEnds) {
    foreach my $pos (keys %{$toRemoveEnds{$chrom}}) {
      delete $BPs{$chrom}->{ends}->{$pos};
    }
  }
  # FIXME: remove empty $chrom

  print "coordType\tchrom\tposition\n";
  foreach my $chrom (sort keys %BPs) {
    foreach my $pos (sort {$a<=>$b} keys %{$BPs{$chrom}->{starts}}) {
      printf "start\t%s\t%d\n", $chrom, $pos;
    }
  }
  foreach my $chrom (sort keys %BPs) {
    foreach my $pos (sort {$a<=>$b} keys %{$BPs{$chrom}->{ends}}) {
      printf "end\t%s\t%d\n", $chrom, $pos;
    }
  }

  if (0!=$bedFormat) {
    my $outfile = sprintf("%s.bed", $alignFile);
    open OUTFILE, ">$outfile" || die "Fail to open $outfile\n$!\n";
    foreach my $chrom (sort keys %BPs) {
      my @rows = ();
      foreach my $pos (keys %{$BPs{$chrom}->{starts}}) {
        push @rows, {pos=>$pos, type=>'b'};
      }
      foreach my $pos (keys %{$BPs{$chrom}->{ends}}) {
        push @rows, {pos=>$pos, type=>'e'};
      }
      @rows = sort { $a->{pos}<=>$b->{pos} || $a->{type} cmp $b->{type} } @rows;
      foreach my $rowRef (@rows) {
        printf OUTFILE "%s\t%d\t%d\t%s\n", $chrom, $rowRef->{pos}-1, $rowRef->{pos}, 
          sprintf("%s:%d_%s", $chrom, $rowRef->{pos}, ("b" eq $rowRef->{type})?"start":"end");
      }
    }
    close OUTFILE;
  }
}

sub getSortableCGIdOld {
  my ($id) = @_;
  my @bits = split(/[\.\/]/, $id);
  my $numBits = scalar(@bits);
  if (2==$numBits) {
    if ($bits[1]=~/^sG\d+\_\d+\_\d+[A-Za-z]*$/) {
      my @parts = split(/\_/, $bits[1]);
      if ($parts[$#parts]=~/\d+[A-Za-z]+$/) {
        my ($num, $words) = $parts[$#parts]=~/(\d+)([A-Za-z]+)$/;
        $parts[$#parts] = sprintf("%02d%s", int($num), $words);
      } else {
        $parts[$#parts] = sprintf("%02d", int($parts[$#parts]));
      }

      $bits[1] = join("_", @parts);
      return join(".", @bits);
    }
  } elsif (3==$numBits) {
    if ($bits[2] eq 'rc') {
      if ($bits[1]=~/^sG\d+\_\d+\_\d+[A-Za-z]*$/) {
        my @parts = split(/\_/, $bits[1]);
        if ($parts[$#parts]=~/\d+[A-Za-z]+$/) {
          my ($num, $words) = $parts[$#parts]=~/(\d+)([A-Za-z]+)$/;
          $parts[$#parts] = sprintf("%02d%s", int($num), $words);
        } else {
          $parts[$#parts] = sprintf("%02d", int($parts[$#parts]));
        }
        $bits[1] = join("_", @parts);
        return sprintf("%s.%s/rc", $bits[0], $bits[1]);
      }
    }
  }
  return $id;
}

sub getSortableCGId {
  my ($id) = @_;

  my @parts = split(/\_/, $id);
  if (scalar(@parts)>1) {
    if ($parts[$#parts]=~/rot\d+[A-Za-z]*/) {
      if ($parts[$#parts - 1]=~/\d+[\/]{0,1}[A-Za-z]+$/) {
        my ($num, $words) = $parts[$#parts - 1]=~/(\d+)([\/]{0,1}[A-Za-z]+)$/;
        $parts[$#parts - 1] = sprintf("%02d%s", int($num), $words);
      } else {
        $parts[$#parts - 1] = sprintf("%02d", int($parts[$#parts - 1]));
      }
    } elsif ($parts[$#parts]=~/\d+[\/]{0,1}[A-Za-z]+$/) {
      my ($num, $words) = $parts[$#parts]=~/(\d+)([\/]{0,1}[A-Za-z]+)$/;
      $parts[$#parts] = sprintf("%02d%s", int($num), $words);
    } else {
      $parts[$#parts] = sprintf("%02d", int($parts[$#parts]));
    }
    return join("_", @parts);
  }
  return $id;
}

sub getReadLenFromCigar {
  my ($cigar) = @_;
  my $readLen = 0;
  my @ops = split(/([MIDNSHP=X])/, $cigar);
  for(my $i=0; $i<scalar(@ops); ++$i) {
    my $count = int($ops[$i]);
    $i++;
    my $op = $ops[$i];
    if ('M' eq $op || '=' eq $op || 'X' eq $op) {
      $readLen += $count;
    } elsif ('I' eq $op || 'S' eq $op) {
      $readLen += $count;
    } elsif ('D' eq $op || 'N' eq $op) {
    } elsif ('H' eq $op) {
      $readLen += $count;
    }
  }
  return $readLen;
}

sub tagBam {
  my $G_USAGE = "
$0 tablebam --bam <bamFile> --alignment <alignmentFile> 
  --bam STR           bam file
  --alignment STR     alignments file
  --exclude STR       query to exclude (multiple)
";

  my $bamFile = undef;
  my $alignFile = undef;
  my @excludeIds = ();
  my $verbose = 0;
	my $help = 0;
  ## my $debugEdge = "";
	
	GetOptions (
	"bam=s"        => \$bamFile,
	"alignment=s"  => \$alignFile,
	"exclude=s"    => \@excludeIds,
	"verbose!"     => \$verbose,
	"help!"        => \$help)
	or die("Error in command line arguments\n$G_USAGE");
	
	die "$G_USAGE" if ($help);

  my %excludeIds = ();
  # grep { $excludeIds{$_} = 1; } @excludeIds;
  grep { $excludeIds{getSortableCGId($_)} = 1; } @excludeIds;

  my %wantedIds = ();
  open INFILE, $alignFile || die "Fail to open $alignFile\n$!\n";
  my @Hs=();
  my $nH=0;
  my %BPs=();
  while (<INFILE>) {
    if (0==$nH) {
      chomp();
      @Hs = split(/\t/);
      $nH=scalar(@Hs);
    } else {
      chomp();
      my @Cs=split(/\t/);
      my %item = ();
      for(my $i=0; $i<$nH; ++$i) {
        $item{$Hs[$i]} = $Cs[$i];
      }
      next if (exists $excludeIds{$item{readId}});
      $wantedIds{$item{readId}}++;
    }
  }
  close INFILE;

  # for recording only
  my @wantedIds = sort keys %wantedIds;
  my $outfile = sprintf("%s.rep.readids", $alignFile);
  open OUTFILE, ">$outfile" || die "Fail to open $outfile\n$!\n";
  print OUTFILE join("\n", @wantedIds), "\n";
  close OUTFILE;
  print "Query IDs written to $outfile\n";

  my $bamPrefix = $bamFile;
  $bamPrefix =~ s/\.bam$//;
  # generate bam by query
  my $queryBam = sprintf("%s.rep.anno.query.bam", $bamPrefix);

  open INFILE, "samtools view -h -F 0x700 $bamFile | " || die "Fail to open $bamFile\n$!\n";
  open OUTFILE, "| samtools view -b -o $queryBam -" || die "Fail to open $queryBam\n$!\n";
  my @readids = ();
  my %reads = ();
  my %refLens = ();
  my $numRefs = 0;
  while (<INFILE>) {
    if (/^@/) {
      if (/^\@SQ\t/) {
        my @bits = split(/\t/);
        my $sqsn = '';
        my $sqln = 0;
        foreach my $bit (@bits) {
          my ($tag, $value) = split(/:/, $bit, 2);
          if ('SN' eq $tag) {
            $sqsn = $value;
          } elsif ('LN' eq $tag) {
            $sqln = int($value);
          }
        }
        if ('' ne $sqsn) {
          $refLens{$sqsn} = $sqln;
          $numRefs++;
        }
      }
      print OUTFILE $_;
      next;
    }
    chomp();
    my @bits = split(/\t/);
    next if (!exists $wantedIds{$bits[0]});
    
    if (!exists $reads{$bits[0]}) {
      push @readids, $bits[0];
      $reads{$bits[0]} = {id=>$bits[0],records=>[]};
    }
    my %record = (read=>$bits[0],flag=>$bits[1],contig=>$bits[2],start=>$bits[3],mapq=>$bits[4],cigar=>$bits[5]);
    $record{samrecord} = \@bits;
    push @{$reads{$bits[0]}->{records}}, \%record;

    # calculate the coordinates
    my ($read_start0, $read_end1, $ref_start0, $ref_end1) = cigarToCoordinates($record{cigar}, ($record{flag} & 0x10), 1);
    $record{read_start0} = $read_start0;
    $record{read_end1} = $read_end1;
    $record{ref_start0} = $ref_start0;
    $record{ref_end1} = $ref_end1;
  }
  close INFILE;
  print "$bamFile read.\n";


  @readids = sort @readids;
  foreach my $readId (@readids) {
    my $readRef = $reads{$readId};
    my @records = sort { $a->{read_start0}<=>$b->{read_start0} || $a->{read_end1}<=>$b->{read_end1} } @{$readRef->{records}};
    my $numAlignments = scalar(@records);
    my $alignId = 0;
    my $readLen = 0;
    my $refLen = 0;
    my $numRecords = scalar(@records);
    my $prevRecRef = undef;
    my $nextRecRef = undef;
    for(my $i=0; $i<$numRecords; ++$i) {
      my $recordRef = $records[$i];
      # TODO: remove unnecessary computations
      # if (0!=$easyCoord && 0==$i) {
      #   last if ($recordRef->{flag} & 0x10);
      # }
      $prevRecRef = $records[$i-1] if ($i>0);
      $nextRecRef = $records[$i+1] if ($i<($numRecords+1));
      $alignId++;
      $refLen = (exists $refLens{$recordRef->{contig}}) ? $refLens{$recordRef->{contig}} : 0;
      if (0==$readLen) {
        $readLen = getReadLenFromCigar($recordRef->{cigar});
      }
      my $overlapStatus = 'n.a.';
      if (defined $prevRecRef) {
        if ($prevRecRef->{read_end1}>=($recordRef->{read_start0}+1)) {
          my $obp = $prevRecRef->{read_end1} - $recordRef->{read_start0};
          $overlapStatus = 'start('.$obp.')';
        }
      }
      if (defined $nextRecRef) {
        if ($recordRef->{read_end1}>=($nextRecRef->{read_start0}+1)) {
          my $obp = $recordRef->{read_end1} - $nextRecRef->{read_start0};
          if ($overlapStatus ne 'n.a.') {
            $overlapStatus .= ',end('.$obp.')';
          } else {
            $overlapStatus = 'end('.$obp.')';
          }
        }
      }

      # inject the tags ZG, ZO replace if necessary
      my $sortableCGId = getSortableCGId($recordRef->{read});
      my $ZG = sprintf("%s", $sortableCGId);
      my $ZGfound = 0;
      my $ZO = sprintf("%s_%02d/%02d", $sortableCGId, $alignId, $numAlignments);
      my $ZOfound = 0;
      my $bitsRef = $recordRef->{samrecord};
      for(my $j=11; $j<scalar(@{$bitsRef}); ++$j) {
        if ("ZG" eq substr($bitsRef->[$j],0,2)) {
          $ZGfound = 1;
          $bitsRef->[$j] = "ZG:Z:" . $ZG;
        } elsif ("ZO" eq substr($bitsRef->[$j],0,2)) {
          $ZOfound = 1;
          $bitsRef->[$j] = "ZO:Z:" . $ZO;
        }
      }
      if (0==$ZGfound) {
        push @{$bitsRef}, "ZG:Z:" . $ZG;
      }
      if (0==$ZOfound) {
        push @{$bitsRef}, "ZO:Z:" . $ZO;
      }
      print OUTFILE join("\t", @{$bitsRef}), "\n";
    }
  }
  close OUTFILE;
  print "$queryBam written.\n";

  # TODO: generate bam by coordinates
  # samtools sort -o A12.sup.sG5_1.CGs.rep.anno.bam -
  my $coordBam = sprintf("%s.rep.anno.bam", $bamPrefix);
  my @commands = ('samtools', 'sort', '-o', $coordBam, $queryBam);
  system(@commands) == 0
      or die "system @commands failed: $?";
  print "$coordBam written.\n";

  # TODO: index coord bam
  # samtools index A12.sup.sG5_1.CGs.rep.anno.bam
  @commands = ('samtools', 'index', $coordBam);
  system(@commands) == 0
      or die "system @commands failed: $?";
  print "$coordBam.bai written.\n";
}