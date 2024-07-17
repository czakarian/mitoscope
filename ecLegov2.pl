#!/usr/bin/env perl

##
## ecLego2 v2.04.01
## Copyright (C) 2019-2024 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##

## TODO: samtools command!!!
## TODO: fix internal "edge_" with data type accordingly! (external fixed)

use strict;
use Getopt::Long;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;
use File::Basename;

my $G_USAGE = "
$0 <command> -h
<command> [sievegraph]
sievegraph  : sieve the assembly graphs for cycles
graphbam    : augment bam file with graph information
initsubgraph: initialize subgraph for checks
spanedge    : span a graph from a given edge

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

# global variables to be localized
# subcommand = sievegraph
my @covThresholds = (
  {'lt'=>3,'width'=>1, label=>'0_chrom'}
  , {'lt'=>5,'width'=>2, label=>'1_low_amp'}
  , {'lt'=>9,'width'=>3, label=>'2_med_amp'}
  , {'width'=>4, label=>'3_high_amp'}
);
my @G_ReportBins = ('0_chrom', '1_low_amp', '2_med_amp', '3_high_amp');
my %G_ReportBins = ();
grep { $G_ReportBins{$_}=1; } @G_ReportBins;

my $debugEdge = "";
# END - global variables to be localized

if ('sievegraph' eq $command) {
	sieveGraph();
} elsif ('graphbam' eq $command) {
	graphBam();
} elsif ('initsubgraph' eq $command) {
	initSubgraph();
} elsif ('spanedge' eq $command) {
	spanEdge();
} else {
	print $G_USAGE;
}

exit 0;

# okie
sub getCoverageBin {
  my ($cov, $cov2cn) = @_;
  my $covcn = 2*$cov/$cov2cn;
  foreach my $thresholdRef (@covThresholds) {
    if (exists $thresholdRef->{lt}) {
      if ($covcn < $thresholdRef->{lt}) {
        return $thresholdRef->{label};
      }
    } else {
      return $thresholdRef->{label};
    }
  }
}

# okie
sub getCoveragePenWidth {
  my ($cov, $cov2cn) = @_;
  my $covcn = 2*$cov/$cov2cn;
  foreach my $thresholdRef (@covThresholds) {
    if (exists $thresholdRef->{lt}) {
      if ($covcn < $thresholdRef->{lt}) {
        return $thresholdRef->{width};
      }
    } else {
      return $thresholdRef->{width};
    }
  }
}

# okie
sub writeReadIds {
  my ($graphRef, $graphEdgeType, $file) = @_;
  open OUTFILE, ">$file" || die "Fail to open $file\n$!\n";
  if (exists $graphRef->{rootnodes}) {
    foreach my $rootNodeRef (@{$graphRef->{rootnodes}}) {
      foreach my $edgeRef (@{$rootNodeRef->{edges}}) {
        printf OUTFILE "%s_%d\n", $graphEdgeType, $edgeRef->{id};
      }
    }
  } elsif (exists $graphRef->{edges}) {
    foreach my $edgeRef (@{$graphRef->{edges}}) {
      print OUTFILE "%s_%d\n", $graphEdgeType, $edgeRef->{id};
    }
  }

  close OUTFILE;
}

# okie
sub isCyclic {
  my ($id, $nodesRef, $nodesAdjRef, $visitedRef, $nodesInStackRef, $endid) = @_;
  if (!exists $nodesInStackRef->{ordered}) {
    $nodesInStackRef->{ordered} = [];
    $nodesInStackRef->{named} = {};
  }
  if (!exists $visitedRef->{$id} || 0==$visitedRef->{$id}) {
    $visitedRef->{$id} = 1;
    push @{$nodesInStackRef->{ordered}}, $id;
    $nodesInStackRef->{named}->{$id} = scalar(@{$nodesInStackRef->{ordered}});
    if (defined $nodesAdjRef->[$id]) {
      for(my $i=0; $i<scalar(@{$nodesAdjRef->[$id]}); ++$i) {
        my $nextChildRef = $nodesAdjRef->[$id]->[$i];
        if ((!exists $visitedRef->{$nextChildRef->{dest}} || 0==$visitedRef->{$nextChildRef->{dest}}) 
          && 0!=isCyclic($nextChildRef->{dest}, $nodesRef, $nodesAdjRef, $visitedRef, $nodesInStackRef, $endid)) {
          $nodesInStackRef->{hits} = [] if (!exists $nodesInStackRef->{hits});
          push @{$nodesInStackRef->{hits}}, $id;
          return 1;
        } elsif (exists $nodesInStackRef->{named}->{$nextChildRef->{dest}}) {
          # back edge
          if (!defined $endid) {
            $nodesInStackRef->{hits} = [] if (!exists $nodesInStackRef->{hits});
            push @{$nodesInStackRef->{hits}}, $id;
            return 1;
          } elsif ($nextChildRef->{dest} == $endid) {
            $nodesInStackRef->{hits} = [] if (!exists $nodesInStackRef->{hits});
            push @{$nodesInStackRef->{hits}}, $id;
            return 1;
          }
        }
      }
    }
  }
  delete $nodesInStackRef->{named}->{$id};
  pop @{$nodesInStackRef->{ordered}};
  return 0;
}

# okie
sub getUniqueSubgraphs {
  my ($rootNodesRef) = @_;
  my $count = 0;
  foreach my $rootNodeRef (@{$rootNodesRef}) {
    if (exists $rootNodeRef->{root}->{$rootNodeRef->{id}}) {
      if (!exists $rootNodeRef->{duplicated} || 0==$rootNodeRef->{duplicated}) {
        $count++;
      }
    }
  }
  return $count;
}

# okie
sub getDisjointGraphRowMetrics {
  my ($graphRef, $cov2cn, $metricRef) = @_;

  # split the data into two set
  my @cmembers = ();
  my @cedges = ();
  my @amembers = ();
  my @aedges = ();
  die "Graph must have attributes cyclicGraphs and acyclicGraphs!\n", Dumper($graphRef), "\n" if (!exists $graphRef->{cyclicGraphs});
  foreach my $rootNodeRef (@{$graphRef->{cyclicGraphs}}) {
    push @cmembers, @{$rootNodeRef->{members}};
    push @cedges, @{$rootNodeRef->{edges}};
  }
  foreach my $rootNodeRef (@{$graphRef->{acyclicGraphs}}) {
    push @amembers, @{$rootNodeRef->{members}};
    push @aedges, @{$rootNodeRef->{edges}};
  }

  if (1==1) {
    my %nodeIds = ();
    grep { $nodeIds{$_} = 1; } @cmembers;
    my @vertices = sort { int($a)<=>int($b) } keys %nodeIds;
    $metricRef->{vertices} = \@vertices;
    $metricRef->{vlen} = scalar(@vertices);

    my %edgeIds = ();
    grep { $edgeIds{$_->{id}} = 1; } @cedges;
    my @edges = sort { int($a)<=>int($b) } keys %edgeIds;
    $metricRef->{edges} = \@edges;
    $metricRef->{elen} = scalar(@edges);
  }

  my %covBins = ();
  grep { $covBins{$_} = {count=>0, covs=>[], lens=>[], lenKb=>0, cov=>0, edgeIds=>[]}; } @G_ReportBins;
  foreach my $edgeRef (@cedges) {
    my $covBin = getCoverageBin($edgeRef->{cov}, $cov2cn);
    my $cn = 2 * int($edgeRef->{cov}) / $cov2cn;
    my $cnLabel = sprintf '%.1fc', $cn;
    $covBins{$covBin}->{count}++;
    if (exists $G_ReportBins{$covBin} && 0!=$G_ReportBins{$covBin}) {
      push @{$covBins{$covBin}->{covs}}, $cnLabel;
      push @{$covBins{$covBin}->{lens}}, $edgeRef->{len};
      push @{$covBins{$covBin}->{edgeIds}}, $edgeRef->{id};
      my $tlen = $edgeRef->{len};
      $tlen =~ s/k$//i;
      $tlen *= 1.0;
      $covBins{$covBin}->{lenKb} += $tlen;
      $covBins{$covBin}->{cov} += ($cn * $tlen);
    }
  }
  grep { $covBins{$_}->{cov} = $covBins{$_}->{cov} / $covBins{$_}->{lenKb} if ($covBins{$_}->{count}>0); } @G_ReportBins;
  $metricRef->{covBins} = \%covBins;

  if (1==1) {
    my %nodeIds = ();
    grep { $nodeIds{$_} = 1; } @amembers;
    my @vertices = sort { int($a)<=>int($b) } keys %nodeIds;
    $metricRef->{avertices} = \@vertices;
    $metricRef->{avlen} = scalar(@vertices);

    my %edgeIds = ();
    grep { $edgeIds{$_->{id}} = 1; } @aedges;
    my @edges = sort { int($a)<=>int($b) } keys %edgeIds;
    $metricRef->{aedges} = \@edges;
    $metricRef->{aelen} = scalar(@edges);
  }
}

# okie
sub writeGraph {
  my ($graphRef, $graphEdgeType, $commentsRef, $cov2cn, $minlen, $file) = @_;
  # write out the graph
  open OUTFILE, ">$file" || die "Fail to open $file\n$!\n";
  print OUTFILE "digraph {\n";
  foreach my $comment (@{$commentsRef}) {
    print OUTFILE "// $comment\n";
  }
  print OUTFILE "nodesep = 1.0;\n";
  print OUTFILE "node [shape = point, label = \"\", height = 0.2, fillcolor=\"lightgrey\", color=\"grey\"];\n";
  print OUTFILE "edge [fontsize=\"20pt\", fontname=\"Arial\"];\n";
  # "0" -> "1" [label = "id 4\l22k 724x", color = "lightgrey", penwidth=3] ;
  foreach my $edgeRef (@{$graphRef->{edges}}) {
    # (src=>$src, dest=>$dest, id=>$id, strand=>$strand, len=>$len, cov=>$cov);
    printf OUTFILE '"%d" -> "%d" [', $edgeRef->{src}, $edgeRef->{dest};
    my $edgeId = $edgeRef->{id};
    $edgeId = '-'.$edgeId if ('-' eq $edgeRef->{strand});
    my $cn = 2 * int($edgeRef->{cov}) / $cov2cn;
    # printf OUTFILE 'label="id %s\\l%s %s"', $edgeId, $edgeRef->{len}, $edgeRef->{cov};
    printf OUTFILE 'label="id %s\\l%s %.1fc"', $edgeId, $edgeRef->{len}, $cn;
    # printf ', color="%s"', $edgeRef->{color};
    my $penWidth = getCoveragePenWidth($edgeRef->{cov}, $cov2cn);
    printf OUTFILE ',penwidth=%d', $penWidth;
    if (exists $edgeRef->{class} && 'acyclic' eq $edgeRef->{class}) {
      printf OUTFILE ',color="grey",fontcolor="grey"', $penWidth;
    }
    print OUTFILE "];\n";
  }
  print OUTFILE "}\n";
  close OUTFILE;

  open DETAILFILE, ">$file.details" || die "Fail to open $file.details\n$!\n";
  foreach my $comment (@{$commentsRef}) {
    print DETAILFILE "// $comment\n";
  }

  # row-based for ranking
  my @dataRows = ();

  my $subgraphCount = 0;
  foreach my $rootNodeRef (@{$graphRef->{rootnodes}}) {
    if (!exists $rootNodeRef->{duplicated} || 0==$rootNodeRef->{duplicated}) {
      $subgraphCount++;
      my %dataRow = (id=>$subgraphCount);
      getDisjointGraphRowMetrics($rootNodeRef, $cov2cn, \%dataRow);

      # FIXME: propagate id for future use
      $rootNodeRef->{subgraphId} = $subgraphCount;

      printf DETAILFILE "sG#%d: |V|=%d, |E|=%d, |linearV|=%d, |linearE|=%d\n", $dataRow{id}, $dataRow{vlen}, $dataRow{elen}, $dataRow{avlen}, $dataRow{aelen};
      printf DETAILFILE "sG#%d V: %s\n", $dataRow{id}, join(",", @{$dataRow{vertices}});
      printf DETAILFILE "sG#%d E: %s\n", $dataRow{id}, join(",", @{$dataRow{edges}});

      push @dataRows, \%dataRow;

      my @covBinClasses = ();
      foreach my $covBin (sort { $a cmp $b } keys %{$dataRow{covBins}}) {
        my $label = sprintf("|%s|=%d", $covBin, $dataRow{covBins}->{$covBin}->{count});
        push @covBinClasses, $label;
      }
      printf DETAILFILE "sG#%d class(E): %s\n", $dataRow{id}, join(", ", @covBinClasses);

      my @covBinCounts = ();
      my @covBinLens = ();
      foreach my $covBin (sort { $a cmp $b } keys %{$dataRow{covBins}}) {
        if (scalar(@{$dataRow{covBins}->{$covBin}->{covs}})>0) {
          printf DETAILFILE "sG#%d cov(E, %s): %d=|{%s}|\n", $dataRow{id}, $covBin, $dataRow{covBins}->{$covBin}->{count}, join(", ", @{$dataRow{covBins}->{$covBin}->{covs}});
          printf DETAILFILE "sG#%d len(E, %s): %d=|{%s}|\n", $dataRow{id}, $covBin, $dataRow{covBins}->{$covBin}->{count}, join(", ", @{$dataRow{covBins}->{$covBin}->{lens}});
        }
      }

      open READFILE, ">$file.subgraph_$subgraphCount.readids" || die "Fail to open $file.subgraph_$subgraphCount.readids\n$!\n";
      foreach my $edgeId (@{$dataRow{edges}}) {
        printf READFILE "%s_%d\n", $graphEdgeType, $edgeId;
      }
      close READFILE;

      open SUBGRAPHFILE, ">$file.subgraph_$subgraphCount.gv" || die "Fail to open $file.subgraph_$subgraphCount.gv\n$!\n";
      print SUBGRAPHFILE "digraph {\n";
      foreach my $comment (@{$commentsRef}) {
        print SUBGRAPHFILE "// $comment\n";
      }
      printf SUBGRAPHFILE "// sG#%d: |V|=%d, |E|=%d, |linearV|=%d, |linearE|=%d\n", $dataRow{id}, $dataRow{vlen}, $dataRow{elen}, $dataRow{avlen}, $dataRow{aelen};
      printf SUBGRAPHFILE "// sG#%d V: %s\n", $dataRow{id}, join(",", @{$dataRow{vertices}});
      printf SUBGRAPHFILE "// sG#%d E:%s\n", $dataRow{id}, join(",", @{$dataRow{edges}});
      printf SUBGRAPHFILE "// sG#%d class(E): %s\n", $dataRow{id}, join(", ", @covBinClasses);
      foreach my $covBin (sort { $a cmp $b } keys %{$dataRow{covBins}}) {
        if (scalar(@{$dataRow{covBins}->{$covBin}->{covs}})>0) {
          printf SUBGRAPHFILE "// sG#%d cov(E, %s): %d=|{%s}|\n", $dataRow{id}, $covBin, $dataRow{covBins}->{$covBin}->{count}, join(", ", @{$dataRow{covBins}->{$covBin}->{covs}});
          printf SUBGRAPHFILE "// sG#%d len(E, %s): %d=|{%s}|\n", $dataRow{id}, $covBin, $dataRow{covBins}->{$covBin}->{count}, join(", ", @{$dataRow{covBins}->{$covBin}->{lens}});
        }
      }
      print SUBGRAPHFILE "nodesep=1.0;\n";
      print SUBGRAPHFILE "node [shape=point,label=\"\",height = 0.2,fillcolor=\"lightgrey\",color=\"grey\"];\n";
      print SUBGRAPHFILE "edge [fontsize=\"20pt\",fontname=\"Arial\"];\n";
      foreach my $edgeRef (@{$rootNodeRef->{edges}}) {
        # (src=>$src, dest=>$dest, id=>$id, strand=>$strand, len=>$len, cov=>$cov);
        printf SUBGRAPHFILE '"%d" -> "%d" [', $edgeRef->{src}, $edgeRef->{dest};
        my $edgeId = $edgeRef->{id};
        $edgeId = '-'.$edgeId if ('-' eq $edgeRef->{strand});
        my $cn = 2 * int($edgeRef->{cov}) / $cov2cn;
        printf SUBGRAPHFILE 'label="id %s\\l%s %.1fc"', $edgeId, $edgeRef->{len}, $cn;
        my $penWidth = getCoveragePenWidth($edgeRef->{cov}, $cov2cn);
        printf SUBGRAPHFILE ',penwidth=%d', $penWidth;
        if (exists $edgeRef->{class} && 'acyclic' eq $edgeRef->{class}) {
          printf SUBGRAPHFILE ',color="grey",fontcolor="grey"', $penWidth;
        }
        print SUBGRAPHFILE "];\n";
      }
      print SUBGRAPHFILE "}\n";

      close SUBGRAPHFILE;
    }
  }
  close DETAILFILE;

  # report tabular overview
  open OVERVIEWFILE, ">$file.overview.xls" || die "Fail to open $file.overview.xls\n$!\n";
  if ($minlen>0) {
    open OVERVIEWPASSFILE, ">$file.min$minlen\kb.overview.xls" || die "Fail to open $file.min$minlen\kb.overview.xls\n$!\n";
  }
  foreach my $comment (@{$commentsRef}) {
    print OVERVIEWFILE "// $comment\n";
    if ($minlen>0) {
      print OVERVIEWPASSFILE "// $comment\n";
    }
  }
  # @dataRows = sort { $a->{id} <=> $b->{id} } @dataRows;
  @dataRows = sort { 
    $b->{covBins}->{'3_high_amp'}->{cov} <=> $a->{covBins}->{'3_high_amp'}->{cov}
    || $b->{covBins}->{'3_high_amp'}->{lenKb} <=> $a->{covBins}->{'3_high_amp'}->{lenKb}
    || $b->{covBins}->{'3_high_amp'}->{count} <=> $a->{covBins}->{'3_high_amp'}->{count}
    || $b->{covBins}->{'2_med_amp'}->{cov} <=> $a->{covBins}->{'2_med_amp'}->{cov}
    || $b->{covBins}->{'2_med_amp'}->{lenKb} <=> $a->{covBins}->{'2_med_amp'}->{lenKb}
    || $b->{covBins}->{'2_med_amp'}->{count} <=> $a->{covBins}->{'2_med_amp'}->{count}
    || $b->{covBins}->{'1_low_amp'}->{cov} <=> $a->{covBins}->{'1_low_amp'}->{cov}
    || $b->{covBins}->{'1_low_amp'}->{lenKb} <=> $a->{covBins}->{'1_low_amp'}->{lenKb}
    || $b->{covBins}->{'1_low_amp'}->{count} <=> $a->{covBins}->{'1_low_amp'}->{count}
    || $b->{covBins}->{'0_chrom'}->{cov} <=> $a->{covBins}->{'0_chrom'}->{cov}
    || $b->{covBins}->{'0_chrom'}->{lenKb} <=> $a->{covBins}->{'0_chrom'}->{lenKb}
    || $b->{covBins}->{'0_chrom'}->{count} <=> $a->{covBins}->{'0_chrom'}->{count}
    || $b->{elen} <=> $a->{elen}
    || $b->{vlen} <=> $a->{vlen}
    || $a->{id} <=> $b->{id} 
    } @dataRows;
  my @cols = ('subgraphId', 'numLinearVertices', 'numLinearEdges', 'numVertices', 'numEdges');
  grep { push @cols, 'num_'.$_; push @cols, 'kbp_'.$_; push @cols, 'avecov_'.$_; } @G_ReportBins;
  push @cols, 'vertices', 'edges';
  grep { push @cols, 'cov_'.$_; push @cols, 'len_'.$_; } @G_ReportBins;
  push @cols, 'linearVertices', 'linearEdges';
  print OVERVIEWFILE join("\t", @cols), "\n";
  foreach my $dataRowRef (@dataRows) {
    @cols = ($dataRowRef->{id}, $dataRowRef->{avlen}, $dataRowRef->{aelen}, $dataRowRef->{vlen}, $dataRowRef->{elen});

    my $passFilter=0;
    foreach my $covBin (@G_ReportBins) {
      if (exists $dataRowRef->{covBins}->{$covBin}) {
        my $covBinRef = $dataRowRef->{covBins}->{$covBin};
        push @cols, $covBinRef->{count}, $covBinRef->{lenKb}, sprintf("%.1f", $covBinRef->{cov});
        if ($minlen>0) {
          $passFilter=1 if ($covBinRef->{lenKb}>=$minlen);
        }
      } else {
        push @cols, 0, 0;
      }
    }

    push @cols, join(";", @{$dataRowRef->{vertices}});
    push @cols, join(";", @{$dataRowRef->{edges}});
    
    foreach my $covBin (@G_ReportBins) {
      if (exists $dataRowRef->{covBins}->{$covBin}) {
        my $covBinRef = $dataRowRef->{covBins}->{$covBin};
        my @elements = ();
        for(my $j=0; $j<$covBinRef->{count}; ++$j) {
          push @elements, sprintf("%s[%d]", $covBinRef->{covs}->[$j], $covBinRef->{edgeIds}->[$j]);
        }
        push @cols, (0==$covBinRef->{count}) ? '.' : join(";", @elements);
        @elements = ();
        for(my $j=0; $j<$covBinRef->{count}; ++$j) {
          push @elements, sprintf("%s[%d]", $covBinRef->{lens}->[$j], $covBinRef->{edgeIds}->[$j]);
        }
        push @cols, (0==$covBinRef->{count}) ? '.' : join(";", @elements);
      } else {
        push @cols, '.', '.';
      }
    }

    push @cols, join(";", @{$dataRowRef->{avertices}});
    push @cols, join(";", @{$dataRowRef->{aedges}});

    print OVERVIEWFILE join("\t", @cols), "\n";
    if ($minlen>0) {
      if (0!=$passFilter) {
        print OVERVIEWPASSFILE join("\t", @cols), "\n";
      }
    }
  }
  close OVERVIEWFILE;
  if ($minlen>0) {
    close OVERVIEWPASSFILE;
  }

  # report html graphvizonline
  open HTMLFILE, ">$file.graphvizonline.html" || die "Fail to open $file.graphvizonline.html\n$!\n";
  print HTMLFILE '<html>', "\n";
  print HTMLFILE '<head>', "\n";
  print HTMLFILE '<title>', $file, '</title>', "\n";
  print HTMLFILE '</head>', "\n";
  print HTMLFILE '<body>', "\n";
  print HTMLFILE '<h1>Metainfo</h1>', "\n";
  foreach my $comment (@{$commentsRef}) {
    print HTMLFILE "$comment", '<br>', "\n";
  }
  print HTMLFILE '<h2>Ranked subgraphs</h2>', "\n";
  print HTMLFILE '<ol>', "\n";
  my $G_VIZURL='https://dreampuf.github.io/GraphvizOnline/#';
  foreach my $dataRowRef (@dataRows) {
    # printf HTMLFILE "<li> subgraph#%d |V|=%d |E|=%d\n", $dataRowRef->{id}, $dataRowRef->{vlen}, $dataRowRef->{elen};
    my $uri = '';
    my $subgraphFile = sprintf("%s.subgraph_%d.gv", $file, $dataRowRef->{id});
    open SUBGRAPHFILE, $subgraphFile || die "Fail to open $subgraphFile\n$!\n";
    while (my $line = <SUBGRAPHFILE>) {
      $line =~ s/([^^A-Za-z0-9\-_.!~*'()])/ sprintf "%%%02x", ord $1 /eg;
      $uri .= $line;
    }
    close SUBGRAPHFILE;

    if ($dataRowRef->{avlen}>0) {
      printf HTMLFILE "<li> <a href=\"%s%s\" target=\"%s\">subgraph#%d</a> |V|=%d |E|=%d |linearV|=%d |linearE|=%d\n", 
        $G_VIZURL, $uri, $subgraphFile, 
        $dataRowRef->{id}, $dataRowRef->{vlen}, $dataRowRef->{elen},
        $dataRowRef->{avlen}, $dataRowRef->{aelen};
    } else {
      printf HTMLFILE "<li> <a href=\"%s%s\" target=\"%s\">subgraph#%d</a> |V|=%d |E|=%d\n", 
        $G_VIZURL, $uri, $subgraphFile, 
        $dataRowRef->{id}, $dataRowRef->{vlen}, $dataRowRef->{elen};
    }
    print HTMLFILE "<div>\n";
    printf HTMLFILE "V=%s<br>", join(";", @{$dataRowRef->{vertices}});
    printf HTMLFILE "E=%s<br>", join(";", @{$dataRowRef->{edges}});
    if ($dataRowRef->{avlen}>0) {
      printf HTMLFILE "linearV=%s<br>", join(";", @{$dataRowRef->{avertices}});
      printf HTMLFILE "linearE=%s<br>", join(";", @{$dataRowRef->{aedges}});
    }
   foreach my $covBin (@G_ReportBins) {
      if (exists $dataRowRef->{covBins}->{$covBin} && $dataRowRef->{covBins}->{$covBin}->{count}>0) {
        my $covBinRef = $dataRowRef->{covBins}->{$covBin};
        printf HTMLFILE "%s: count=%d totalEdgesLen=%.1f Kb averageCov=%.1f<br>", 
          $covBin, $covBinRef->{count}, $covBinRef->{lenKb}, $covBinRef->{cov};

        my @elements = ();
        for(my $j=0; $j<$covBinRef->{count}; ++$j) {
          push @elements, sprintf("%s[%d]", $covBinRef->{covs}->[$j], $covBinRef->{edgeIds}->[$j]);
        }
        push @cols, join(";", @elements);
        printf HTMLFILE "cov(%s)=%s<br>", $covBin, join(";", @elements);

        @elements = ();
        for(my $j=0; $j<$covBinRef->{count}; ++$j) {
          push @elements, sprintf("%s[%d]", $covBinRef->{lens}->[$j], $covBinRef->{edgeIds}->[$j]);
        }
        push @cols, join(";", @elements);
        printf HTMLFILE "len(%s)=%s<br>", $covBin, join(";", @elements);
      }
    }

    # printf HTMLFILE "[<a href=\"%s%s\" target=\"%s\">view subgraph</a>]<br>", $G_VIZURL, $uri, $subgraphFile;
    print HTMLFILE "&nbsp;<br></div>\n";
  }
  print HTMLFILE '</ul>', "\n";
  print HTMLFILE '</body>', "\n";
  print HTMLFILE '</html>', "\n";
  close 
}

# okie
sub removeEdgesInsufficientCoverage {
  my ($graphRef, $coverage, $verbose) = @_;
  my @cnaEdges = ();
  foreach my $edgeRef (@{$graphRef->{edges}}) {
    # (src=>$src, dest=>$dest, id=>$id, strand=>$strand, len=>$len, cov=>$cov);
    if ($edgeRef->{cov}>=$coverage) {
      push @cnaEdges, $edgeRef;
      if (0!=$verbose) {
        printf "DEBUG: edge_%d (cov:%d) retained.\n", $edgeRef->{id}, $edgeRef->{cov};
      }
    } else {
      if (0!=$verbose) {
        printf "DEBUG: edge_%d (cov:%d<%d) dropped.\n", $edgeRef->{id}, $edgeRef->{cov}, $coverage;
      }
    }
  }
  $graphRef->{edges} = \@cnaEdges;
}

# okie
sub readGVGraphFile {
  my ($file, $graphRef) = @_;
  open INFILE, $file || die "Fail to open $file\n$!\n";
  my @edges = ();
  $graphRef->{edges} = \@edges;
  $graphRef->{edgeType} = 'edge';
  my $filename = basename($file);
  if ($filename =~ /assembly/i) {
    $graphRef->{edgeType} = 'contig';
  }
  
  $graphRef->{edges} = \@edges;
  while (<INFILE>) {
    next if (/^#/);
    chomp();
    # "0" -> "1" [label = "id 4\l22k 724x", color = "red" , penwidth = 3] ;
    # "2" -> "3" [label = "id -4\l22k 724x", color = "red" , penwidth = 3] ;
    my ($src, $dest, $remaining) = $_ =~ /\"(\d+)\"\s+\-\>\s+\"(\d+)\"\s+\[(.+)\]/;
    next if ((!defined $src) || (!defined $dest) || (!defined $remaining));
    $remaining =~ s/^\s+//;$remaining =~ s/\s+$//;
    my $id=''; my $strand=''; my $len=''; my $cov='';
    foreach my $element (split(/\s*,\s*/, $remaining)) {
      my ($key, $value) = split(/\s*\=\s*/, $element, 2);
      if (defined $key && uc($key) eq "LABEL") {
        $value =~ s/^\s*\"//;$value =~ s/\"\s*$//;
        my @bits = split(/\\l|\s+/, $value);
        if ('-' eq substr($bits[1],0,1)) {
          $id = -1 * $bits[1];
          $strand = '-';
        } else {
          $id = 1 * $bits[1];
          $strand = '+';
        }
        $len = $bits[2];
        $cov = $bits[3];
        last;
      }
    }
    my %edge = (src=>$src, dest=>$dest, id=>$id, strand=>$strand, len=>$len, cov=>$cov);
    push @edges, \%edge;
  }
}

#####
# refactoring
#####

sub subGraphDFS {
  my ($id, $nodesRef, $nodesAdjRef, $rootId, $visitedRef) = @_;
  if (!exists $visitedRef->{$id} || 0==$visitedRef->{$id}) {
    $visitedRef->{$id} = 1;
    if (!exists $nodesRef->{$id}) {
      # print STDERR "subGraphDFS: looking for $id under root $rootId, but $id does not exists!\n";
    } else {
      push @{$nodesRef->{$rootId}->{members}}, $nodesRef->{$id}->{id};
      $nodesRef->{$id}->{root}->{$rootId} = $nodesRef->{$rootId}->{id};
      if (defined $nodesAdjRef->[$id]) {
        my @ascEdgeRefs = sort { $a->{id}<=>$b->{id} || $a->{strand} cmp $b->{strand} } @{$nodesAdjRef->[$id]};
        for(my $i=0; $i<scalar(@ascEdgeRefs); ++$i) {
          my $nextChildRef = $ascEdgeRefs[$i];
          subGraphDFS($nextChildRef->{dest}, $nodesRef, $nodesAdjRef, $rootId, $visitedRef);
        }
      }
    }
  }
}

sub isRCGraphOriginal {
  my ($membersNodeRef, $nodesAdjRef, $visitedEdgesRef) = @_;
  foreach my $nodeId (@{$membersNodeRef}) {
    if (defined $nodesAdjRef->[$nodeId]) {
      foreach my $edgeRef (@{$nodesAdjRef->[$nodeId]}) {
        return 1 if (exists $visitedEdgesRef->{$edgeRef->{id}});
      }
    }
  }
  return 0;
}

sub isRCGraph {
  my ($membersNodeRef, $nodesAdjRef, $visitedEdgesRef) = @_;
  foreach my $nodeId (@{$membersNodeRef}) {
    if (defined $nodesAdjRef->[$nodeId]) {
      foreach my $edgeRef (@{$nodesAdjRef->[$nodeId]}) {
        return ($nodeId, $edgeRef->{id}) if (exists $visitedEdgesRef->{$edgeRef->{id}});
      }
    }
  }
  return (-1, -1);
}

sub recordVisitedEdges {
  my ($edgesRef, $visitedEdgesRef) = @_;
  grep { 
    $visitedEdgesRef->{$_->{id}} = 1;
  } @{$edgesRef};
}

sub regroupDisjointGraph {
  my ($graphRef, $subgraphId, $verbose) = @_;

   # build the graph for bidirectional traversal
  my $numEdges = scalar(@{$graphRef->{edges}});
  my %graphNodes = ();
  foreach my $edgeRef (@{$graphRef->{edges}}) {
    if (!exists $graphNodes{$edgeRef->{src}}) {
      $graphNodes{$edgeRef->{src}} = {id=>$edgeRef->{src}, self=>{}, next=>{}, prev=>{}}
    }
    if (!exists $graphNodes{$edgeRef->{dest}}) {
      $graphNodes{$edgeRef->{dest}} = {id=>$edgeRef->{dest}, self=>{}, next=>{}, prev=>{}}
    }
    if ($edgeRef->{src} == $edgeRef->{dest}) {
      my $graphNodeRef = $graphNodes{$edgeRef->{src}};
      $graphNodeRef->{self}->{$edgeRef->{dest}} = {id=>$edgeRef->{dest}, idRef=>$graphNodes{$edgeRef->{src}}, edgeRef=>[]} if (!exists $graphNodeRef->{self}->{$edgeRef->{dest}});
      push @{$graphNodeRef->{self}->{$edgeRef->{dest}}->{edgeRefs}}, $edgeRef;
    } else {
      my $graphNodeRef = $graphNodes{$edgeRef->{src}};
      $graphNodeRef->{next}->{$edgeRef->{dest}} = {id=>$edgeRef->{dest}, idRef=>$graphNodes{$edgeRef->{dest}}, edgeRefs=>[]} if (!exists $graphNodeRef->{next}->{$edgeRef->{dest}});
      push @{$graphNodeRef->{next}->{$edgeRef->{dest}}->{edgeRefs}}, $edgeRef;
      $graphNodeRef = $graphNodes{$edgeRef->{dest}};
      $graphNodeRef->{prev}->{$edgeRef->{src}} = {id=>$edgeRef->{src}, idRef=>$graphNodes{$edgeRef->{src}}, edgeRefs=>[]} if (!exists $graphNodeRef->{prev}->{$edgeRef->{src}});
      push @{$graphNodeRef->{prev}->{$edgeRef->{src}}->{edgeRefs}}, $edgeRef;
    }
  }
  my $numVertices = scalar(keys %graphNodes);
  # remove vertex without out-bound edge
  my %removedNodes = ();
  my @removedEdges = ();
  my %toCheckNodes = ();
  while (my ($nodeId, $nodeRef) = each %graphNodes) {
    $toCheckNodes{$nodeId} = $nodeRef;
  }
  while (scalar(keys %toCheckNodes)>0) {
    my @noOutBoundNodes = ();
    while (my ($nodeId, $nodeRef) = each %toCheckNodes) {
      if (0==scalar(keys %{$nodeRef->{next}})) {
        push @noOutBoundNodes, $nodeRef;
      } elsif (0==scalar(keys %{$nodeRef->{prev}})) {
        push @noOutBoundNodes, $nodeRef;
      }
    }
    %toCheckNodes = ();
    foreach my $nodeRef (@noOutBoundNodes) {
      $removedNodes{$nodeRef->{id}} = $nodeRef;
      if (0==scalar(keys %{$nodeRef->{next}})) {
        while (my ($prevNodeId, $prevNodeRef) = each %{$nodeRef->{prev}}) {
          push @removedEdges, @{$prevNodeRef->{edgeRefs}};
          delete $prevNodeRef->{idRef}->{next}->{$nodeRef->{id}};
          $toCheckNodes{$prevNodeId} = $prevNodeRef->{idRef};
        }
      }
      if (0==scalar(keys %{$nodeRef->{prev}})) {
        while (my ($nextNodeId, $nextNodeRef) = each %{$nodeRef->{next}}) {
          push @removedEdges, @{$nextNodeRef->{edgeRefs}};
          delete $nextNodeRef->{idRef}->{prev}->{$nodeRef->{id}};
          $toCheckNodes{$nextNodeId} = $nextNodeRef->{idRef};
        }
      }
    }
  }

  # let's make sure every nodes are invovled in a cyclic graph
  if (1==1) {
    my @edges = ();
    while (my ($nodeId, $nodeRef) = each %graphNodes) {
      next if (exists $removedNodes{$nodeId});
      while (my ($selfNodeId, $selfNodeRef) = each %{$nodeRef->{self}}) {
        push @edges, @{$selfNodeRef->{edgeRefs}};
      }
      while (my ($nextNodeId, $nextNodeRef) = each %{$nodeRef->{next}}) {
        push @edges, @{$nextNodeRef->{edgeRefs}};
      }
    }
    my @graphNodes = ();
    my %graphNodes = ();
    my @graphNodeAdj = ();
    foreach my $edgeRef (@edges) {
      if (!defined $graphNodeAdj[$edgeRef->{src}]) {
        $graphNodes{$edgeRef->{src}} = {id=>$edgeRef->{src}, members=>[], root=>{}};
        push @graphNodes, $graphNodes{$edgeRef->{src}};
        $graphNodeAdj[$edgeRef->{src}] = [];
      }
      if (!defined $graphNodeAdj[$edgeRef->{dest}]) {
        $graphNodes{$edgeRef->{dest}} = {id=>$edgeRef->{dest}, members=>[], root=>{}};
        push @graphNodes, $graphNodes{$edgeRef->{dest}};
        $graphNodeAdj[$edgeRef->{dest}] = [];
      }
      push @{$graphNodeAdj[$edgeRef->{src}]}, $edgeRef;
    }
    for(my $i=0; $i<scalar(@graphNodes); ++$i) {
      my %visited = ();
      subGraphDFS($graphNodes[$i]->{id}, \%graphNodes, \@graphNodeAdj, $graphNodes[$i]->{id}, \%visited);
    }

    # identify subgraph that contain cycle(s)
    my @cycles = ();
    my @noncycles = ();
    for(my $i=0; $i<scalar(@graphNodes); ++$i) {
      # check only from the root node!
      my %visited = ();
      my %nodesInStack = ();
      if (0!=isCyclic($graphNodes[$i]->{id}, \%graphNodes, \@graphNodeAdj, \%visited, \%nodesInStack, $graphNodes[$i]->{id})) {
        push @cycles, $graphNodes[$i];
        $graphNodes[$i]->{containCycle} = 1;
        my $minOffset = 0;
        my $minId = $nodesInStack{hits}->[$minOffset];
        for(my $j=1; $j<scalar(@{$nodesInStack{hits}}); ++$j) {
          if ($nodesInStack{hits}->[$j]<$minId) {
            $minId = $nodesInStack{hits}->[$j];
            $minOffset = $j;
          }
        }
        my @reorderedHits = ();
        for(my $j=$minOffset; $j<scalar(@{$nodesInStack{hits}}); ++$j) {
          push @reorderedHits, $nodesInStack{hits}->[$j];
        }
        for(my $j=0; $j<$minOffset; ++$j) {
          push @reorderedHits, $nodesInStack{hits}->[$j];
        }
        $graphNodes[$i]->{cyclicpath} = \@reorderedHits;
        if (0!=$verbose) {
          printf "DEBUG: %d is cylic in hits[%s] members[%s]\n", $graphNodes[$i]->{id}, join(",", @reorderedHits), join(",", @{$graphNodes[$i]->{members}});
        }
      } else {
        push @noncycles, $graphNodes[$i];
        $graphNodes[$i]->{containCycle} = 0;
        if (0!=$verbose) {
          printf "DEBUG: %d is acylic in members[%s]\n", $graphNodes[$i]->{id}, join(",", @{$graphNodes[$i]->{members}});
        }
      }
    }

    # summarize the results
    # 1) acyclic - all in a single graph or multiple?
    # 2) cyclic - single, nested or multiple?
    # a) remove cyclic edges into a group
    # b) remain will be acyclic edges
    my @cyclicEdges = ();
    my @acyclicEdges = ();
    if (1==1) {
      my %cyclicNodes = ();
      for(my $i=0; $i<scalar(@cycles); ++$i) {
        $cyclicNodes{$cycles[$i]->{id}}++;
      }
      foreach my $edgeRef (@{$graphRef->{edges}}) {
        if ($edgeRef->{src}==$edgeRef->{dest}) {
          push @cyclicEdges, $edgeRef;
        } elsif (exists $cyclicNodes{$edgeRef->{src}} && exists $cyclicNodes{$edgeRef->{dest}}) {
          push @cyclicEdges, $edgeRef;
        } else {
          push @acyclicEdges, $edgeRef;
        }
      }
    }

    if (1==1) {
      my @cNodes = ();
      my %cNodes = ();
      my @cNodeAdj = ();
      foreach my $edgeRef (@cyclicEdges) {
        if (!defined $cNodeAdj[$edgeRef->{src}]) {
          $cNodes{$edgeRef->{src}} = {id=>$edgeRef->{src}, members=>[], root=>{}};
          push @cNodes, $cNodes{$edgeRef->{src}};
          $cNodeAdj[$edgeRef->{src}] = [];
        }
        if (!defined $cNodeAdj[$edgeRef->{dest}]) {
          $cNodes{$edgeRef->{dest}} = {id=>$edgeRef->{dest}, members=>[], root=>{}};
          push @cNodes, $cNodes{$edgeRef->{dest}};
          $cNodeAdj[$edgeRef->{dest}] = [];
        }
        push @{$cNodeAdj[$edgeRef->{src}]}, $edgeRef;
      }
      my %visited = ();
      for(my $i=0; $i<scalar(@cNodes); ++$i) {
        if (!exists $visited{$cNodes[$i]->{id}} || 0==$visited{$cNodes[$i]->{id}}) {
          subGraphDFS($cNodes[$i]->{id}, \%cNodes, \@cNodeAdj, $cNodes[$i]->{id}, \%visited);
        }
      }
      %visited = (); # clear out
      # collect the root node(s) for cyclic graph
      my @rootCNodesIndex = ();
      for(my $i=0; $i<scalar(@cNodes); ++$i) {
        if (exists $cNodes[$i]->{root}->{$cNodes[$i]->{id}}) {
          push @rootCNodesIndex, $i;
        }
      }
      # let's reduce the duplication (reduce future work)
      my @nrRootCNodesIndex = ();
      my %visitedEdges = ();
      for(my $i=0; $i<scalar(@rootCNodesIndex); ++$i) {
        my $rootNodeRef = $cNodes[$rootCNodesIndex[$i]];
        my ($nodeIdInRCGraph, $edgeIdInRCGraph) = isRCGraph($rootNodeRef->{members}, \@cNodeAdj, \%visitedEdges);
        if (-1==$edgeIdInRCGraph) {
          $rootNodeRef->{edges} = [];
          my %namedMembers = (); grep { $namedMembers{$_}=1; } @{$rootNodeRef->{members}};
          foreach my $nodeId (@{$rootNodeRef->{members}}) {
            if (defined $cNodeAdj[$nodeId]) {
              foreach my $edgeRef (@{$cNodeAdj[$nodeId]}) {
                if (exists $namedMembers{$edgeRef->{src}} && exists $namedMembers{$edgeRef->{dest}}) {
                  push @{$rootNodeRef->{edges}}, $edgeRef;
                }
              }
            }
          }
          recordVisitedEdges($rootNodeRef->{edges}, \%visitedEdges);
          push @nrRootCNodesIndex, $rootCNodesIndex[$i];
        } else {
          # TODO : report $nodeIdInRCGraph, $edgeIdInRCGraph
          $rootNodeRef->{duplicated} = 1;
        }
      }
      %visitedEdges = (); # clear out
      # keep the results
      $graphRef->{cyclicGraphs} = [] if (!exists $graphRef->{cyclicGraphs});
      for(my $i=0; $i<scalar(@nrRootCNodesIndex); ++$i) {
        push @{$graphRef->{cyclicGraphs}}, $cNodes[$nrRootCNodesIndex[$i]];
      }
    }

    if (1==1) {
      my @aNodes = ();
      my %aNodes = ();
      my @aNodeAdj = ();
      foreach my $edgeRef (@acyclicEdges) {
        if (!defined $aNodeAdj[$edgeRef->{src}]) {
          $aNodes{$edgeRef->{src}} = {id=>$edgeRef->{src}, members=>[], root=>{}};
          push @aNodes, $aNodes{$edgeRef->{src}};
          $aNodeAdj[$edgeRef->{src}] = [];
        }
        if (!defined $aNodeAdj[$edgeRef->{dest}]) {
          $aNodes{$edgeRef->{dest}} = {id=>$edgeRef->{dest}, members=>[], root=>{}};
          push @aNodes, $aNodes{$edgeRef->{dest}};
          $aNodeAdj[$edgeRef->{dest}] = [];
        }
        push @{$aNodeAdj[$edgeRef->{src}]}, $edgeRef;
      }
      my %visited = ();
      for(my $i=0; $i<scalar(@aNodes); ++$i) {
        if (!exists $visited{$aNodes[$i]->{id}} || 0==$visited{$aNodes[$i]->{id}}) {
          subGraphDFS($aNodes[$i]->{id}, \%aNodes, \@aNodeAdj, $aNodes[$i]->{id}, \%visited);
        }
      }
      %visited = (); # clear out
      # collect the root node(s) for acyclic graph
      my @rootANodesIndex = ();
      for(my $i=0; $i<scalar(@aNodes); ++$i) {
        if (exists $aNodes[$i]->{root}->{$aNodes[$i]->{id}}) {
          push @rootANodesIndex, $i;
        }
      }
      # let's reduce the duplication (reduce future work)
      my @nrRootANodesIndex = ();
      my %visitedEdges = ();
      for(my $i=0; $i<scalar(@rootANodesIndex); ++$i) {
        my $rootNodeRef = $aNodes[$rootANodesIndex[$i]];
        my ($nodeIdInRCGraph, $edgeIdInRCGraph) = isRCGraph($rootNodeRef->{members}, \@aNodeAdj, \%visitedEdges);
        if (-1==$edgeIdInRCGraph) {
          $rootNodeRef->{edges} = [];
          my %namedMembers = (); grep { $namedMembers{$_}=1; } @{$rootNodeRef->{members}};
          foreach my $nodeId (@{$rootNodeRef->{members}}) {
            if (defined $aNodeAdj[$nodeId]) {
              foreach my $edgeRef (@{$aNodeAdj[$nodeId]}) {
                if (exists $namedMembers{$edgeRef->{src}} && exists $namedMembers{$edgeRef->{dest}}) {
                  push @{$rootNodeRef->{edges}}, $edgeRef;
                }
              }
            }
          }
          recordVisitedEdges($rootNodeRef->{edges}, \%visitedEdges);
          push @nrRootANodesIndex, $rootANodesIndex[$i];
        } else {
          # TODO: report $nodeIdInRCGraph, $edgeIdInRCGraph
          $rootNodeRef->{duplicated} = 1;
        }
      }
      %visitedEdges = (); # clear out
      # keep the results
      $graphRef->{acyclicGraphs} = [] if (!exists $graphRef->{acyclicGraphs});
      for(my $i=0; $i<scalar(@nrRootANodesIndex); ++$i) {
        push @{$graphRef->{acyclicGraphs}}, $aNodes[$nrRootANodesIndex[$i]];
      }
    }
  }

  # let's report what have been removed, and the remaining cycle
  printf "# partitionEdges - %d --> cg(%d) + ag(%d)\n", $subgraphId, 
    scalar(@{$graphRef->{cyclicGraphs}}), scalar(@{$graphRef->{acyclicGraphs}});

  my %nodesClass = ();
  my %edgesClass = ();
  for(my $i=0; $i<scalar(@{$graphRef->{cyclicGraphs}}); ++$i) {
    my $rootNodeRef = $graphRef->{cyclicGraphs}->[$i];
    printf "// cyclic %d.%d\n", $subgraphId, $i+1;
    foreach my $edgeRef (@{$rootNodeRef->{edges}}) {
      my $edgeId = $edgeRef->{id};
      $edgeId = '-'.$edgeId if ('-' eq $edgeRef->{strand});
      printf "\"%d\" -> \"%d\" [label=\"id %s\\l%s %d\"];\n", 
        $edgeRef->{src}, $edgeRef->{dest}, $edgeId, $edgeRef->{len}, $edgeRef->{cov};

      $nodesClass{$edgeRef->{src}} |= 2;
      $nodesClass{$edgeRef->{dest}} |= 2;
      $edgesClass{$edgeId} |= 2;
    }
    printf "// END - cyclic %d.%d\n", $subgraphId, $i+1;
  }

  for(my $i=0; $i<scalar(@{$graphRef->{acyclicGraphs}}); ++$i) {
    my $rootNodeRef = $graphRef->{acyclicGraphs}->[$i];
    printf "// acyclic %d.%d\n", $subgraphId, $i+1;
    foreach my $edgeRef (@{$rootNodeRef->{edges}}) {
      my $edgeId = $edgeRef->{id};
      $edgeId = '-'.$edgeId if ('-' eq $edgeRef->{strand});
      printf "\"%d\" -> \"%d\" [label=\"id %s\\l%s %d\",color=\"grey\",fontcolor=\"grey\"];\n", 
        $edgeRef->{src}, $edgeRef->{dest}, $edgeId, $edgeRef->{len}, $edgeRef->{cov};

      $nodesClass{$edgeRef->{src}} |= 1;
      $nodesClass{$edgeRef->{dest}} |= 1;
      $edgesClass{$edgeId} |= 1;
    }
    printf "// END - acyclic %d.%d\n", $subgraphId, $i+1;
  }
  my $nodesLen = scalar(keys %nodesClass);
  my $cnodesOnlyLen = 0;
  my $anodesOnlyLen = 0;
  my $sharedNodesLen = 0;
  foreach my $nodeClass (values %nodesClass) {
    if (1==$nodeClass) {
      $anodesOnlyLen++;
    } elsif (2==$nodeClass) {
      $cnodesOnlyLen++;
    } elsif (3==$nodeClass) {
      $sharedNodesLen++;
    }
  }
  my $edgesLen = scalar(keys %edgesClass);
  my $cedgesOnlyLen = 0;
  my $aedgesOnlyLen = 0;
  my $sharedEdgesLen = 0;
  foreach my $edgeClass (values %edgesClass) {
    if (1==$edgeClass) {
      $aedgesOnlyLen++;
    } elsif (2==$edgeClass) {
      $cedgesOnlyLen++;
    } elsif (3==$edgeClass) {
      $sharedEdgesLen++;
    }
  }
  printf "// tally - %d |V|=%d=c%d+a%d+s%d%s, |E|=%d=c%d+a%d%s\n",
    $subgraphId,
    $nodesLen, $cnodesOnlyLen, $anodesOnlyLen, $sharedNodesLen, 
    ($nodesLen==($cnodesOnlyLen+$anodesOnlyLen+$sharedNodesLen)) ? '' : ' [*delta='.($nodesLen-($cnodesOnlyLen+$anodesOnlyLen+$sharedNodesLen)).']',
    $edgesLen, $cedgesOnlyLen+$sharedEdgesLen, $aedgesOnlyLen+$sharedEdgesLen,
    ($edgesLen==($cedgesOnlyLen+$sharedEdgesLen+$aedgesOnlyLen+$sharedEdgesLen)) ? '' : ' [*delta='.($edgesLen-($cedgesOnlyLen+$sharedEdgesLen)-($aedgesOnlyLen+$sharedEdgesLen)).']';
}

sub collectDisjointGraphs {
  my ($graphRef, $verbose) = @_;
  # identify disjoint subgraphs via DFS
  # setup
  my @graphNodes = ();
  my %graphNodes = ();
  my @graphNodeAdj = ();
  foreach my $edgeRef (@{$graphRef->{edges}}) {
    # (src=>$src, dest=>$dest, id=>$id, strand=>$strand, len=>$len, cov=>$cov);
    # printf '"%d" -> "%d" [', $edgeRef->{src}, $edgeRef->{dest};
    if (!defined $graphNodeAdj[$edgeRef->{src}]) {
      $graphNodes{$edgeRef->{src}} = {id=>$edgeRef->{src}, members=>[], root=>{}};
      push @graphNodes, $graphNodes{$edgeRef->{src}};
      $graphNodeAdj[$edgeRef->{src}] = [];
    }
    if (!defined $graphNodeAdj[$edgeRef->{dest}]) {
      $graphNodes{$edgeRef->{dest}} = {id=>$edgeRef->{dest}, members=>[], root=>{}};
      push @graphNodes, $graphNodes{$edgeRef->{dest}};
      $graphNodeAdj[$edgeRef->{dest}] = [];
    }
    push @{$graphNodeAdj[$edgeRef->{src}]}, $edgeRef;
  }
  # traversal
  my %visited = ();
  @graphNodes = sort { 
    scalar(@{$graphNodeAdj[$b->{id}]}) <=> scalar(@{$graphNodeAdj[$a->{id}]}) || $a->{id} <=> $b->{id}
  } @graphNodes;
  for(my $i=0; $i<scalar(@graphNodes); ++$i) {
    if (!exists $visited{$graphNodes[$i]->{id}} || 0==$visited{$graphNodes[$i]->{id}}) {
      subGraphDFS($graphNodes[$i]->{id}, \%graphNodes, \@graphNodeAdj, $graphNodes[$i]->{id}, \%visited);
    }
  }
  %visited = (); # clear out

  # locate each root node of each disjoint graph
  my @rootNodesIndex = ();
  for(my $i=0; $i<scalar(@graphNodes); ++$i) {
    if (exists $graphNodes[$i]->{root}->{$graphNodes[$i]->{id}}) {
      push @rootNodesIndex, $i;
    }
  }
  # let's reduce the duplication (reduce future work)
  my @nrRootNodesIndex = ();
  my %visitedEdges = ();
  my $numDuplicates = 0;
  for(my $i=0; $i<scalar(@rootNodesIndex); ++$i) {
    my $rootNodeRef = $graphNodes[$rootNodesIndex[$i]];
    my ($nodeIdInRCGraph, $edgeIdInRCGraph) = isRCGraph($rootNodeRef->{members}, \@graphNodeAdj, \%visitedEdges);
    if (-1==$edgeIdInRCGraph) {
      $rootNodeRef->{edges} = [];
      my %namedMembers = (); grep { $namedMembers{$_}=1; } @{$rootNodeRef->{members}};
      foreach my $nodeId (@{$rootNodeRef->{members}}) {
        if (defined $graphNodeAdj[$nodeId]) {
          foreach my $edgeRef (@{$graphNodeAdj[$nodeId]}) {
            if (exists $namedMembers{$edgeRef->{src}} && exists $namedMembers{$edgeRef->{dest}}) {
              push @{$rootNodeRef->{edges}}, $edgeRef;
            }
          }
        }
      }
      recordVisitedEdges($rootNodeRef->{edges}, \%visitedEdges);
      push @nrRootNodesIndex, $rootNodesIndex[$i];
    } else {
      # TODO: report $nodeIdInRCGraph, $edgeIdInRCGraph
      if (0!=$verbose) {
        printf "DEBUG: Skipping root node %d, node %d edge %d visisted previously\n", 
          $rootNodeRef->{id}, $nodeIdInRCGraph, $edgeIdInRCGraph;
      }
      $rootNodeRef->{duplicated} = 1;
      $numDuplicates++;
    }
  }
  %visitedEdges = (); # clear out
  my $numRootNodes = scalar(@rootNodesIndex);
  printf STDERR "%d disjoint graphs detected.\n", $numRootNodes;
  printf STDERR "%d (%.2f%%) duplicated disjoint graphs.\n", $numDuplicates, $numDuplicates*100.0/(0==$numRootNodes?1:$numRootNodes);

  ##### FIXME: DEBUG - find where is the missing edges
  my %wantedEdges = ();
  if ('' ne $debugEdge) {
    grep { $wantedEdges{$_} = $_; } split(/,/, $debugEdge);
    for(my $i=0; $i<scalar(@nrRootNodesIndex); ++$i) {
      my $rootNodeRef = $graphNodes[$nrRootNodesIndex[$i]];
      foreach my $edgeRef (@{$rootNodeRef->{edges}}) {
        if (exists $wantedEdges{$edgeRef->{id}}) {
          printf "DEBUG: after NR, edge %s%d found in rootNode($nrRootNodesIndex[$i])\n", $edgeRef->{strand}, $edgeRef->{id};
        }
      }
    }
  }
  ##### FIXME: DEBUG - find where is the missing edges
  

  # TODO: let's further partition each disjoint graph if possible
  $graphRef->{cyclicGraphs} = {type=>'cyclic', rootnodes=>[], edges=>[]};
  $graphRef->{acyclicGraphs} = {type=>'acyclic', rootnodes=>[], edges=>[]};
  $graphRef->{disjointcyclicGraphs} = {type=>'disjointcyclic', rootnodes=>[]};
  my $subgraphCount = 0;
  for(my $i=0; $i<scalar(@nrRootNodesIndex); ++$i) {
    my $rootNodeRef = $graphNodes[$nrRootNodesIndex[$i]];

    ##### FIXME: DEBUG - find where is the missing edges
    if ('' ne $debugEdge) {
      my $found = 0;
      foreach my $edgeRef (@{$rootNodeRef->{edges}}) {
        if (exists $wantedEdges{$edgeRef->{id}}) {
          printf "DEBUG: partition, edge %s%d found in rootNode($nrRootNodesIndex[$i])\n", $edgeRef->{strand}, $edgeRef->{id};
          $found = 1;
        }
      }
      if (0==$found) {
          print "DEBUG: partition, edge $debugEdge NOT found in rootNode($nrRootNodesIndex[$i])\n";
      }
    }
    ##### FIXME: DEBUG - find where is the missing edges

    if (scalar(@{$rootNodeRef->{edges}})>0) {
      $subgraphCount++;
      regroupDisjointGraph($rootNodeRef, $subgraphCount, $verbose);

      for(my $j=0; $j<scalar(@{$rootNodeRef->{cyclicGraphs}}); ++$j) {
        my $subrootNodeRef = $rootNodeRef->{cyclicGraphs}->[$j];
        $subrootNodeRef->{disjoinGraphRef} = $rootNodeRef;
        push @{$graphRef->{cyclicGraphs}->{rootnodes}}, $subrootNodeRef;
        foreach my $edgeRef (@{$subrootNodeRef->{edges}}) {
          $edgeRef->{class} = 'cyclic';
          push @{$graphRef->{cyclicGraphs}->{edges}}, $edgeRef;
        }
        ##### FIXME: DEBUG - find where is the missing edges
        if ('' ne $debugEdge) {
          foreach my $edgeRef (@{$subrootNodeRef->{edges}}) {
            if (exists $wantedEdges{$edgeRef->{id}}) {
              print "DEBUG: partition+cyclic, edge $edgeRef->{strand}$edgeRef->{id} found in rootNode($nrRootNodesIndex[$i]) cyclic#$j\n";
            }
          }
        }
        ##### FIXME: DEBUG - find where is the missing edges
      }
      for(my $j=0; $j<scalar(@{$rootNodeRef->{acyclicGraphs}}); ++$j) {
        my $subrootNodeRef = $rootNodeRef->{acyclicGraphs}->[$j];
        $subrootNodeRef->{disjoinGraphRef} = $rootNodeRef;
        push @{$graphRef->{acyclicGraphs}->{rootnodes}}, $subrootNodeRef;
        foreach my $edgeRef (@{$subrootNodeRef->{edges}}) {
          $edgeRef->{class} = 'acyclic';
          push @{$graphRef->{acyclicGraphs}->{edges}}, $edgeRef;
        }
        ##### FIXME: DEBUG - find where is the missing edges
        if ('' ne $debugEdge) {
          foreach my $edgeRef (@{$subrootNodeRef->{edges}}) {
            if (exists $wantedEdges{$edgeRef->{id}}) {
              print "DEBUG: partition+acyclic, edge $edgeRef->{strand}$edgeRef->{id} found in rootNode($nrRootNodesIndex[$i]) acyclic#$j\n";
            }
          }
        }
        ##### FIXME: DEBUG - find where is the missing edges
      }

      $rootNodeRef->{disjointId} = $subgraphCount;
      if (scalar(@{$rootNodeRef->{cyclicGraphs}})>0) {
        push @{$graphRef->{disjointcyclicGraphs}->{rootnodes}}, $rootNodeRef;
      }
    } else {
      printf STDERR "WARNING: Node %d has no edge.\n", $rootNodeRef->{id};
    }
  }
}

sub graphBamFile {
  my ($bamFile, $graphRef, $graphEdgeType, $cov2cn, $file) = @_;

  # set up the data for working with bam file!
  # edges list; each edge has its subgraph_Id, edge_id, amplClass
  # edge should have auxiliary info: prevEdges, nextEdges
  my %edgesInfo = ();
  foreach my $rootNodeRef (@{$graphRef->{rootnodes}}) {
    my %verticesOut = ();
    my %verticesIn = ();
    foreach my $edgeRef (@{$rootNodeRef->{edges}}) {
      # to record subgraph ID [ZG:i:]
      # to record edge [ZE:Z:0000_0000]
      # to record amplifcation classification [ZA:Z:]
      my $edgeId = sprintf("%s_%d", $graphEdgeType, $edgeRef->{id});
      my %item = (id=>$edgeId, edgeRef=>$edgeRef, rootNodeRef=>$rootNodeRef, tags=>{});
      $item{tags}->{ZG} = sprintf("ZG:Z:%s", $rootNodeRef->{subgraphId});
      my $uedgeId = sprintf("%04d_%04d", $rootNodeRef->{subgraphId}, $edgeRef->{id});
      $item{uid} = $uedgeId;
      $item{tags}->{ZE} = sprintf("ZE:Z:%s", $uedgeId);
      my $covBin = getCoverageBin($edgeRef->{cov}, $cov2cn);
      $item{tags}->{ZA} = sprintf("ZA:Z:%s", $covBin);
      my $cn = 2 * int($edgeRef->{cov}) / $cov2cn;
      $item{tags}->{ZD} = sprintf("ZD:f:%.1f", $cn);
      $edgesInfo{$edgeId} = \%item;

      $verticesOut{$edgeRef->{src}} = {} if (!exists $verticesOut{$edgeRef->{src}});
      $verticesOut{$edgeRef->{src}}->{$uedgeId}++;
      $verticesIn{$edgeRef->{dest}} = {} if (!exists $verticesIn{$edgeRef->{dest}});
      $verticesIn{$edgeRef->{dest}}->{$uedgeId}++;
    }
    # to record prevEdge(s) and nextEdge(s) [ZP:Z:],  [ZN:Z:]
    foreach my $edgeRef (@{$rootNodeRef->{edges}}) {
      my $edgeId = sprintf("%s_%d", $graphEdgeType, $edgeRef->{id});
      my $edgeInfoRef = $edgesInfo{$edgeId};

      my @bits = ();
      my $selfLoop = 0;
      if (exists $verticesIn{$edgeRef->{src}}) {
        while (my ($uedgeId, $count) = each %{$verticesIn{$edgeRef->{src}}}) {
          if ($edgeInfoRef->{uid} eq $uedgeId) {
            $selfLoop++;
            push @bits, 'SL'.$uedgeId;
          } else {
            push @bits, $uedgeId;
          }
        }
      }
      $edgeInfoRef->{tags}->{ZP} = sprintf("ZP:Z:%d%s;%s", scalar(@bits), (0==$selfLoop)?'':'wSL', join(";", @bits));

      @bits = ();
      $selfLoop = 0;
      if (exists $verticesOut{$edgeRef->{dest}}) {
        while (my ($uedgeId, $count) = each %{$verticesOut{$edgeRef->{dest}}}) {
          if ($edgeInfoRef->{uid} eq $uedgeId) {
            $selfLoop++;
            push @bits, 'SL'.$uedgeId;
          } else {
            push @bits, $uedgeId;
          }
        }
      }
      $edgeInfoRef->{tags}->{ZN} = sprintf("ZN:Z:%d%s;%s", scalar(@bits), (0==$selfLoop)?'':'wSL', join(";", @bits));
    }
  }

  # stream bam headers
  # read bam record and inject tags
  # index new bam file
  open INFILE, "samtools view -h $bamFile | " || die "Fail to open $bamFile\n$!\n";
  open OUTFILE, "|samtools view -b - >  $file" || die "Fail to open $file\n$!\n";
  while (<INFILE>) {
    if (/^@/) {
      print OUTFILE $_;
    } else {
      chomp();
      my @bits = split(/\t/);
      if (exists $edgesInfo{$bits[0]}) {
        # to inject
        my $edgeRef = $edgesInfo{$bits[0]};
        my %tags = ();
        while (my ($tag, $value) = each %{$edgeRef->{tags}}) {
          $tags{$tag} = $value;
        }
        for(my $i=11; $i<scalar(@bits); ++$i) {
          my ($tag, $tagType, $tagValue) = split(/:/,$bits[$i],3);
          if (exists $tags{$tag}) {
            $bits[$i] = $tags{$tag};
            delete $tags{$tag}
          }
        }
        while (my ($tag, $value) = each %{$edgeRef->{tags}}) {
          push @bits, $value;
        }
        print OUTFILE join("\t", @bits), "\n";
      } else {
        print OUTFILE $_,"\n";
      }
    }
  }
  close OUTFILE;
  close INFILE;
  system("samtools index $file");
}

sub sieveGraph {
  my $G_USAGE = "
$0 sievegraph --gv <gvFile> --diploidcov <diploidCoverage> --oprefix <outputPrefix>
  --gv STR            assembly graphs (gv format)
  --diploidcov INT    average diploid coverage
  --oprefix STR       prefix to output filename
";

  my $gvFile = undef;
  my $cov = undef;
  my $prefix = undef;
  my $bamFile = undef;
  my $minlen = 0; # min total edges length (kb) for an amp class
	my $verbose = 0;
	my $help = 0;
  ## my $debugEdge = "";
	
	GetOptions (
	"gv=s"         => \$gvFile,
	"diploidcov=i" => \$cov,
	"bam=s"        => \$bamFile,
	"oprefix=s"    => \$prefix,
	"minlenkb=i"     => \$minlen,
	"verbose!"     => \$verbose,
	"edge=s"       => \$debugEdge,
	"help!"        => \$help)
	or die("Error in command line arguments\n$G_USAGE");
	
	die "$G_USAGE" if ($help);

  my %graph = ();
  readGVGraphFile($gvFile, \%graph);

  # remove the edge failing the coverage
  removeEdgesInsufficientCoverage(\%graph, $cov, $verbose);

  # clean up the graph
  collectDisjointGraphs(\%graph, $verbose);

  my $oFile = '';
  my @comments = ();
  #
  $oFile = sprintf("%s.cn%d.disjointcyclic.gv", $prefix, $cov);
  @comments = ();
  push @comments, "disjointcyclic subgraphs";
  push @comments, sprintf("file = %s", $oFile);
  push @comments, sprintf("|subgraphs| = %d", getUniqueSubgraphs($graph{disjointcyclicGraphs}->{rootnodes}));
  push @comments, sprintf("data = %s", $gvFile);
  push @comments, sprintf("cn >= %d", $cov);
  writeGraph($graph{disjointcyclicGraphs}, $graph{edgeType}, \@comments, $cov, $minlen, $oFile);
  # 
  $oFile = sprintf("%s.cn%d.disjointcyclic.readids", $prefix, $cov);
  writeReadIds($graph{disjointcyclicGraphs}, $graph{edgeType}, $oFile);

  # optional; generate a new bam file!
  if (defined $bamFile) {
    $oFile = sprintf("%s.cn%d.disjointcyclic.bam", $prefix, $cov);
    graphBamFile($bamFile, $graph{disjointcyclicGraphs}, $graph{edgeType}, $cov, $oFile);
  }
}

sub initSubgraph {
  my $G_USAGE = "
$0 initsubgraph --summary <summaryFile>
  --summary STR       summary file
";

  my $summaryFile = undef;
  my $bamFile = undef;
	my $verbose = 0;
	my $help = 0;
  my $minKb = 25;
	
	GetOptions (
	"summary=s"    => \$summaryFile,
	"bam=s"        => \$bamFile,
	"verbose!"     => \$verbose,
	"help!"        => \$help)
	or die("Error in command line arguments\n$G_USAGE");
	
	die "$G_USAGE" if ($help);

  # read each line and set up the data file
  my @putatives = ();
  my $numSubgraphs = 0;
  open INFILE, $summaryFile || die "Fail to open $summaryFile\n$!\n";
  my $colsLine = '';
  my @cols = ();
  my $numCols = 0;
  my $outputPrefix = 'noname';
  while (<INFILE>) {
    if (/^\/\//) {
      # // file = /net/nwgc/vol1/techdev/GBMTW/data/B168/B168.fq.gz-ecLegoV3/B168.fq.gz_candidate_assembly.ci5000/B168.cn37.disjointcyclic.gv
      if (/^\/\/\s*file\s*=\s*/) {
        $outputPrefix = $';
        chomp($outputPrefix);
      }
    } else {
      if (0==$numCols) {
        chomp();
        @cols = split(/\t/);
        $numCols = scalar(@cols);
      } else {
        chomp();
        my @values = split(/\t/);
        my %item = ();
        for(my $i=0; $i<$numCols; ++$i) {
          if (defined $values[$i]) {
            $item{$cols[$i]} = $values[$i];
          }
        }
        $numSubgraphs++;
        if ($item{kbp_1_low_amp}>=$minKb || $item{kbp_2_med_amp}>=$minKb || $item{kbp_3_high_amp}>=$minKb) {
          # generate the script to refine a subgraph!
          # TODO: extract out the modules and cluster jobs nature
          my @size = ({type=>'kbp_1_low_amp',value=>$item{kbp_1_low_amp}}
            ,{type=>'kbp_2_med_amp',value=>$item{kbp_2_med_amp}}
            ,{type=>'kbp_3_high_amp',value=>$item{kbp_3_high_amp}});
          @size = sort { $b->{value}<=>$a->{value} } @size;
          $item{'_sort_value'} = $size[0]->{value};
          $item{'_sort_type'} = $size[0]->{type};
          push @putatives, \%item;
        }
      }
    }
  }
  close INFILE;

  # sort and start with the largest (then higher copy number)
  @putatives = sort { $b->{_sort_value}<=>$a->{_sort_value} } @putatives;
  print "#\n";
  printf "# %d/%d of the subgraphs contain putative circular genomes.\n", scalar(@putatives), $numSubgraphs;
  my $fastq = 'noname';
  my $fastqDir = '/';
  my $sample = 'noname';
  if ($outputPrefix=~/\-ecLegoV3\//) {
    $fastq = $`;
    my @bits = split(/\//, $fastq);
    pop @bits;
    $fastqDir = join("\/", @bits);
    if (scalar(@bits)>=2) {
      $sample = $bits[$#bits];
    }
  }
  my $numPutatives = scalar(@putatives);
  if ($numPutatives>0) {
    printf "cd %s \n", $fastqDir;
  }
  for(my $i=0; $i<$numPutatives; ++$i) {
    my $subgraphRef = $putatives[$i];
    printf "# %d. subgraph_%d (%sL: %s kbp, %sM: %s kbp, %sH: %s kbp)\n", 
      ($i+1), $subgraphRef->{subgraphId}, 
      ($subgraphRef->{kbp_1_low_amp}>=$minKb) ? '*' : '', $subgraphRef->{kbp_1_low_amp}, 
      ($subgraphRef->{kbp_2_med_amp}>=$minKb) ? '*' : '', $subgraphRef->{kbp_2_med_amp}, 
      ($subgraphRef->{kbp_3_high_amp}>=$minKb) ? '*' : '', $subgraphRef->{kbp_3_high_amp};
    print "/net/nwgc/vol1/techdev/GBMTW/pipeline/ecLegoV3_s01.refine.sh \\\n";
    print "$fastq \\\n";
    print "/net/nwgc/vol1/techdev/GBMTW/pipeline/genomes/t2tv2.fasta  \\\n";
    printf "%d 2>&1 | tee %s.ecLegoV3.%d_%d_subgraph_%d.run.log\n", 
      $subgraphRef->{subgraphId}, $sample, ($i+1), $numPutatives, $subgraphRef->{subgraphId};
    print "#\n";
  }
}

sub growGraphFromEdge {
  my ($graphRef, $edge, $idonly, $verbose) = @_;
  # identify disjoint subgraphs via DFS
  # setup
  my @graphNodes = ();
  my %graphNodes = ();
  my @graphNodeAdj = ();
  foreach my $edgeRef (@{$graphRef->{edges}}) {
    # (src=>$src, dest=>$dest, id=>$id, strand=>$strand, len=>$len, cov=>$cov);
    # printf '"%d" -> "%d" [', $edgeRef->{src}, $edgeRef->{dest};
    if (!defined $graphNodeAdj[$edgeRef->{src}]) {
      $graphNodes{$edgeRef->{src}} = {id=>$edgeRef->{src}, members=>[], root=>{}};
      push @graphNodes, $graphNodes{$edgeRef->{src}};
      $graphNodeAdj[$edgeRef->{src}] = [];
    }
    if (!defined $graphNodeAdj[$edgeRef->{dest}]) {
      $graphNodes{$edgeRef->{dest}} = {id=>$edgeRef->{dest}, members=>[], root=>{}};
      push @graphNodes, $graphNodes{$edgeRef->{dest}};
      $graphNodeAdj[$edgeRef->{dest}] = [];
    }
    push @{$graphNodeAdj[$edgeRef->{src}]}, $edgeRef;
  }
  # traversal
  my %visited = ();
  # find the nodes related to the edges and restrict to these traversal only!!!
  ### my @edgeIds = split(/,/, $edge);
  ### my %edgeIds = ();
  ### grep { $edgeIds{$_} = $_; } @edgeIds;
  my %wantedEdges = ();
  if ('' ne $debugEdge) {
    grep { $wantedEdges{$_} = $_; } split(/,/, $edge);
  }
  @graphNodes = ();
  my %nodesSetup = ();
  if ('' ne $debugEdge) {
    foreach my $edgeRef (@{$graphRef->{edges}}) {
      if (exists $wantedEdges{$edgeRef->{id}}) {
        if (!exists $nodesSetup{$edgeRef->{src}}) {
          push @graphNodes, $graphNodes{$edgeRef->{src}};
        }
        if (!exists $nodesSetup{$edgeRef->{dest}}) {
          push @graphNodes, $graphNodes{$edgeRef->{dest}};
        }
      }
    }
  }
  @graphNodes = sort { 
    scalar(@{$graphNodeAdj[$b->{id}]}) <=> scalar(@{$graphNodeAdj[$a->{id}]}) || $a->{id} <=> $b->{id}
  } @graphNodes;
  for(my $i=0; $i<scalar(@graphNodes); ++$i) {
    if (!exists $visited{$graphNodes[$i]->{id}} || 0==$visited{$graphNodes[$i]->{id}}) {
      subGraphDFS($graphNodes[$i]->{id}, \%graphNodes, \@graphNodeAdj, $graphNodes[$i]->{id}, \%visited);
    }
  }
  %visited = (); # clear out

  # locate each root node of each disjoint graph
  my @rootNodesIndex = ();
  for(my $i=0; $i<scalar(@graphNodes); ++$i) {
    if (exists $graphNodes[$i]->{root}->{$graphNodes[$i]->{id}}) {
      push @rootNodesIndex, $i;
    }
  }

  # just report everything... no attempt for RC & other logics
  # FIXME: Temp hack
  my @nrRootNodesIndex = ();
  my %visitedEdges = ();
  my $numDuplicates = 0;
  for(my $i=0; $i<scalar(@rootNodesIndex); ++$i) {
    my $rootNodeRef = $graphNodes[$rootNodesIndex[$i]];
    my ($nodeIdInRCGraph, $edgeIdInRCGraph) = isRCGraph($rootNodeRef->{members}, \@graphNodeAdj, \%visitedEdges);
    if (-1==$edgeIdInRCGraph) {
      $rootNodeRef->{edges} = [];
      my %namedMembers = (); grep { $namedMembers{$_}=1; } @{$rootNodeRef->{members}};
      foreach my $nodeId (@{$rootNodeRef->{members}}) {
        if (defined $graphNodeAdj[$nodeId]) {
          foreach my $edgeRef (@{$graphNodeAdj[$nodeId]}) {
            if (exists $namedMembers{$edgeRef->{src}} && exists $namedMembers{$edgeRef->{dest}}) {
              push @{$rootNodeRef->{edges}}, $edgeRef;
            }
          }
        }
      }
      recordVisitedEdges($rootNodeRef->{edges}, \%visitedEdges);
      push @nrRootNodesIndex, $rootNodesIndex[$i];
    } else {
      # TODO: report $nodeIdInRCGraph, $edgeIdInRCGraph
      if (0!=$verbose) {
        printf "DEBUG: Skipping root node %d, node %d edge %d visisted previously\n", 
          $rootNodeRef->{id}, $nodeIdInRCGraph, $edgeIdInRCGraph;
      }
      $rootNodeRef->{duplicated} = 1;
      $numDuplicates++;
    }
  }
  %visitedEdges = (); # clear out

  # FIXME: too much ... let's just dump out the results for now!
  if (0==$idonly) {
    my $subgraphCount = 0;
    for(my $i=0; $i<scalar(@nrRootNodesIndex); ++$i) {
      my $rootNodeRef = $graphNodes[$nrRootNodesIndex[$i]];
      if (scalar(@{$rootNodeRef->{edges}})>0) {
        $subgraphCount++;
        # let's report plainly
        printf "# START : sG#%d (%d edges)\n", $subgraphCount, scalar(@{$rootNodeRef->{edges}});
        foreach my $edgeRef (@{$rootNodeRef->{edges}}) {
          printf "\"%d\" -> \"%d\" [label=\"id %s%d\\l%s %s\",penwidth=3];\n",
            $edgeRef->{src}, $edgeRef->{dest},
            $edgeRef->{strand}, $edgeRef->{id}, 
            $edgeRef->{len}, $edgeRef->{cov};
        }
        printf "# END : --- sG#%d (%d edges)\n", $subgraphCount, scalar(@{$rootNodeRef->{edges}});
      }
    }
  } else {
    my %seenEdges = (); # clear out
    for(my $i=0; $i<scalar(@nrRootNodesIndex); ++$i) {
      my $rootNodeRef = $graphNodes[$nrRootNodesIndex[$i]];
      if (scalar(@{$rootNodeRef->{edges}})>0) {
        foreach my $edgeRef (@{$rootNodeRef->{edges}}) {
          if (!exists $seenEdges{$edgeRef->{id}}) {
            $seenEdges{$edgeRef->{id}} = $edgeRef;
          }
        }
      }
    }
    foreach my $edgeId (sort { $a <=> $b} keys %seenEdges) {
      printf "edge_%d\n", $edgeId;
    }
  }
}

sub spanEdge {
  my $G_USAGE = "
$0 spanedge --gv <gvFile> --diploidcov <diploidCoverage> --oprefix <outputPrefix>
  --gv STR            assembly graphs (gv format)
  --diploidcov INT    average diploid coverage
  --edge STR          edge to seed
";

  my $gvFile = undef;
  my $cov = undef;
  my $verbose = 0;
  my $idonly = 0;
	my $help = 0;
	
	GetOptions (
	"gv=s"         => \$gvFile,
	"diploidcov=i" => \$cov,
	"verbose!"     => \$verbose,
	"edge=s"       => \$debugEdge,
	"idonly!"      => \$idonly,
	"help!"        => \$help)
	or die("Error in command line arguments\n$G_USAGE");
	
	die "$G_USAGE" if ($help);

  my %graph = ();
  readGVGraphFile($gvFile, \%graph);

  # remove the edge failing the coverage
  removeEdgesInsufficientCoverage(\%graph, $cov, $verbose) if ($cov>0);

  # collect only related edges
  growGraphFromEdge(\%graph, $debugEdge, $idonly, $verbose);
}

