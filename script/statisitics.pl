#! usr/bin/perl
open LIST, '<statistic.txt' or die "Cannot read from statistic.txt:$!";
open STDOUT, '>RESULT.txt' or die "Cannot write to RESULT.txt:$!";
chomp(my @list = <LIST>);
my %count;
my $name;
foreach $name (@list) {
    $count{$name} += 1;
}
foreach $name (sort keys %count) {
    print "$name\t$count{$name}\n";
}
