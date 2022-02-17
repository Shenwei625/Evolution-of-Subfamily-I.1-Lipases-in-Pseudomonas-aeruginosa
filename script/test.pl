#! usr/bin/perl
open LIST, '<statistic.txt' or die "Cannot write to ls.out:$!";
open STDOUT, '>RESULT.txt' or die "Cannot write to ls.out:$!";
chomp(my @list = <LIST>);
my %count;
my $name;
foreach $name (@list) {
    $count{$name} += 1;
}
foreach $name (sort keys %count) {
    print "$name\t$count{$name}\n";
}
