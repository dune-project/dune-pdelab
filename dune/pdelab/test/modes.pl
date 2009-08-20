#/usr/bin/perl

my $max = 5e9;
my $binsize = 2.9966e+07/8;
my $sigma = $binsize;

my @freq;

for(my $n = 1; ; ++$n) {
    my $m;
    for($m = 1; ; ++$m) {
        my $p;
        for($p = 0; ; ++$p) {
            my $f = 3e8*sqrt($n**2+$m**2+$p**2)/2;
            last if $f > $max;
            push @freq, $f;
        }
        last if $p == 0;
    }
    last if $m == 1;
}

for(my $pos = $binsize / 2; $pos < $max; $pos += $binsize) {
    my $amp = 0;
    for my $freq (@freq) {
        $amp += exp(-($pos-$freq)**2/2/$sigma**2);
    }
    print "$pos\t$amp\n";
}
