set boxwidth 0.9 absolute
set style fill   solid 1.00 border lt -1
set key fixed right top vertical Right noreverse noenhanced autotitle nobox
set style increment default
set style histogram clustered gap 1 title textcolor lt -1
set style data histograms
set yrange [10:62]

p for [i=0:4] sprintf("%03.0f/profile.dat", i) u 2 every ::2 ti sprintf("%03.0f", i)
