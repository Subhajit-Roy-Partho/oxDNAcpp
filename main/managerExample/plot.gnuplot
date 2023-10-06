set term png size 1200,900
set xlabel 'Time(SU)'
set logscale x
set yrange [0:-1.4]
set ylabel 'Energy(SU)'
set title '4flat'
set output './plots/4flat-r4.png'
plot '/scratch/sroy85/seeding/4flat/replica-4/T-0.122/energy.dat' w l title 'T-0.122',\
	'/scratch/sroy85/seeding/4flat/replica-4/T-0.123/energy.dat' w l title 'T-0.123',\
	'/scratch/sroy85/seeding/4flat/replica-4/T-0.124/energy.dat' w l title 'T-0.124',\
	'/scratch/sroy85/seeding/4flat/replica-4/T-0.125/energy.dat' w l title 'T-0.125',\
	'/scratch/sroy85/seeding/4flat/replica-4/T-0.126/energy.dat' w l title 'T-0.126',\
	'/scratch/sroy85/seeding/4flat/replica-4/T-0.127/energy.dat' w l title 'T-0.127',\
	'/scratch/sroy85/seeding/4flat/replica-4/T-0.128/energy.dat' w l title 'T-0.128',\
	'/scratch/sroy85/seeding/4flat/replica-4/T-0.129/energy.dat' w l title 'T-0.129'
