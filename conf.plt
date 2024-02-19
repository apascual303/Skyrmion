set pm3d map
set palette rgb 35,13,10
set xrange[0:80]
set yrange[0:80]
set size ratio 1
set cbrange[-1:1]
do for[n=0:100]{
	splot"WT_".n.".txt" matrix
}