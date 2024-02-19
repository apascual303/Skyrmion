set pm3d map
set xrange[0:80]
set yrange[0:80]
set size ratio 1
set cbrange[-1:1]
do for[n=0:100]{
	splot"alpha_".n.".txt" matrix
}