
#rm *.png
#rm frames/*.png
#rm movie2.*.pov 
f_n="10" #file name doesn't work if more than 2 characters long...fix maybe?
t1=0.0937 #radius of hole from original pov file
t2=0.0852 #radius of electron for original pov file
npi=156
#step=1000

#steps=`ls -l movie.*.pov | wc | awk '{print $1}'`
steps=1
#this copies all movie.*.pov to movie2.*.pov
#ls -l movie.*.pov | awk -v n=1 '{print "cp",$9,"movie2."n".pov" ; n++}' > temp.sh
#bash temp.sh
#i replaced "movie2.$s"".pov" with please.pov 
for((s = 1 ; s <= $steps ; s++)) 
do
    echo $s
    #-F makes the field separator what you choose (in this case, a comma)
    #shit1 labels hole density lines, shit2 labels electron density lines
    #c1, c2 +1 after each of those lines
    #prints the rest of the file exactly the same
    awk -v c1=0 -v c2=0 -v t1=$t1 -v t2=$t2 -F "," '{if ($4==t1) {print "shit1,"c1","$1","$2","$3","$4","$5","$6","$7","$8; c1++} else {if ($4==t2){print "shit2,"c2","$1","$2","$3","$4","$5","$6","$7","$8 ; c2++} else print}}' $f_n.pov > ptemp1
    
    for ((l = 0 ; l < $npi ; l++)) 
    do 
	cp ptemp1 ptemp2 
	n1=`awk -v n=$l '{if ($2==n) print $5}' atemp_$s` #n1 = hole density
	n2=`awk -v n=$l '{if ($2==n) print $6}' atemp_$s` #n2 = electron density
	echo $n1 $n2
	echo $n1 $n2 > temp
	test=`awk -v n=0 '{if ($1>$2) n=1}END{print n}' temp` #test = 1 if hole density is larger than electron density

	if [ $test -eq 1 ] 
	then
	    size=1
	    awk -v l=$l -v n1=$n1 -v s=$size -F "," '{if ($1=="shit1" && $2==l) print $3","$4","$5","s*n1","$7","$8","$9",0.500>)" ; else print}' ptemp2 > ptemp3
	    awk -v l=$l -v n2=$n2 -v s=$size -F "," '{if ($1=="shit2" && $2==l) print $3","$4","$5","s*n2","$7","$8","$9",0.000>)" ; else print}' ptemp3 > ptemp1
	else
	    awk -v l=$l -v n1=$n1 -v s=$size -F "," '{if ($1=="shit1" && $2==l) print $3","$4","$5","s*n1","$7","$8","$9",0.000>)" ; else print}' ptemp2 > ptemp3
	    awk -v l=$l -v n2=$n2 -v s=$size -F "," '{if ($1=="shit2" && $2==l) print $3","$4","$5","s*n2","$7","$8","$9",0.500>)" ; else print}' ptemp3 > ptemp1

	fi
    done
    
    cp ptemp1 $f_n.pov
    povray +W360 +H292 -Display=0 +UA $f_n.pov
    #povray  +W360 +H566 -Display=0 +UA state23-new.pov
    #for ((k=0 ; k < 1 ; k++)) 
    #do
	#cp temp.png mtemp_$step.png
	#step=$((step+1))
    #done
    #mv temp.png frames/frame_$s.png
done

#mencoder mf://*.png -mf fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o output.avi -ffourcc DX50

#cp output.avi frames/.

rm ptemp*
cp $f_n.pov $f_n.pov.temp
head -20 gap.dat > headgap
cat $f_n.pov.temp qcff_parm.inpt start.inpt headgap > $f_n.pov
rm $f_n.pov.temp
