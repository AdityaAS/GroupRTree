#path1=../Dataset/
path1=/Users/adityasharma/PersonalProjects/ADAPTREPLICA/datasets
GMINENTRIES=4
GMAXENTRIES=8
GAUXMINENTRIES=4
GAUXMAXENTRIES=8
RMINENTRIES=4
RMAXENTRIES=8
Epsilon=$2
MINPOINTS=$3
UNDEFINED=10000000

echo "=========================================="
echo "Following code was profiled on : `date`"
echo "Script for Profiling serial Code using Vampir"
echo "\nFiles given->\n\t$1\n\t$2\n\t$3"
echo "=========================================="
echo
echo "\tCompiling Code..."
make clean
make -f $4
if [ $? -eq 1 ] ;	then
	echo "\tError occured during compilation\nTerminating script\n"
	exit 1
else
	echo "\n\tCompilation done successfully.\n"
fi
#executable file name
exename=output
#temp=./temp
#rm -rf $temp
#mkdir $temp

#export VT_ON=no
#export VT_UNIFY=no
#export VT_MAX_FLUSHES=10
#export VT_BUFFER_SIZE=256M
export VT_MODE=STAT
#export VT_PFORM_GDIR=$temp
#export VT_PFORM_LDIR=$temp
#export VT_FILE_PREFIX=$1_$2
#nm output>output.nm
#export VT_GNU_NMFILE=output.nm
echo "\t\tRunning Code..."

#$1 is input file and $2 is output file to which the program output is to be written
mpirun -hostfile hostlist1 -np $5 ./output $path1$1 output_$1\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS\_nodes=$5\_threads=$6.txt $GMINENTRIES $GMAXENTRIES $GAUXMINENTRIES $GAUXMAXENTRIES $RMINENTRIES $RMAXENTRIES $Epsilon $MINPOINTS $UNDEFINED cluster_$1\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS\_nodes=$5\_threads=$6.txt $6
#mpirun -pernode ./output $path1$1 output_$1\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS\_nodes=$5\_threads=$6.txt $GMINENTRIES $GMAXENTRIES $GAUXMINENTRIES $GAUXMAXENTRIES $RMINENTRIES $RMAXENTRIES $Epsilon $MINPOINTS $UNDEFINED cluster_$1\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS\_nodes=$5\_threads=$6.txt $6
#mpirun -np $5 valgrind --tool=memcheck --leak-check=full ./output $path1$1 output_$1\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS\_nodes=$5\_threads=$6.txt $GMINENTRIES $GMAXENTRIES $GAUXMINENTRIES $GAUXMAXENTRIES $RMINENTRIES $RMAXENTRIES $Epsilon $MINPOINTS $UNDEFINED cluster_$1\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS\_nodes=$5\_threads=$6.txt $6
#valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all ./output $path1$1 output_$1\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS.txt $GMINENTRIES $GMAXENTRIES $GAUXMINENTRIES $GAUXMAXENTRIES $RMINENTRIES $RMAXENTRIES $CELLSIZE $Epsilon $MINPOINTS $UNDEFINED cluster_$1\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS.txt

#valgrind --tool=memcheck --leak-check=full ./output $path1$1 output_$1\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS.txt $GMINENTRIES $GMAXENTRIES $GAUXMINENTRIES $GAUXMAXENTRIES $RMINENTRIES $RMAXENTRIES $CELLSIZE $Epsilon $MINPOINTS $UNDEFINED cluster_$1\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS.txt

#gprof -z output > gprof_dataset-$2\_m=$m\_M=$M\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS.txt 

if [ $? -ne 0 ] ;	then
	echo "\tError occured. Terminating script...\n"
	exit 1
fi

echo "Exeuction Done on : `date`"
#cd $temp
#vtunify --stats $temp/$1_$2
#cd ..
#cp $temp/$1_$2.prof.txt Profiling_result[dataset-$1\_Gm=$GMINENTRIES\_GM=$GMAXENTRIES\_GAm=$GAUXMINENTRIES\_GAM=$GAUXMAXENTRIES\_Rm=$RMINENTRIES\_RM=$RMAXENTRIES\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS\_nodes=$5\_threads=$6.txt].txt
cp output.prof.txt Profiling_result[dataset-$1\_Gm=$GMINENTRIES\_GM=$GMAXENTRIES\_GAm=$GAUXMINENTRIES\_GAM=$GAUXMAXENTRIES\_Rm=$RMINENTRIES\_RM=$RMAXENTRIES\_EPS=$Epsilon\_MINPOINTS=$MINPOINTS\_nodes=$5\_threads=$6.txt].txt
#rm -rf $temp

echo "Unification completed on : `date`"
