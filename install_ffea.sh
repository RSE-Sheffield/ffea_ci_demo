#!/bin/bash
doinstall=${1:-0}
dotest=${2:-1}
restrict_tests=${3:-1}

if [[ $doinstall -eq 1 ]];
then
	echo "MAKE AND INSTALL"
	rm -r $FFEA_BUILD
	mkdir $FFEA_BUILD
	cd $FFEA_BUILD
	cmake ../ffea -DCMAKE_INSTALL_PREFIX=$FFEA_BUILD -DPYTHON_EXECUTABLE=/home/ryan/Software/anaconda3/envs/py2/bin/python
	make -j2
	make install -j2
elif [[ $doinstall -eq 0 ]];
then
	echo "MAKE"
	cd $FFEA_BUILD
	make -j2
fi

if [[ $dotest -eq 1 ]];
then
	echo "TEST"
	if [[ $restrict_tests -eq 1 ]];
	then
		echo "SPECIFIC"
		# read specific tests from file
		ctest -j2 --verbose -R rod_neighbour_list_construction
	else
		echo "ALL"
		ctest -j2
	fi
fi
echo "DONE"