all:
	cd lib/kent/src/lib && ${MAKE}
	cd code/makeBlocks && ${MAKE}
	cd code && ${MAKE}

clean:
	cd lib/kent/src/lib && ${MAKE} clean
	cd code/makeBlocks && ${MAKE} clean
	cd code && ${MAKE} clean
