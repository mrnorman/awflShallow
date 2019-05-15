
CC := mpic++
# CFLAGS := -DARRAY_DEBUG -O0 -g -I${PNETCDF_PATH}/include
CFLAGS := -O3 -I${PNETCDF_PATH}/include
LDFLAGS := -L${PNETCDF_PATH}/lib -lpnetcdf

all:
	${CC} ${CFLAGS} driver.cpp -o awfl_sw ${LDFLAGS}

clean:
	rm -f *.gch *.o *.dat awfl_sw
