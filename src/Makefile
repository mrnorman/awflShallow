
CC := mpiCC
CFLAGS := -O0 -g
INCLUDE := -I${PNETCDF_PATH}/include
LDFLAGS := -L${PNETCDF_PATH}/lib -lpnetcdf -L${LAPACK_PATH}/lib
OMPFLAGS := -mp
ACCFLAGS := -ta=tesla,pinned,cc50,cuda9.0,ptxinfo -Minfo=accel

all:
	${CC} ${INCLUDE} ${CFLAGS} -o sw_cart cfl.cpp timeStepping.cpp haloExchange.cpp transform.cpp finalize.cpp init.cpp io.cpp main.cpp ${LDFLAGS}

clean:
	rm -f *.o sw_cart
