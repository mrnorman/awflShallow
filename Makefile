
CC := mpiCC
CFLAGS := -O3
INCLUDE := -I${PNETCDF_PATH}/include
LDFLAGS := -L${PNETCDF_PATH}/lib -lpnetcdf -L${LAPACK_PATH}/lib
OMPFLAGS := -mp
ACCFLAGS := -ta=tesla,pinned,cc50,cuda9.0,ptxinfo -Minfo=accel

all: 
	${CC} ${INCLUDE} ${CFLAGS} -o sw_cart sw_cart.cpp ${LDFLAGS}
