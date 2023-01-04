CC = gcc
CFLAG = -Wall -fopenmp
link = -lfftw3_threads -lfftw3 -lm -lpthread -lhdf5
look = -L/usr/nld/hdf5-1.8.9/lib -I /usr/nld/hdf5-1.8.9/include/
all = main.o time_march.o init_flow.o find_vel.o storage.o

at_desktop : $(all)
	$(CC) -o $@ $(all) $(look) $(link) $(CFLAG)
	
main.o : main.c
	$(CC) -c -o $@ main.c $(look) $(CFLAG)
	
time_march.o : time_march.c
	$(CC) -c -o $@ time_march.c $(look) $(CFLAG)
	
init_flow.o : init_flow.c
	$(CC) -c -o $@ init_flow.c $(look) $(CFLAG)
	
find_vel.o : find_vel.c
	$(CC) -c -o $@ find_vel.c $(look) $(CFLAG)
	
storage.o: storage.c
	$(CC) -c -o $@ storage.c $(look) $(CFLAG)


clean:
	rm *.o

clean_data:
	rm field/* spectra/*
