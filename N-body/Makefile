all: particles particles_serial

particles: particles.cpp 
	mpicxx -o particles particles.cpp 

particles_serial: particles_serial.cpp
	icpc -o particles_serial particles_serial.cpp

clean:
	rm -f particles
