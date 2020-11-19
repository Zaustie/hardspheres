all: rule1


rule1:
		g++ main.cpp -o spheres -lglut -lGLU -lGL -lSDL2 -O3
clean:
	    rm spheres

debug:
		g++ main.cpp -o spheres -g -lglut -lGLU -lGL -lSDL2 -O3
