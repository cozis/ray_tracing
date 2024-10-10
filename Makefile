
ifeq ($(OS),Windows_NT)
    EXT = .exe
	CFLAGS = -O2 -DNDEBUG -I3p/glad/include -I3p/glfw-3.4.bin.WIN64/include -L3p/glfw-3.4.bin.WIN64/lib-mingw-w64
	LDFLAGS = -lglfw3 -lopengl32 -lgdi32
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        EXT =
		CFLAGS = -O2 -DNDEBUG -I3p/glad/include
		LDFLAGS = -lglfw -lm
    endif
    ifeq ($(UNAME_S),Darwin)
        EXT =
		CFLAGS = 
		LDFLAGS = 
    endif
endif

all: ray_trace$(EXT)

ray_trace$(EXT): Makefile $(wildcard src/*.c src/*.h)
	gcc -o $@ src/main.c src/utils.c src/camera.c src/mesh.c src/vector.c src/thread.c src/sync.c src/clock.c src/profile.c 3p/glad/src/glad.c -std=c11 $(CFLAGS) $(LDFLAGS)

clean:
	rm ray_trace ray_trace.exe
