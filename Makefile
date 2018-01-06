all:

clean:
	find . -name '*.so' -exec rm {} \; &
	find . -name '*.o' -exec rm {} \;

shlib:
	find 'src/algorithms' -name '*.c' | \
		sed -E 's/(.*).c/\1/' | xargs -L1 basename | \
			xargs -L1 -I {} gcc --shared -fPIC -o inst/jl/{}.so src/algorithms/{}.c
