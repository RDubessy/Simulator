all:
	cd src && make all
install:
	make all
	cp bin/simulator /usr/local/bin/
clean:
	cd src && make clean
	cd doc && rm -rf html
	cd bin && rm -rf *
