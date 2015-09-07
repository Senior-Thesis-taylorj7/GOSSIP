default: dpAlign.c dpAlign_wrap.c dpAlign.i
	swig -python dpAlign.i
	gcc -c -fPIC dpAlign.c dpAlign_wrap.c -I/usr/include/python2.7
	ld -shared dpAlign.o dpAlign_wrap.o -o _dpAlign.so

install:
	cp ./dpAlign.py ..
	cp ./_dpAlign.so ..

clean:
	rm -rf *.o _dpAlign.so dpAlign.py dpAlign.pyc
