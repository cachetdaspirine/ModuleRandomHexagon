This cpp program is a module that can be easily use into a python program using ctypes module.
The files System_c.cpp/.h is a c version of the cpp program that constitute the core of the shared library.
the shared library is called lib.so, it is the only file needed to use the library onto python.
To compile the program into a library simply run make, if instead you would like to create a simple executable cpp programm (this option can be used to catch errors) then run make EXEC which will create a program called Aight.
A python file called System.py is a simplified interface with the library using only python. because ctypes require c/c++ knowledge, you can just look at System.py.
A jupyter file called TestPythonObject.ipynb shows example of the use of the system.py python object.
