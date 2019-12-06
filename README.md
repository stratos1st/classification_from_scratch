# project2
HOW to RUN unit_testing_main:

FIRST DOWNLOAD/INSTALL THE LIBRARY
-----------------------------
sudo apt-get install libgtest-dev
sudo apt-get install cmake
cd /usr/src/gtest
sudo cmake CMakeLists.txt
sudo make
sudo cp *.a /usr/lib

---------------------------
cmake CMakeLists.txt
make -f Makefile
./executeTests
-----------------------------
