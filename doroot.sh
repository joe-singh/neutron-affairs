set -x
g++ -O2 -pipe -Wall -W -Woverloaded-virtual -fPIC -Iinclude  -pthread \
-I $ROOTSYS/include -o $1.o -c $1.cxx


g++ -O2 -m32 $1.o -L$ROOTSYS/lib -lCore -lCint -lRIO -lNet -lHist -lGraf \
-lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore \
-lThread -pthread -lm -ldl -rdynamic  -o $1
