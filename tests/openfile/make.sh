clang -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes \
-I/opt/baltrad/rave/include -I/opt/baltrad/HLHDF/include \
-L/opt/baltrad/rave/lib     -L/opt/baltrad/HLHDF/lib  \
-L/opt/local/lib/proj47/lib -I/opt/local/lib/proj47/include \
./openfile.c \
-o openfile -lravetoolbox -lhlhdf -lproj
