include ../../makeconf

INCPATH = -I../../include

CFLAGS = $(MOREFLAGS) $(WARNLEVEL) -pedantic $(INCPATH)

SOURCES.C    =  jmesh.cpp

OBJECTS.O    =  $(SOURCES.C:.cpp=.o)

../OBJECTS/jmesh.o: $(OBJECTS.O)
	$(LD) -r -o ../OBJECTS/jmesh.o $(OBJECTS.O)

depend:
	makedepend -- $(CFLAGS) -- $(SOURCES.C)
	$(RM) -f Makefile.bak

clean :
	$(RM) -f *.o ../OBJECTS/jmesh.o

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

# DO NOT DELETE

jmesh.o: ../../include/jmesh.h ../../include/j_mesh.h /usr/include/stdio.h
jmesh.o: /usr/include/features.h /usr/include/sys/cdefs.h
jmesh.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
jmesh.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
jmesh.o: /usr/include/bits/typesizes.h /usr/include/libio.h
jmesh.o: /usr/include/_G_config.h /usr/include/wchar.h
jmesh.o: /usr/include/bits/wchar.h /usr/include/gconv.h
jmesh.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
jmesh.o: /usr/include/math.h /usr/include/bits/huge_val.h
jmesh.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
jmesh.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
jmesh.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
jmesh.o: /usr/include/bits/posix2_lim.h ../../include/list.h
jmesh.o: ../../include/point.h ../../include/vertex.h ../../include/edge.h
jmesh.o: ../../include/triangle.h ../../include/matrix.h ../../include/tin.h
jmesh.o: /usr/include/stdlib.h /usr/include/sys/types.h /usr/include/time.h
jmesh.o: /usr/include/endian.h /usr/include/bits/endian.h
jmesh.o: /usr/include/sys/select.h /usr/include/bits/select.h
jmesh.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
jmesh.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
jmesh.o: /usr/include/alloca.h /usr/include/string.h
