include ../Make.inc

OBJCpp = CppDynamicBuffer.o
OBJC = CDynamicBuffer.o
OBJF = FDynamic.o

default:	CppDynamicBuffer.exe CDynamicBuffer.exe FDynamic.exe

CppDynamicBuffer.exe: $(OBJCpp)
	$(CXX) -o $@ $< $(LINK) 

CDynamicBuffer.exe: $(OBJC)
	$(CC) -o $@ $<   $(LINK)

FDynamic.exe: $(OBJF)
	$(F77) -o $@ $<  $(LINK)

cleanall:
	@rm -f *.o *~ *.exe