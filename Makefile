exec1   = shuriken_ab12
exec2   = shuriken_ab3
exec3   = shuriken_bigten
exec4   = shuriken_mcfr
cc      = g++
opt     = -O3
cflags  = -std=c++17 -fpermissive $(opt)

main1   = main_ab12.cpp
main2   = main_ab3.cpp
main3   = main_bigten.cpp
main4   = main_mcfr.cpp
objects = $(patsubst %.cpp,%.o,$(filter-out $(main1) $(main2) $(main3) $(main4), $(wildcard *.cpp)))

.PHONY : all clean

all :	$(objects)
	@rm -f $(exec1)
	@$(MAKE) $(exec1)
	@rm -f $(exec2)
	@$(MAKE) $(exec2)
	@rm -f $(exec3)
	@$(MAKE) $(exec3)
	@rm -f $(exec4)
	@$(MAKE) $(exec4)

%.o : %.cpp
	$(cc) $(cflags) -c $<

$(exec1) : $(main1)
	$(cc) $(cflags) $(objects) $< -o $@

$(exec2) : $(main2)
	$(cc) $(cflags) $(objects) $< -o $@

$(exec3) : $(main3)
	$(cc) $(cflags) $(objects) $< -o $@

$(exec4) : $(main4)
	$(cc) $(cflags) $(objects) $< -o $@

clean :
	rm -f $(objects) $(exec)
