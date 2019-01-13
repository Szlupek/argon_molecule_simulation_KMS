CC		=g++
#kompilator uzywany do kompilacji klas C++
LD		=g++
#kompilator uzywany do kompilacji finalnej aplikacji
COPTS		=
#flagi uzywane przy kompilacji klas C++
LDOPTS		=
#flagi uzywane przy kompilacji finalnej aplikacji

# *** C++ files ***
SOURCES		=argon.cpp argon_functions.cpp
#pliki zrodlowe
OBJECTS		= $(SOURCES:.cpp=.o)
#pliki utworzone przez kompilator, nazwy generowane na podstawie nazw plikow zrodlowych [nie trzeba zmieniac]

# *** Dictionary classes ***
HEADERS		=plik.h
#pliki naglowkowe
EXECUTABLE	=program
#nazwa pliku wykonywalnego

all:		$(EXECUTABLE)
#to zrobi komenda 'make' [nie trzeba zmieniac]

$(EXECUTABLE): $(OBJECTS)
		$(LD) -o $@ $^ $(LDOPTS)
#to jest mechanizm kompilujacy finalna aplikacje [nie trzeba zmieniac]

# *** C++ files ***
.cpp.o:		$(CC) -o $@ $^ -c $(COPTS)
#to jest mechanizm kompilujacy klasy C++ [nie trzeba zmieniac]

clean:;		@rm -f $(OBJECTS) $(EXECUTABLE) *.o *.d
#to zrobi komenda 'make clean' [nie trzeba zmieniac]
