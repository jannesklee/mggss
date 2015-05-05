TARGET = fdgss.out
CC = gcc
CFLAGS = -std=c99 -pedantic -Wall -Wextra -Wdeclaration-after-statement -O2
DEBUG = -g
SRCDIR = src
OBJDIR = obj
#BINDIR = bin

SOURCES = $(wildcard $(SRCDIR)/**/*.c $(SRCDIR)/*.c)
DEPS = $(wildcard $(SRCDIR)/**/*.h $(SRCDIR)/*.h)
OBJECTS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o,$(SOURCES))

# not used yet
#TEST_SRC = $(wildcard tests/*_tests.cpp)
#TESTS = $(patsubst %.cpp,%,$(TEST_SRC))

# compiling
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	LC_MESSAGES=C $(CC) -o $@ $< -c $(CFLAGS)

# linking
$(TARGET): $(OBJECTS)
	LC_MESSAGES=C $(CC) -o $@ $^ $(CFLAGS)

all: main.out

# cleaning procedure
.PHONY: clean
clean:
	rm $(OBJECTS) $(TARGET)

