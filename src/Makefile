FC = gfortran
FFLAGS = -O2 -Wall -std=f2008

# Directories
SRCDIR = .
BUILDDIR = build
BINDIR = bin

# Create directories
$(shell mkdir -p $(BUILDDIR) $(BINDIR))

# Source files
SOURCES = kerr_geometry.f90 integrators.f90 trajectory.f90 kerr_simulator.f90
OBJECTS = $(addprefix $(BUILDDIR)/, $(SOURCES:.f90=.o))

# Target executable
TARGET = $(BINDIR)/kerr_sim

# Main target
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)

# Object file rules
$(BUILDDIR)/kerr_geometry.o: $(SRCDIR)/kerr_geometry.f90
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILDDIR)

$(BUILDDIR)/integrators.o: $(SRCDIR)/integrators.f90
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILDDIR)

$(BUILDDIR)/trajectory.o: $(SRCDIR)/trajectory.f90 $(BUILDDIR)/kerr_geometry.o $(BUILDDIR)/integrators.o
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILDDIR)

$(BUILDDIR)/kerr_simulator.o: $(SRCDIR)/kerr_simulator.f90 $(BUILDDIR)/kerr_geometry.o $(BUILDDIR)/trajectory.o
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILDDIR)

# Clean rule
clean:
	rm -rf $(BUILDDIR) $(BINDIR) trajectory.dat

# Run target
run: $(TARGET)
	./$(TARGET)

# Default target
all: $(TARGET)

.PHONY: clean run all