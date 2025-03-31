#!/bin/bash
# Cleanup script for INCA project

# Remove all module files
echo "Removing module files..."
rm -f *.mod

# Remove object files
echo "Removing object files..."
rm -f obj/*.o

# Remove executable
echo "Removing executable..."
rm -f roda.x

# Remove backup files
echo "Removing backup files..."
rm -f *~

# Remove log files
echo "Removing log files..."
rm -f *.log

echo "Cleanup complete!"
