#!/bin/bash

# Example script to test pager functionality in Docker container
# This script demonstrates basic pager module functionality

echo "=== Pager Docker Container Test ==="
echo ""

# Test Python import
echo "Testing Python environment..."
python3 -c "
import sys
sys.path.append('/opt/pager')
import pager
import pandas as pd
import os

print('✓ Python 3 available')
print('✓ pager module imported successfully')
print('✓ pandas available:', pd.__version__)
print('✓ Working directory:', os.getcwd())
"

echo ""
echo "=== Environment Variables ==="
echo "PAGEDIR: $PAGEDIR"
echo "TEISERDIR: $TEISERDIR"

echo ""
echo "=== Available Scripts ==="
ls -la /opt/pager/*.sh

echo ""
echo "=== Test Data ==="
if [ -f "/opt/pager/test/input.txt" ]; then
    echo "Test input data found:"
    head -n 5 /opt/pager/test/input.txt
    echo "... (showing first 5 lines)"
else
    echo "No test data found in /opt/pager/test/"
fi

echo ""
echo "=== Docker Container Ready ==="
echo "To get started:"
echo "1. Mount your data: docker run -v /path/to/data:/data pager"
echo "2. Mount iPAGE/TEISER tools for full functionality"
echo "3. Use 'python3 -c \"import pager; help(pager)\"' to explore the module"