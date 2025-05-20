#!/bin/bash

# Define paths
SOURCE_DIR="docs/source"
BUILD_DIR="docs/_build/html"
OUTPUT_DIR="docs"

# Clean previous build
echo "Cleaning old HTML files in $OUTPUT_DIR..."
find "$OUTPUT_DIR" -maxdepth 1 -type f -name '*.html' -delete
rm -rf "$OUTPUT_DIR/_static" "$OUTPUT_DIR/_sources" "$OUTPUT_DIR/genindex.html" "$OUTPUT_DIR/objects.inv" "$OUTPUT_DIR/search.html" "$OUTPUT_DIR/searchindex.js"
rm -rf "$BUILD_DIR"
# Build the docs
echo "Building documentation..."
sphinx-build -b html "$SOURCE_DIR" "$BUILD_DIR"

# Copy to docs/ root
echo "Copying built HTML to $OUTPUT_DIR..."
cp -r "$BUILD_DIR"/* "$OUTPUT_DIR"

echo "Done! Your documentation is ready for GitHub Pages."
