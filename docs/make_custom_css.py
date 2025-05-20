import os

# Create the _static directory if it doesn't exist
if not os.path.exists('_static'):
    os.makedirs('_static')

# Create a starter custom.css file for Sphinx documentation using the sphinx_rtd_theme
custom_css_content = """
/* _static/custom.css */

/* General body styling */
body {
    font-family: 'Segoe UI', 'Helvetica Neue', sans-serif;
    font-size: 16px;
    line-height: 1.6;
    background-color: #fdfdfd;
    color: #333;
}

/* Content width and padding */
.wy-nav-content {
    max-width: 1000px;
    margin: auto;
    padding: 2em;
}

/* Sidebar styling */
.wy-side-nav-search {
    background-color: #2980b9;
    color: white;
}

.wy-menu-vertical {
    background-color: #f7f7f7;
}

.wy-menu-vertical a {
    font-weight: 500;
    color: #2980b9;
}

.wy-menu-vertical a:hover {
    color: #1a5276;
}

/* Link styling */
a {
    color: #2980b9;
}

a:hover {
    color: #1a5276;
}

/* Header styling */
h1, h2, h3, h4, h5, h6 {
    color: #2c3e50;
    font-weight: 600;
}

/* Code block styling */
code {
    background-color: #f8f8f8;
    padding: 2px 4px;
    border-radius: 4px;
}

/* Table styling */
table {
    width: 100%;
    border-collapse: collapse;
}

table, th, td {
    border: 1px solid #ddd;
}

th, td {
    padding: 8px;
    text-align: left;
}

th {
    background-color: #f2f2f2;
}
"""

# Write the custom CSS content to a file
with open('_static/custom.css', 'w') as f:
    f.write(custom_css_content)

print("The custom.css file has been created successfully.")

