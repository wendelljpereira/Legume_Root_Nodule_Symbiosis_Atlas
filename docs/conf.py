project = "Legume Root Nodule Symbiosis Atlas"
author = "Wendell Pereira and collaborators"
copyright = "2026, Wendell Pereira and collaborators"
release = "0.1.0-preview"

extensions = [
    "myst_parser",
    "sphinx_copybutton",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "attrs_inline",
]

html_theme = "sphinx_rtd_theme"
html_title = "Legume RNS Atlas"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_show_sourcelink = False

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

nitpicky = False
