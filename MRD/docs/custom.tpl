\documentclass[paper=letter, fontsize=12pt]{article}
\usepackage[magyar]{babel} 
\usepackage{amsmath,amsfonts,amsthm} % Math packages||nem tudom kell-e
\usepackage[utf8]{inputenc}
\usepackage{blindtext}
\usepackage{graphicx}
\usepackage{amsmath}
\usepackage{physics}
\usepackage{mathtools}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{verbatim}
\usepackage[sc]{mathpazo} % Use the Palatino font
\usepackage[T1]{fontenc} % Use 8-bit encoding that has 256 glyphs
\linespread{1.05} % Line spacing - Palatino needs more space between lines
\usepackage{microtype} % Slightly tweak font spacing for aesthetics
\usepackage[left=24mm,hmarginratio=1:1,top=24mm,bottom=24mm,columnsep=10pt]{geometry} % Document margins
\usepackage{url}
\usepackage{xfrac}
\usepackage{nicefrac}
\usepackage{booktabs} % Horizontal rules in tables
\usepackage{float} % Required for tables and figures in the multi-column environment - they need to be placed in specific locations with the [H] 
\usepackage{fancyhdr} % Headers and footers
\pagestyle{fancy} % All pages have headers and footers
\fancyhead{} % Blank out the default header
\fancyfoot{} % Blank out the default footer
\fancyfoot[R]{\thepage} % Custom footer text
\usepackage[nottoc]{tocbibind}
\usepackage{multirow}
\usepackage{pdfpages}

\usepackage{lmodern}
\usepackage{bm}
\usepackage{microtype}
\usepackage{hyperref}
\setlength{\parindent}{0pt}
\setlength{\parskip}{1.2ex}

\hypersetup
       {   pdfauthor = { {{{:author}}} },
           pdftitle={ {{{:title}}} },
           colorlinks=TRUE,
           linkcolor=black,
           citecolor=blue,
           urlcolor=blue
       }

{{#:title}}
\title{ {{{ :title }}} }
{{/:title}}

{{#:author}}
\author{ {{{ :author }}} }
{{/:author}}

{{#:date}}
\date{ {{{ :date }}} }
{{/:date}}

{{ :highlight }}

\begin{document}

{{#:title}}\maketitle{{/:title}}

{{{ :body }}}

\end{document}